#!/usr/bin/env python3
"""Aliasing-Test fuer Wellenlaengen-Sweeps (EXFO JSON).

Ein Peak im Reflektogramm kann aus drei Gruenden falsch sein:
  (a) er ist gefaltet   -> Laenge falsch
  (b) er ist kein Peak, sondern ein Berg -> Aufloesung nicht erreicht
  (c) beides

Dieses Skript prueft alle drei und sagt ausdruecklich, wenn es die Frage
NICHT entscheiden kann.

  python3 alias_check.py datei.json                 # eine Messung bewerten
  python3 alias_check.py datei_A.json datei_B.json  # Paartest (der harte Test)

Der Paartest ist das einzige narrensichere Verfahren: zwei Sweeps derselben
Optik mit VERSCHIEDENEM delta-lambda. Bleibt der Peak an derselben Stelle
(in cm), ist die Messung gueltig. Wandert er, ist sie gefaltet.
"""
from __future__ import annotations

import sys
import json
from pathlib import Path

import numpy as np

N_GROUP = 1.468
C0 = 299_792_458.0
K_BANDS = 8


# ------------------------------------------------------------------ laden
def load(path: Path):
    raw = json.load(path.open("r", encoding="utf-8"))
    block = raw["data"][0]
    info = (block.get("Device Description") or {}).get("info", "")
    arr = {k: np.asarray(v, float) for k, v in block.items()
           if isinstance(v, list) and v and isinstance(v[0], (int, float))}
    wl_key = next(k for k in arr if k.lower().startswith("wavelength"))
    wl = arr.pop(wl_key)

    step = np.median(np.diff(wl))
    keep = np.ones(wl.size, bool)
    for i in np.where(np.abs(np.diff(wl) - step) > 50 * abs(step))[0]:
        keep[i + 1] = False
    wl, arr = wl[keep], {k: v[keep] for k, v in arr.items()}

    ports = [k for k in arr if k.lower().startswith("ch")]
    if len(ports) >= 2 and np.corrcoef(arr[ports[0]], arr[ports[1]])[0, 1] < -0.3:
        a, b = arr[ports[0]], arr[ports[1]]
        sig, chan = (a - b) / (a + b), f"balanciert {ports[0]}/{ports[1]}"
    else:
        sig, chan = arr[ports[0]], ports[0]
    return wl, sig, chan, info


def spectrum(sig, lam):
    """Reflektogramm: Amplitude gegen Reflektorabstand [cm]."""
    w = 501 | 1
    env = np.convolve(np.pad(sig, (w // 2, w // 2), mode="edge"),
                      np.ones(w) / w, mode="valid")[: sig.size]
    v = (sig - env) * np.hanning(sig.size)
    k = 1.0 / lam
    o = np.argsort(k)
    ku = np.linspace(k[o][0], k[o][-1], k.size)
    vu = np.interp(ku, k[o], v[o])
    S = np.abs(np.fft.rfft(vu))
    f = np.fft.rfftfreq(vu.size, d=ku[1] - ku[0])
    return f / N_GROUP * 1e-9 * 50, S          # cm, Amplitude


def band_freq(sig, d_lam, f_min):
    w = 301 | 1
    env = np.convolve(np.pad(sig, (w // 2, w // 2), mode="edge"),
                      np.ones(w) / w, mode="valid")[: sig.size]
    v = (sig - env) * np.hanning(sig.size)
    S = np.abs(np.fft.rfft(v))
    f = np.fft.rfftfreq(v.size, d=d_lam)
    m = f > f_min
    return float((f[m] * S[m] ** 2).sum() / (S[m] ** 2).sum())


# ---------------------------------------------------------------- analyse
def analyse(path: Path, verbose=True):
    wl, sig, chan, info = load(path)
    lam = np.linspace(wl[0], wl[-1], wl.size)
    d_lam = lam[1] - lam[0]
    lam0, span = float(lam.mean()), float(lam[-1] - lam[0])

    B = C0 * span * 1e-9 / (lam0 * 1e-9) ** 2
    d_res_cm = C0 / (2 * N_GROUP * B) * 100
    d_max = lam0 ** 2 / (4 * N_GROUP * d_lam) * 1e-9 * 100

    d, S = spectrum(sig, lam)
    P = S ** 2
    m = d > 0.3
    dd, PP = d[m], P[m]

    # Rauschboden aus dem obersten Viertel, das der Berg meist nicht erreicht
    floor = np.median(PP[dd > 0.85 * d_max])
    sig_only = np.maximum(PP - floor, 0)
    pk = float(dd[np.argmax(sig_only)])

    # Breite bei halber Leistung (nur der Berg, Rauschboden abgezogen)
    thr = sig_only.max() / 2
    above = dd[sig_only > thr]
    fwhm = float(above.max() - above.min()) if above.size else float("nan")
    broadening = fwhm / d_res_cm

    # Drift der Fringe-Frequenz
    f_app = band_freq(sig, d_lam, f_min=0.05 / d_lam)
    n = sig.size
    L, F = [], []
    for j in range(K_BANDS):
        s = slice(j * n // K_BANDS, (j + 1) * n // K_BANDS)
        L.append(float(lam[s].mean()))
        F.append(band_freq(sig[s], d_lam, f_min=0.3 * f_app))
    L, F = np.array(L), np.array(F)
    x = 1.0 / L ** 2
    slope = float(np.polyfit(x, F, 1)[0])
    r = float(np.corrcoef(x, F)[0, 1])

    drift_meas = abs(F[-1] - F[0]) / F.mean()
    drift_theo = 1.0 - (lam[0] / lam[-1]) ** 2      # ~6.25 % bei 1520..1570 nm
    drift_ratio = drift_meas / drift_theo
    opd_chirp_cm = abs(slope) / N_GROUP * 1e-9 * 100

    res = dict(name=path.name, info=info, chan=chan, n=lam.size, d_lam_pm=d_lam * 1e3,
               span=span, B_THz=B / 1e12, d_res_cm=d_res_cm, d_max=d_max,
               peak=pk, fwhm=fwhm, broadening=broadening,
               slope=slope, r=r, drift_ratio=drift_ratio,
               opd_chirp_cm=opd_chirp_cm, L=L, F=F, floor_snr=sig_only.max() / floor)

    if verbose:
        report(res)
    return res


def report(res):
    print(f"Datei          : {res['name']}")
    print(f"Info-Feld      : {res['info']}")
    print(f"Kanal          : {res['chan']}")
    print(f"Punkte         : {res['n']}     delta-lambda = {res['d_lam_pm']:.4f} pm")
    print(f"Bandbreite     : {res['span']:.3f} nm = {res['B_THz']:.3f} THz")
    print(f"Aufloesung     : {res['d_res_cm']*1e4:.1f} um")
    print(f"Nyquist-Grenze : {res['d_max']:.2f} cm")
    print()
    print(f"Maximum bei    : {res['peak']:.2f} cm")
    print(f"Breite (FWHM)  : {res['fwhm']:.2f} cm  =  {res['broadening']:.0f} x Aufloesung")
    print(f"Berg / Rauschen: {res['floor_snr']:.1f} x")
    print()

    # --- Test 1: ist es ueberhaupt ein Peak? ---
    if res["broadening"] > 20:
        print("[X] KEIN PEAK. Der Berg ist "
              f"{res['broadening']:.0f}-fach breiter als die Aufloesung.")
        print("    Eine Laengenangabe daraus ist bedeutungslos, egal ob gefaltet oder nicht.")
    elif res["broadening"] > 5:
        print(f"[!] Peak deutlich verbreitert ({res['broadening']:.0f}x). Laenge grob brauchbar,")
        print("    die 16 um Aufloesung erreichst du so nicht.")
    else:
        print(f"[OK] Sauberer Peak ({res['broadening']:.0f}x Aufloesung).")
    print()

    # --- Test 2: Drift-Auswertung, nur wenn zulaessig ---
    print("Drift der Fringe-Frequenz ueber den Sweep:")
    for l, f in zip(res["L"], res["F"]):
        print(f"   {l:8.2f} nm : {f:9.2f} cycles/nm")
    print(f"   Fit gegen 1/lambda^2: Steigung {res['slope']:+.4g}   (r = {res['r']:+.4f})")
    print(f"   Drift-Betrag        : {res['drift_ratio']:.2f} x der Theorie "
          f"(ein Einzelziel muss 1.00 geben)")
    print()

    if abs(res["r"]) < 0.8 or not (0.6 < res["drift_ratio"] < 1.7):
        print("[?] DRIFT-AUSWERTUNG NICHT ZULAESSIG.")
        print("    Der Drift-Betrag passt nicht zu einem einzelnen Reflektor -- vermutlich")
        print("    mehrere ueberlagerte (gefaltete) Komponenten. Ich kann aus DIESER Datei")
        print("    allein nicht sagen, ob und wie oft gefaltet ist.")
        print("    -> Paartest fahren (siehe unten).")
        return

    parity = "GERADE (0, 2, 4, ...)" if res["slope"] > 0 else "UNGERADE (1, 3, ...)"
    print(f"Faltungs-Ordnung ist {parity}")
    print(f"OPD aus dem Drift-Betrag : {res['opd_chirp_cm']:.1f} cm  ->  "
          f"Abstand {res['opd_chirp_cm']/2:.1f} cm")
    if res["opd_chirp_cm"] / 2 > res["d_max"]:
        print(f"    {res['opd_chirp_cm']/2:.1f} cm > Nyquist {res['d_max']:.2f} cm  "
              f"==>  GEFALTET.")
    elif res["slope"] > 0:
        print("    innerhalb Nyquist und gerade Ordnung  ==>  Ordnung 0, Messung gueltig.")
    else:
        print("    ungerade Ordnung  ==>  GEFALTET.")


def pair(a: Path, b: Path):
    ra, rb = analyse(a, verbose=False), analyse(b, verbose=False)
    print("=" * 68)
    print("PAARTEST")
    print("=" * 68)
    for r in (ra, rb):
        print(f"  {r['name'][:44]:44s} dl={r['d_lam_pm']:6.3f} pm  "
              f"Nyquist {r['d_max']:6.2f} cm  Maximum {r['peak']:6.2f} cm")
    tol = 3 * max(ra["d_res_cm"], rb["d_res_cm"]) + 0.05 * max(ra["fwhm"], rb["fwhm"])
    diff = abs(ra["peak"] - rb["peak"])
    print(f"\n  Verschiebung: {diff:.2f} cm   (Toleranz {tol:.2f} cm)")
    if abs(ra["d_lam_pm"] - rb["d_lam_pm"]) < 1e-3:
        print("\n  [?] Beide Messungen haben dasselbe delta-lambda. Der Paartest braucht")
        print("      zwei VERSCHIEDENE Schrittweiten derselben Optik.")
        return
    if diff <= tol:
        print("\n  [OK] Peak bleibt stehen  ==>  NICHT gefaltet. Die Laenge ist echt.")
    else:
        print("\n  [X] Peak wandert  ==>  mindestens eine der beiden ist GEFALTET.")
        # gemeinsame Kandidaten
        def cands(r, kmax=6):
            P = 2 * r["d_max"]
            out = set()
            for k in range(kmax):
                out.add(round(k * P + r["peak"], 2))
                if k >= 1:
                    out.add(round(k * P - r["peak"], 2))
            return sorted(x for x in out if x > 0)
        ca, cb = cands(ra), cands(rb)
        both = [(x, y) for x in ca for y in cb if abs(x - y) < 2.0]
        print("\n      Mit beiden Messungen vereinbare wahre Abstaende:")
        for x, y in both[:12]:
            print(f"        {(x+y)/2:8.2f} cm")
        print("\n      Warnung: ist ein delta-lambda ein ganzzahliges Vielfaches des")
        print("      anderen, bleiben mehrere Kandidaten uebrig -- das ist mathematisch")
        print("      unvermeidbar. Nimm ein Paar wie 1 pm / 1.5 pm, oder geh feiner.")


if __name__ == "__main__":
    args = [Path(a) for a in sys.argv[1:]]
    if not args:
        print(__doc__)
        sys.exit(1)
    if len(args) == 1:
        analyse(args[0])
    else:
        pair(args[0], args[1])
