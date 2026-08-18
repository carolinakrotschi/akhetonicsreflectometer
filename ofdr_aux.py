#!/usr/bin/env python3
"""
OFDR-Verarbeitung fuer den FREILAUF-Modus mit Aux-MZI (4 Kanaele).

Das ist die Weiterentwicklung von ofdr_process.py fuer den neuen Aufbau.
Unterschied in einem Satz: die Frequenzachse kommt nicht mehr aus der
Wellenlaengentabelle des EXFO (die es im Single-Trigger-Modus nicht mehr
gibt), sondern aus der Phase des Aux-MZI.

    Ch1 / Ch2 : Mess-Interferometer (komplementaer)
    Ch3 / Ch4 : Aux-MZI = Lineal   (komplementaer)

Was aus ofdr_process.py UNVERAENDERT uebernommen ist:
    balancierte Subtraktion, Fensterung, FFT, Peak-Liste, Breitenpruefung
Was WEGFAELLT:
    to_uniform_nu()   -- keine Wellenlaengentabelle mehr
    phase_correct()   -- keine Selbstreferenzierung, also auch kein
                         diagnostic/cosmetic-Konflikt mehr
Was NEU ist:
    resample_on_aux() -- eine Funktion, ~25 Zeilen

WICHTIG: die Aux-Phase wird NICHT geglaettet. Sie muss die schnelle
Abstimmwelligkeit mittragen, sonst korrigiert sie sie nicht weg.

Kalibrierung von tau_aux (einmal, im ALTEN Triggermodus mit
Wellenlaengentabelle):
    tau_aux = gesamte entwickelte Aux-Phase / (2*pi * Frequenzspanne)
Uebergib das Ergebnis mit --tau-aux-ns. Ohne Angabe wird es aus --dl
geschaetzt (tau_aux = n_g * dL / c), was fuer Positionen auf ~1 % genau ist.

Beispiele:
    python ofdr_aux.py testdata.npz --dl 4
    python ofdr_aux.py scan.json --tau-aux-ns 19.587 --zmax 2.5
"""

import argparse
import json
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks
from scipy.signal.windows import kaiser, hann, blackmanharris

C = 299_792_458.0
NG = 1.468


# ---------------------------------------------------------------- laden
def load(path):
    """-> (ch1, ch2, ch3, ch4, meta). Akzeptiert .npz und .json."""
    if path.endswith(".npz"):
        d = np.load(path)
        meta = {k: d[k].item() for k in
                ("step_us", "lam0_nm", "speed_nms") if k in d}
        if "truth_z" in d:
            meta["truth_z"] = d["truth_z"]
            meta["truth_db"] = d["truth_db"]
        return (d["ch1"], d["ch2"], d["ch3"], d["ch4"], meta)

    with open(path) as f:
        doc = json.load(f)
    e = doc["data"][0]
    want = ["Ch1 [mW]", "Ch2 [mW]", "Ch3 [mW]", "Ch4 [mW]"]
    missing = [k for k in want if k not in e]
    if missing:
        sys.exit(f"Kanaele fehlen in der Datei: {missing}\n"
                 f"vorhandene Schluessel: {list(e)}")
    ch = [np.asarray(e[k], float) for k in want]
    meta = {k: v for k, v in (e.get("Device Description") or {}).items()}
    return (*ch, meta)


# ------------------------------------------------- balancierte Subtraktion
def balanced(a, b, nseg=32):
    """P = a - g*b, g abschnittsweise per Median (identisch zu
    ofdr_process.balanced_subtract)."""
    n = len(a)
    seg = max(n // nseg, 256)
    g = np.empty(n)
    for i in range(0, n, seg):
        s = slice(i, min(i + seg, n))
        m = np.median(b[s])
        g[s] = np.median(a[s]) / m if m > 0 else 1.0
    g = uniform_filter1d(g, seg)
    return a - g * b, g


def analytic(x):
    n = len(x)
    X = np.fft.fft(x)
    X[n // 2 + 1:] = 0.0
    X[1:n // 2] *= 2.0
    return np.fft.ifft(X)


# ------------------------------------------------------ der neue Kern
def resample_on_aux(meas, aux, tau_aux, trim=0.01):
    """Messsignal auf gleichmaessige Aux-Phase umtasten.

    Kern der Sache: die Aux-Phase ist phi(t) = 2*pi*tau_aux*nu(t) + const.
    Gleichmaessige Schritte in phi sind also gleichmaessige Schritte in nu --
    egal wie ungleichmaessig der Laser tatsaechlich gefahren ist. Damit ist
    die Achse fertig und die FFT erlaubt.

    Rueckgabe: (y, dnu, span_nu, diag)
    """
    n = len(aux)
    an = analytic(aux)
    phi = np.unwrap(np.angle(an))
    if phi[-1] < phi[0]:                 # nu faellt, wenn lambda steigt
        phi = -phi

    # Raender wegschneiden: dort laeuft der Laser noch an, und die
    # Hilbert-Transformation hat Kantenartefakte.
    k = max(1, int(trim * n))
    sl = slice(k, n - k)
    phi = phi[sl]
    meas = meas[sl]
    amp = np.abs(an)[sl]

    # Monotonie erzwingen (Rauschen kann winzige Rueckschritte machen).
    # Die Stufe muss gross gegen die Rechengenauigkeit von phi sein, sonst
    # bleiben Duplikate stehen -- daher relativ zum mittleren Schritt.
    eps = 1e-6 * (phi[-1] - phi[0]) / len(phi)
    phi = np.maximum.accumulate(phi) + np.arange(len(phi)) * eps

    m = len(phi)
    phi_u = np.linspace(phi[0], phi[-1], m)
    y = PchipInterpolator(phi, meas)(phi_u)

    dphi = phi_u[1] - phi_u[0]
    dnu = dphi / (2 * np.pi * tau_aux)
    span_nu = (phi[-1] - phi[0]) / (2 * np.pi * tau_aux)
    fringes = (phi[-1] - phi[0]) / (2 * np.pi)
    diag = dict(fringes=fringes, pts_per_fringe=m / fringes,
                amp_min=amp.min() / amp.mean(), n_used=m)
    return y, dnu, span_nu, diag


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scan")
    p.add_argument("--tau-aux-ns", type=float, default=None,
                   help="kalibrierte Aux-Verzoegerung in ns (bevorzugt)")
    p.add_argument("--dl", type=float, default=None,
                   help="Aux-Armdifferenz in m (Notloesung, wenn tau_aux "
                        "noch nicht kalibriert ist)")
    p.add_argument("--window", default="kaiser",
                   choices=["hann", "blackmanharris", "kaiser"])
    p.add_argument("--kaiser-beta", type=float, default=12.0)
    p.add_argument("--zmax", type=float, default=None)
    p.add_argument("--peak-floor-db", type=float, default=-45.0)
    p.add_argument("--trim", type=float, default=0.01,
                   help="Anteil, der an jedem Rand verworfen wird")
    p.add_argument("--out", default=None)
    a = p.parse_args()

    if a.tau_aux_ns is not None:
        tau_aux = a.tau_aux_ns * 1e-9
        src = "kalibriert"
    elif a.dl is not None:
        tau_aux = NG * a.dl / C
        src = f"geschaetzt aus dL = {a.dl} m"
    else:
        sys.exit("--tau-aux-ns oder --dl angeben")

    ch1, ch2, ch3, ch4, meta = load(a.scan)
    n = len(ch1)
    print(f"{a.scan}: {n:,} Punkte x 4 Kanaele")
    print(f"tau_aux = {tau_aux*1e9:.4f} ns  ({src})"
          f"  -> Aux erscheint bei z = {C*tau_aux/(2*NG):.3f} m")

    meas, gm = balanced(ch1, ch2)
    aux, ga = balanced(ch3, ch4)
    print(f"balancierte Subtraktion: g_mess {gm.min():.3f}..{gm.max():.3f}, "
          f"g_aux {ga.min():.3f}..{ga.max():.3f}")

    y, dnu, span_nu, diag = resample_on_aux(meas, aux, tau_aux, a.trim)
    m = len(y)
    print(f"Aux: {diag['fringes']:,.0f} Fringes, "
          f"{diag['pts_per_fringe']:.1f} Punkte/Fringe, "
          f"Amplitudenminimum {diag['amp_min']*100:.0f} % vom Mittel")
    if diag["pts_per_fringe"] < 4:
        print("  WARNUNG: unter 4 Punkte pro Aux-Fringe -- langsamer sweepen "
              "oder kuerzeres Aux.")
    if diag["amp_min"] < 0.2:
        print("  WARNUNG: Aux-Kontrast bricht irgendwo ein (Polarisations-"
              "fading?). Phase dort unzuverlaessig.")

    # Restbasislinie entfernen (wie in ofdr_process.py)
    t = np.linspace(-1, 1, m)
    y = y - np.polyval(np.polyfit(t, y, 5), t)

    dz_bin = C / (2 * NG * span_nu)
    z_nyq = C / (4 * NG * dnu)
    lam_mid = 1565e-9
    print(f"Span {span_nu/1e12:.3f} THz "
          f"(~{span_nu*lam_mid**2/C*1e9:.1f} nm), Punktabstand "
          f"{dnu/1e6:.2f} MHz")
    print(f"  AUFLOESUNG {dz_bin*1e6:.2f} um   REICHWEITE {z_nyq:.3f} m   "
          f"Zellen {z_nyq/dz_bin:,.0f}")
    print(f"  !! kein Reflektor jenseits {z_nyq:.2f} m, sonst aliasiert er")

    win = {"hann": hann(m), "blackmanharris": blackmanharris(m),
           "kaiser": kaiser(m, a.kaiser_beta)}[a.window]
    R = np.abs(np.fft.rfft(y * win))
    z = np.arange(len(R)) * C / (2 * NG * dnu * m)
    db = 20 * np.log10(R / R.max() + 1e-15)

    i = int(np.argmax(R))
    half = R[i] / np.sqrt(2)
    l = r = i
    while l > 0 and R[l] > half:
        l -= 1
    while r < len(R) - 1 and R[r] > half:
        r += 1
    width = (r - l) * (z[1] - z[0])
    wlim = {"hann": 1.6, "blackmanharris": 2.7, "kaiser": 2.6}[a.window]
    print(f"\nHauptpeak {z[i]*1000:.4f} mm, -3 dB Breite {width*1e6:.1f} um "
          f"(Fenstergrenze ~{wlim*dz_bin*1e6:.1f} um)")
    if width > 2 * wlim * dz_bin:
        print("  WARNUNG: Hauptpeak > 2x Fenstergrenze -- Frequenzachse "
              "verdaechtig. Erst untersuchen, dann Positionen glauben.")

    zmax = a.zmax if a.zmax else z[-1]
    pk, _ = find_peaks(db, height=a.peak_floor_db,
                       distance=max(3, int(200e-6 / (z[1] - z[0]))))
    print(f"\nPeaks ueber {a.peak_floor_db:.0f} dB:")
    for j in pk:
        if z[j] <= zmax:
            h = z[j] / z[i] if z[i] > 0 else 0
            tag = ("   <- Harmonische, kein Reflektor"
                   if 1.5 < h < 6 and abs(h - round(h)) < 0.03 else "")
            print(f"   {z[j]*1000:10.4f} mm  {db[j]:6.1f} dB{tag}")

    # Selbsttest gegen bekannte Wahrheit (nur bei synthetischen Daten)
    if "truth_z" in meta:
        print("\n--- Vergleich mit bekannter Wahrheit ---")
        worst = 0.0
        for zt, dbt in zip(meta["truth_z"], meta["truth_db"]):
            w = (z > zt - 20 * dz_bin) & (z < zt + 20 * dz_bin)
            if not w.any():
                print(f"   z_wahr {zt:.4f} m  -- ausserhalb des Bereichs")
                continue
            jj = int(np.argmax(np.where(w, R, 0)))
            err = (z[jj] - zt) * 1e6
            worst = max(worst, abs(err))
            print(f"   z_wahr {zt:7.4f} m ({dbt:5.1f} dB) -> gefunden "
                  f"{z[jj]:7.4f} m, Fehler {err:+7.1f} um "
                  f"({err/(dz_bin*1e6):+.2f} Zellen), {db[jj]:6.1f} dB")
        print(f"   groesster Fehler: {worst:.1f} um "
              f"(eine Zelle = {dz_bin*1e6:.1f} um)")
        print("   -> BESTANDEN" if worst < 2 * dz_bin * 1e6
              else "   -> DURCHGEFALLEN, Pipeline pruefen")

    prefix = a.out or a.scan.rsplit(".", 1)[0]
    keep = z <= zmax
    np.savetxt(f"{prefix}_reflectogram.csv",
               np.column_stack([z[keep], db[keep]]), delimiter=",",
               header="distance_m,amplitude_dB", comments="")
    print(f"\ngeschrieben: {prefix}_reflectogram.csv")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.plot(z[keep], db[keep], lw=0.6)
        ax.set_xlabel("Entfernung (m, einfacher Weg / Reflexionskonvention)")
        ax.set_ylabel("Amplitude (dB rel. Maximum)")
        ax.set_title(f"{a.scan} | aux-referenziert, {a.window}, "
                     f"dz {dz_bin*1e6:.1f} um, Nyquist {z_nyq:.2f} m")
        ax.grid(alpha=0.3)
        ax.set_ylim(max(-110, db[keep].min() - 5), 5)
        fig.tight_layout()
        fig.savefig(f"{prefix}_reflectogram.png", dpi=150)
        print(f"geschrieben: {prefix}_reflectogram.png")
    except Exception as e:
        print(f"(Plot uebersprungen: {e})")


if __name__ == "__main__":
    main()
