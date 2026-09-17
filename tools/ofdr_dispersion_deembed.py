#!/usr/bin/env python3
"""Dispersions-De-Embedding nach Zhao et al., IEEE PTL 29(16), 1379 (2017).

  D. Zhao, D. Pustakhod, K. Williams, X. Leijtens, "High Resolution Optical
  Frequency Domain Reflectometry for Analyzing Intra-Chip Reflections",
  doi:10.1109/LPT.2017.2723242

Das Paper misst Reflexionen INNERHALB eines InP-Chips und holt sich die dafuer
noetige Aufloesung, indem es die Gruppenlaufzeitdispersion des Wellenleiters
aus den Daten selbst bestimmt und herausrechnet. Vier Schritte (Abschnitt II):

  1. FFT nicht ueber das ganze Spektrum, sondern ueber ein gleitendes
     Rechteckfenster [lam - w/2, lam + w/2]. Das gibt ein Reflektogramm JE
     WELLENLAENGE und damit die Gruppenlaufzeit tau_f, tau_b von Vorder- und
     Rueckfacette als Funktion von lam.
  2. Gruppenindex nach Gl. (3):   ng(lam) = c * (tau_b - tau_f) / L_chip
  3. ng(lam) linear fitten, Gl. (4):  ng(lam) ~ a - b*lam0 - c*lam0*(lam-lam0)
     mit a = neff(lam0), b = neff'(lam0), c = neff''(lam0); daraus Gl. (5).
  4. Das Spektrum auf die Variable neff(lam)/lam umtasten und erst dann
     FFT -- damit verschwindet die dispersionsbedingte Verbreiterung.

Was hier ANDERS ist als im Paper, und warum
-------------------------------------------
Bei Zhao bilden die beiden Chipfacetten selbst das Interferometer: der ganze
Pfad liegt im dispersiven Wellenleiter. Dann linearisiert eine einzige
Umtastung auf u = neff(lam)/lam die Phase 4*pi*d*u fuer JEDE Tiefe d
gleichzeitig -- der ganze Chip wird auf einmal scharf.

Hier liegt zwischen LO und Chip rund 1.6 m Faser. Die Phase hat zwei Terme
mit verschiedener lam-Abhaengigkeit:

    phi(lam) = (2*pi/lam) * [ n_faser(lam)*Lambda_faser + 2*n_chip(lam)*d ]

Der Faseranteil ist bereits exakt kompensiert -- nicht durch dieses Skript,
sondern durch das Aux-Resampling: die Aux ist dieselbe Fasersorte, ihre Phase
ist (2*pi/lam)*n_faser(lam)*Lambda_aux, also STRENG PROPORTIONAL zum
Faseranteil der Messstrecke. Auf der Aux-Phasenachse ist jeder reine
Faserreflektor dispersionsfrei, unabhaengig davon, wie verschieden die beiden
Armlaengen sind. Genau das zeigen die Daten: Facette A steht ueber den ganzen
Sweep still, Facette B wandert.

Wuerde man nun wie im Paper global auf u = n_chip(lam)/lam umtasten, wuerde
der Faserterm nichtlinear -- und zwar auf einer 1.6-m-Strecke, also um
Zentimeter. Deshalb wird hier stattdessen der DIFFERENZANTEIL als Phase
abgezogen:

    tau_chip(lam) = 2*d*ng(lam)/c
    dphi(nu)      = 2*pi * Integral [ tau_chip(nu') - tau_chip(nu0) ] dnu'
    y_korrigiert  = analytisch(y) * exp(+i*dphi)

Das ist mathematisch dieselbe Operation wie die Umtastung des Papers, nur auf
den Chipanteil beschraenkt. Der Preis: sie ist auf EINE Tiefe d scharf
gestellt (HANDOVER.md Abschnitt 6, Punkt 5 sagt genau das voraus). Deshalb
faehrt `--focus-scan` d durch und sucht das Optimum -- ein Autofokus, der
zugleich das Produkt d*dng/dlam misst, unabhaengig vom linearen Fit.

Kontrollen, ohne die kein Befund gilt
-------------------------------------
* Facette A ist die Negativkontrolle der Schritte 1-3: sie liegt bei d = 0
  und darf NICHT dispersiv wandern. Tut sie es doch, misst man einen Artefakt.
* Das Spiegelfenster (gleicher Abstand auf der Faserseite von Facette A) und
  ein zweiter Scan auf einem leeren Kanal (`--control`) sagen, ob ein Peak
  Kanalinhalt ist oder Instrument.

Aufruf
    python tools/ofdr_dispersion_deembed.py \
        --scan raw_data/2026-09-17-10-33_2680_fiber48_20dbsplitteramanfang.npz \
        --control raw_data/2026-09-17-10-33_2680_fiber1_20dbsplitteramanfang.npz \
        --l-chip-um 7144.09 \
        --marks-csv results/2026-09-16/ligentec_kanal48_struktur_bridged.csv \
        --tag 2680_fiber48_zhao --out-dir results/2026-09-17
"""

import argparse
import csv as _csv
import datetime
import json
import os
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import find_peaks, zoom_fft
from scipy.signal.windows import blackmanharris, hann, kaiser

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import process_reflectogram_aux as pra  # noqa: E402

ROOT = os.path.dirname(_HERE)
C, NG = pra.C, pra.NG
NG_MARKS_CSV = 1.90      # der Gruppenindex, mit dem die Marken-CSV gerechnet ist


# ------------------------------------------------------------------ laden
def prep(path, trim=0.01):
    """Aux-resampelter Record MIT mitgefuehrter Wellenlaengenachse.

    `process_reflectogram_aux.resample_on_aux` gibt lam nicht zurueck, das
    gleitende Fenster von Schritt 1 braucht sie aber fuer jeden Abtastpunkt.
    Deshalb hier dieselbe Rechnung, nur mit lam auf demselben Phasengitter
    interpoliert. Alles andere (balanced, analytic, Trim, Monotonie) kommt
    unveraendert aus dem Hauptmodul.
    """
    d = np.load(path, allow_pickle=True)
    meta = json.loads(str(d["meta"])) if "meta" in d.files else {}
    crop = meta.get("crop")
    sl0 = slice(int(crop[0]), int(crop[1])) if crop else slice(None)
    lam_raw = d["wavelength_nm"][sl0]
    n = len(lam_raw)
    tau_aux = meta["wavelength_axis_aux"]["tau_aux_implied_s"]

    aux, _ = pra.balanced(d["ch2"][sl0], d["ch4"][sl0])
    meas, _ = pra.balanced(d["ch1"][sl0], d["ch3"][sl0])

    phi = np.unwrap(np.angle(pra.analytic(aux)))
    if phi[-1] < phi[0]:
        phi = -phi
    k = max(1, int(trim * n))
    sl = slice(k, n - k)
    phi, meas, lam_raw = phi[sl], meas[sl], lam_raw[sl]
    eps = 1e-6 * (phi[-1] - phi[0]) / len(phi)
    phi = np.maximum.accumulate(phi) + np.arange(len(phi)) * eps

    m = len(phi)
    phi_u = np.linspace(phi[0], phi[-1], m)
    y = PchipInterpolator(phi, meas)(phi_u)
    lam = PchipInterpolator(phi, lam_raw)(phi_u)
    t = np.linspace(-1, 1, m)
    y = y - np.polyval(np.polyfit(t, y, 5), t)

    span_nu = (phi[-1] - phi[0]) / (2 * np.pi * tau_aux)
    return dict(path=path, y=y, lam=lam, m=m, nu=C / (lam * 1e-9),
                dnu=(phi_u[1] - phi_u[0]) / (2 * np.pi * tau_aux),
                tau_aux=tau_aux, span_nu=span_nu,
                dz_cell=C / (2 * NG * span_nu))


def make_window(name, m, beta):
    return {"hann": hann, "blackmanharris": blackmanharris,
            "kaiser": lambda k: kaiser(k, beta)}[name](m)


def zspec(sig, dnu, zlo, zhi, nz, win=None):
    """|FFT| auf einem feinen z-Gitter [zlo, zhi] per chirp-z.

    Ein gezoomtes Spektrum statt einer genullten FFT: bei 2.4 Mio. Punkten
    kostet 0.5-um-Aufloesung sonst eine 200-Mio.-Punkt-FFT.
    """
    w = 1.0 if win is None else win
    s = len(sig) if win is None else win.sum()
    f_lo = 2 * NG * dnu * zlo / C
    f_hi = 2 * NG * dnu * zhi / C
    X = np.abs(zoom_fft(sig * w, [f_lo, f_hi], m=nz, fs=1.0)) / s
    return np.linspace(zlo, zhi, nz), X


def full_spec(rec, win, pad=2):
    """Gewoehnliche rFFT ueber die ganze Achse -- fuer die Uebersichtspanels.

    `pad` interpoliert nur die Peakform, es verbessert die Aufloesung nicht.
    Ohne Polsterung liegt das Bingitter genau auf dem Zellmass und der
    Facettenpeak wird um bis zu 0.4 dB zu niedrig abgelesen -- das reicht,
    um die Uebersichtspanels gegen die Chippanels zu verschieben.
    """
    n = pad * rec["m"]
    R = np.abs(np.fft.rfft(rec["y"] * win, n=n)) / win.sum()
    z = np.arange(len(R)) * C / (2 * NG * rec["dnu"] * n)
    return z, R


# ------------------------------------------- Schritt 1: gleitendes Fenster
def track_facets(rec, za_win, zb_win, w_nm, step_nm, snr_db):
    """Gruppenlaufzeit beider Facetten je Wellenlaengenfenster.

    Rechteckfenster wie im Paper. Zurueck kommt eine Tabelle
    (lam_c, z_A, z_B, snr_A, snr_B) plus die Maske der brauchbaren Fenster.
    `snr` ist der Peak ueber dem Median desselben Suchfensters -- ein
    lokaler Kontrast, kein globaler Rauschboden, denn genau dort steht der
    Multipath-Wald.
    """
    lam, y, dnu = rec["lam"], rec["y"], rec["dnu"]
    centers = np.arange(lam[0] + w_nm / 2, lam[-1] - w_nm / 2 + 1e-9, step_nm)
    rows = []
    for lc in centers:
        s = (lam >= lc - w_nm / 2) & (lam <= lc + w_nm / 2)
        i0 = int(np.argmax(s))
        i1 = int(len(s) - np.argmax(s[::-1]))
        seg = y[i0:i1]
        zs, snrs = [], []
        for lo, hi in (za_win, zb_win):
            nz = max(400, int((hi - lo) / 0.25e-6))
            zz, X = zspec(seg, dnu, lo, hi, nz)
            j = int(np.argmax(X))
            zs.append(zz[j])
            snrs.append(20 * np.log10(X[j] / np.median(X)))
        # Spaltenordnung: lam, z_A, z_B, snr_A, snr_B -- ng_fit und die
        # Kontrastschwelle lesen genau diese Reihenfolge.
        rows.append([lc] + zs + snrs)
    r = np.array(rows)
    return r, r[:, 4] >= snr_db


def ng_fit(r, good, l_chip, lam0):
    """Gl. (3) bis (5): ng(lam), linearer Fit, neff''(lam0), D."""
    ng = (r[:, 2] - r[:, 1]) * NG / l_chip               # Gl. (3)
    S, A = np.polyfit(r[good, 0] - lam0, ng[good], 1)    # Gl. (4)
    resid = ng[good] - (A + S * (r[good, 0] - lam0))
    return dict(ng=ng, A=A, S=S, resid_rms=float(resid.std()),
                neff2=-S / lam0,                         # Gl. (5): c = neff''
                D_ps_nm_km=S / C * 1e15,
                lam_lo=float(r[good, 0].min()),
                lam_hi=float(r[good, 0].max()), n_used=int(good.sum()))


# ----------------------------------------------- Schritt 4: De-Embedding
def phase_ramp(rec, A, S, lam0):
    """Das Phasenintegral aus der Kopfzeile, normiert auf 1 m Fokustiefe.

    Bezug ist tau(lam0), nicht der Mittelwert -- dann bleibt die Position
    eines Peaks die, die er bei lam0 hat, und die Skala heisst weiterhin
    "Gruppenindex bei lam0".
    """
    tau1 = (2.0 / C) * (A + S * (rec["lam"] - lam0))     # pro Meter Tiefe
    tau0 = (2.0 / C) * A
    f = 0.5 * (tau1[1:] + tau1[:-1]) - tau0
    return 2 * np.pi * np.concatenate(([0.0],
                                       np.cumsum(f * np.diff(rec["nu"]))))


def deembed(ya, ramp, d_focus):
    return ya if d_focus == 0.0 else ya * np.exp(1j * d_focus * ramp)


def width_db(zz, X, j, lev=3.0):
    """-lev-dB-Breite um den Index j, in um. NaN, wenn der Peak nicht auf
    beiden Seiten unter die Schwelle faellt -- dann ist es kein Peak."""
    th = X[j] * 10 ** (-lev / 20.0)
    L = np.where(X[:j] < th)[0]
    R = np.where(X[j:] < th)[0]
    if not len(L) or not len(R):
        return np.nan
    return (zz[j + R[0]] - zz[L[-1]]) * 1e6


# ------------------------------------------------------------------ plots
def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_ng(r, good, fit, lam0, path, title):
    plt = _plt()
    fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True,
                           gridspec_kw=dict(height_ratios=[2, 1]))
    ng = fit["ng"]
    ax[0].plot(r[~good, 0], ng[~good], ".", ms=3.5, color="0.72",
               label="verworfen (Facette B unter der Kontrastschwelle)")
    ax[0].plot(r[good, 0], ng[good], ".", ms=5, color="tab:blue",
               label="benutzt")
    xs = np.linspace(fit["lam_lo"], fit["lam_hi"], 50)
    ax[0].plot(xs, fit["A"] + fit["S"] * (xs - lam0), "-", lw=2,
               color="magenta",
               label=("linearer Fit, Gl. (4):  ng = %.5f %+.3e (lam-%.0f)"
                      % (fit["A"], fit["S"], lam0)))
    ax[0].set_ylabel("Gruppenindex  ng(lam)")
    ax[0].set_ylim(np.percentile(ng[good], 1) - 0.02,
                   np.percentile(ng[good], 99) + 0.02)
    ax[0].legend(fontsize=8, loc="upper right")
    ax[0].grid(alpha=0.3)
    ax[0].set_title(title, fontsize=10.5)
    ax[1].plot(r[:, 0], r[:, 4], ".-", ms=3, lw=0.7, color="tab:red",
               label="Facette B")
    ax[1].plot(r[:, 0], r[:, 3], ".-", ms=3, lw=0.7, color="tab:green",
               label="Facette A (Kontrolle, d = 0)")
    ax[1].set_xlabel("Fenstermitte lambda [nm]")
    ax[1].set_ylabel("Peak ueber\nlokalem Median [dB]")
    ax[1].grid(alpha=0.3)
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote: %s" % path)


def zoom_ladder(z_end, inner_lo, inner_hi, k=4):
    """k Fenster, geometrisch von der ganzen Achse bis inner.

    Geometrisch und nicht linear: von 3 m auf 15 mm ist der lineare
    Mittelwert 1.5 m und zeigte zweimal fast dasselbe Bild.
    """
    spans = np.geomspace(z_end, inner_hi - inner_lo, k)
    mid = 0.5 * (inner_lo + inner_hi)
    out = [(0.0, z_end)]
    for sp in spans[1:]:
        lo = max(mid - sp / 2, 0.0)
        out.append((lo, min(lo + sp, z_end)))
    return out


def plot_ladder(panels, za, marks, path, title, ymin):
    """Vier Zoomstufen, x = absolute Distanz vom Instrument."""
    plt = _plt()
    n = len(panels)
    fig, axes = plt.subplots(n, 1, figsize=(13, 2.5 * n + 1.0))
    axes = np.atleast_1d(axes)
    for ax, (win, traces) in zip(axes, panels):
        span = win[1] - win[0]
        for lbl, z, db, col, lw in traces:
            s = (z >= win[0]) & (z <= win[1])
            ax.plot(z[s] * 1e3, db[s], lw=lw, color=col, label=lbl)
        ax.axvline(za * 1e3, color="k", lw=0.9, alpha=0.75)
        if span < 0.060:
            # Einzelmarken erst, wenn sie nicht zu einem Strich verschmelzen
            for dA, lbl in marks:
                ax.axvline((za + dA * 1e-3) * 1e3, color="0.55", lw=0.7,
                           ls=":", zorder=0)
        else:
            ax.axvspan(za * 1e3, (za + 0.0148) * 1e3, color="0.85",
                       zorder=0, lw=0)
        ax.set_xlim(win[0] * 1e3, win[1] * 1e3)
        ax.set_ylim(ymin, 3)
        ax.grid(alpha=0.28)
        ax.set_ylabel("dB")
        ax.text(0.004, 0.95, "Fenster %.1f mm" % (span * 1e3),
                transform=ax.transAxes, va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="gray"))
        ax.legend(fontsize=7.5, loc="upper right", ncol=2)
    axes[-1].set_xlabel("absolute Distanz vom Instrument [mm] "
                        "(einfacher Weg / Reflexionskonvention)")
    fig.suptitle(title, fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote: %s" % path)


def plot_map(centers, M, zs, za, path, title):
    plt = _plt()
    db = 20 * np.log10(M / M.max() + 1e-18)
    fig, ax = plt.subplots(figsize=(12, 6))
    im = ax.pcolormesh((zs - za) * 1e3, centers, db, cmap="magma",
                       vmin=np.percentile(db, 55), vmax=db.max(),
                       shading="auto")
    ax.set_xlabel("Abstand hinter Facette A [mm]")
    ax.set_ylabel("Fenstermitte lambda [nm]")
    ax.set_title(title, fontsize=10.5)
    fig.colorbar(im, ax=ax, label="dB")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote: %s" % path)


# ------------------------------------------------------------------- main
def build_argparser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scan", required=True)
    p.add_argument("--control", default=None,
                   help="zweiter Scan, leerer Kanal -- die Nullkontrolle")
    p.add_argument("--l-chip-um", type=float, required=True,
                   help="Wellenleiterlaenge Facette A -> Facette B [um]; "
                        "geht nur in die Skala von ng ein (Gl. 3), nicht in "
                        "die Korrektur selbst")
    p.add_argument("--lam0", type=float, default=1580.0)
    p.add_argument("--window-nm", type=float, default=20.0,
                   help="Breite des gleitenden Rechteckfensters (Paper: 20)")
    p.add_argument("--step-nm", type=float, default=1.0,
                   help="Schrittweite des Fensters (Paper: 0.2)")
    p.add_argument("--snr-db", type=float, default=14.0,
                   help="Mindestkontrast von Facette B, damit ein Fenster in "
                        "den Fit eingeht")
    p.add_argument("--back-search", default="5:14",
                   help="Suchbereich fuer Facette B hinter Facette A [mm]")
    p.add_argument("--focus-scan", default="0:10:0.25",
                   help="Autofokus LO:HI:STEP in mm Chiptiefe")
    p.add_argument("--chip-mm", type=float, default=14.0,
                   help="Tiefe des Chipfensters hinter Facette A [mm]")
    p.add_argument("--guard-mm", type=float, default=2.2,
                   help="Nahzone um Facette A, die ausgespart bleibt")
    p.add_argument("--marks-csv", default=None)
    p.add_argument("--window", default="kaiser",
                   choices=["hann", "blackmanharris", "kaiser"])
    p.add_argument("--kaiser-beta", type=float, default=12.0)
    p.add_argument("--trim", type=float, default=0.01)
    p.add_argument("--tag", default="zhao")
    p.add_argument("--out-dir", default=None)
    p.add_argument("--ymin", type=float, default=-85.0)
    return p


def chip_view(rec, ramp, d_focus, za, chip_mm, win):
    """Reflektogramm vor und nach der Korrektur, feines Gitter um den Chip."""
    ya = pra.analytic(rec["y"])
    lo, hi = za - chip_mm * 1e-3, za + chip_mm * 1e-3
    nz = int((hi - lo) / 0.5e-6) + 1
    zz, X0 = zspec(ya, rec["dnu"], lo, hi, nz, win)
    _, X1 = zspec(deembed(ya, ramp, d_focus), rec["dnu"], lo, hi, nz, win)
    # Das analytische Signal traegt die doppelte Amplitude der rFFT des
    # reellen Interferogramms (die negativen Frequenzen sind hineingefaltet).
    # Ohne dieses Halbieren staenden die Chippanels 6.02 dB ueber den
    # Uebersichtspanels aus `full_spec` -- auf EINER Achse ein echter Fehler.
    X0, X1 = X0 / 2.0, X1 / 2.0
    return dict(za=za, zz=zz, X0=X0, X1=X1, ref=X0.max())


def main():
    a = build_argparser().parse_args()
    out_dir = a.out_dir or os.path.join(ROOT, "results",
                                        datetime.date.today().isoformat())
    os.makedirs(out_dir, exist_ok=True)
    lam0, l_chip = a.lam0, a.l_chip_um * 1e-6
    bs_lo, bs_hi = (float(v) * 1e-3 for v in a.back_search.split(":"))
    fs_lo, fs_hi, fs_st = (float(v) for v in a.focus_scan.split(":"))

    print("Zhao et al. 2017 -- Dispersions-De-Embedding")
    print("scan    : %s" % a.scan)
    print("control : %s" % (a.control or "-- keine --"))

    rec = prep(a.scan, a.trim)
    win = make_window(a.window, rec["m"], a.kaiser_beta)
    print("\n%d Punkte, %.3f .. %.3f nm, dnu %.4f MHz, Zelle %.2f um, "
          "Nyquist %.3f m"
          % (rec["m"], rec["lam"][0], rec["lam"][-1], rec["dnu"] / 1e6,
             rec["dz_cell"] * 1e6, C / (4 * NG * rec["dnu"])))

    zfull, Rfull = full_spec(rec, win)
    za = zfull[int(np.argmax(Rfull))]
    print("Facette A bei %.4f mm" % (za * 1e3))

    zz, X = zspec(rec["y"], rec["dnu"], za + bs_lo, za + bs_hi,
                  int((bs_hi - bs_lo) / 0.5e-6), win)
    zb = zz[int(np.argmax(X))]
    print("Facette B bei %.4f mm  (%.4f mm hinter A)  ->  ng(grob) = %.4f"
          % (zb * 1e3, (zb - za) * 1e3, (zb - za) * NG / l_chip))

    # ---------------------------------------------------- Schritt 1 bis 3
    r, good = track_facets(rec, (za - 200e-6, za + 200e-6),
                           (zb - 400e-6, zb + 400e-6),
                           a.window_nm, a.step_nm, a.snr_db)
    fit = ng_fit(r, good, l_chip, lam0)
    za_drift = float(np.ptp(r[good, 1]) * 1e6)
    zb_drift = float(np.ptp(r[good, 2]) * 1e6)
    print("\nSchritt 1-3 (Gl. 3-5): %d von %d Fenstern brauchbar "
          "(%.1f .. %.1f nm)"
          % (fit["n_used"], len(r), fit["lam_lo"], fit["lam_hi"]))
    print("  ng(%.0f nm)   = %.5f" % (lam0, fit["A"]))
    print("  dng/dlam     = %.4e /nm   (Residuum rms %.2e)"
          % (fit["S"], fit["resid_rms"]))
    print("  neff''(lam0) = %.4e /nm^2" % fit["neff2"])
    print("  D            = %.0f ps/(nm*km)" % fit["D_ps_nm_km"])
    print("  Wanderung ueber das benutzte Band:  Facette B %.1f um,  "
          "Facette A %.1f um (Kontrolle, sollte klein sein)"
          % (zb_drift, za_drift))

    with open(os.path.join(out_dir, "%s_ng_vs_lambda.csv" % a.tag), "w",
              newline="") as f:
        f.write("lambda_nm,z_facetA_mm,z_facetB_mm,dz_mm,ng,"
                "snr_A_dB,snr_B_dB,used\n")
        for k in range(len(r)):
            f.write("%.3f,%.6f,%.6f,%.6f,%.6f,%.2f,%.2f,%d\n"
                    % (r[k, 0], r[k, 1] * 1e3, r[k, 2] * 1e3,
                       (r[k, 2] - r[k, 1]) * 1e3, fit["ng"][k],
                       r[k, 3], r[k, 4], int(good[k])))

    # -------------------------------------------- Schritt 4 und Autofokus
    ramp = phase_ramp(rec, fit["A"], fit["S"], lam0)
    ya = pra.analytic(rec["y"])
    focs = np.arange(fs_lo, fs_hi + 1e-9, fs_st) * 1e-3
    print("\nAutofokus ueber %d Tiefen (%.2f .. %.2f mm)"
          % (len(focs), fs_lo, fs_hi))
    rows, best = [], None
    for dfoc in focs:
        z2, X2 = zspec(deembed(ya, ramp, dfoc), rec["dnu"],
                       zb - 300e-6, zb + 300e-6, 1201, win)
        j = int(np.argmax(X2))
        rows.append((dfoc, X2[j], z2[j], width_db(z2, X2, j)))
        if best is None or X2[j] > best[1]:
            best = rows[-1]
    d_opt = best[0]
    S_af = d_opt * fit["S"] / l_chip if l_chip else np.nan
    print("  Optimum bei d = %.2f mm:  Facette B %.4f mm hinter A,  "
          "-3 dB = %.1f um" % (d_opt * 1e3, (best[2] - za) * 1e3, best[3]))
    print("  Der Autofokus misst das Produkt d*dng/dlam. Mit L_chip = "
          "%.1f um folgt dng/dlam = %.4e /nm, D = %.0f ps/(nm*km)"
          % (a.l_chip_um, S_af, S_af / C * 1e15))

    with open(os.path.join(out_dir, "%s_autofocus.csv" % a.tag), "w",
              newline="") as f:
        f.write("d_focus_mm,peak_rel_dB,z_mm,w3_um\n")
        ref0 = rows[0][1]
        for dfoc, pk, z, w3 in rows:
            f.write("%.4f,%.4f,%.6f,%.2f\n"
                    % (dfoc * 1e3, 20 * np.log10(pk / ref0), z * 1e3, w3))

    # -------------------------------------------------- Chip, vor / nach
    res = chip_view(rec, ramp, d_opt, za, a.chip_mm, win)
    ctl = None
    if a.control:
        rc = prep(a.control, a.trim)
        wc = make_window(a.window, rc["m"], a.kaiser_beta)
        zc, Rc = full_spec(rc, wc)
        zac = zc[int(np.argmax(Rc))]
        ctl = chip_view(rc, phase_ramp(rc, fit["A"], fit["S"], lam0),
                        d_opt, zac, a.chip_mm, wc)
        print("\nNullkontrolle: Facette A bei %.4f mm" % (zac * 1e3))

    zz, ref = res["zz"], res["ref"]
    d0 = 20 * np.log10(res["X0"] / ref + 1e-18)
    d1 = 20 * np.log10(res["X1"] / ref + 1e-18)
    dA = (zz - za) * 1e3
    dc, dAc = None, None
    if ctl:
        dc = 20 * np.log10(ctl["X1"] / ctl["ref"] + 1e-18)
        dAc = (ctl["zz"] - ctl["za"]) * 1e3

    def at(arr, ax, x, halfw=0.030):
        s = (ax > x - halfw) & (ax < x + halfw)
        return float(arr[s].max()) if s.sum() else np.nan

    pk, _ = find_peaks(d1, distance=10)
    pk = pk[(dA[pk] >= a.guard_mm) & (dA[pk] <= a.chip_mm)]
    pk = pk[np.argsort(d1[pk])[::-1][:18]]
    csv_path = os.path.join(out_dir, "%s_peaks.csv" % a.tag)
    print("\nPeaks im Chipfenster nach De-Embedding "
          "(alle dB relativ zur eigenen Facette A)")
    print("   dA [mm]   unkorr    fokus  Gewinn  -3dB um |  Spiegel  "
          "Kontrolle |  Ueberschuss Spiegel / Kontrolle")
    with open(csv_path, "w", newline="") as f:
        f.write("dA_mm,z_abs_mm,uncorrected_dB,deembedded_dB,gain_dB,w3_um,"
                "mirror_dB,control_dB,excess_mirror_dB,excess_control_dB\n")
        for j in sorted(pk, key=lambda k: dA[k]):
            mir = at(d1, dA, -dA[j])
            cv = at(dc, dAc, dA[j]) if ctl else np.nan
            w3 = width_db(zz, res["X1"], j)
            f.write("%.4f,%.4f,%.3f,%.3f,%.3f,%.2f,%.3f,%.3f,%.3f,%.3f\n"
                    % (dA[j], zz[j] * 1e3, d0[j], d1[j], d1[j] - d0[j], w3,
                       mir, cv, d1[j] - mir, d1[j] - cv))
            print("  %8.3f  %7.2f  %7.2f  %+6.2f  %6.1f  |  %7.2f  %8.2f  |"
                  "   %+6.2f / %+6.2f"
                  % (dA[j], d0[j], d1[j], d1[j] - d0[j], w3, mir, cv,
                     d1[j] - mir, d1[j] - cv))
    print("wrote: %s" % csv_path)

    # ------------- greift die Korrektur den Peak an oder auch den Sockel?
    chip = (dA >= a.guard_mm) & (dA <= a.chip_mm)
    mirr = (dA <= -a.guard_mm) & (dA >= -a.chip_mm)
    print("\nChipfenster gegen Spiegelfenster (%.1f .. %.1f mm, beide gleich "
          "breit)" % (a.guard_mm, a.chip_mm))
    print("                       max     Median (Sockel)")
    for lbl, dd in (("unkorrigiert", d0), ("de-embedded ", d1)):
        print("  %s  Chip    %7.2f  %7.2f"
              % (lbl, dd[chip].max(), np.median(dd[chip])))
        print("  %s  Spiegel %7.2f  %7.2f"
              % (" " * len(lbl), dd[mirr].max(), np.median(dd[mirr])))
    print("  Ueberschuss unkorrigiert: max %+6.2f   Sockel %+6.2f"
          % (d0[chip].max() - d0[mirr].max(),
             np.median(d0[chip]) - np.median(d0[mirr])))
    print("  Ueberschuss de-embedded : max %+6.2f   Sockel %+6.2f"
          % (d1[chip].max() - d1[mirr].max(),
             np.median(d1[chip]) - np.median(d1[mirr])))

    with open(os.path.join(out_dir, "%s_kennzahlen.csv" % a.tag), "w",
              newline="") as f:
        f.write("groesse,wert\n")
        for k, v in (("facetA_mm", za * 1e3), ("facetB_mm", zb * 1e3),
                     ("dz_AB_mm", (zb - za) * 1e3),
                     ("L_chip_um", a.l_chip_um),
                     ("lam0_nm", lam0), ("ng_lam0", fit["A"]),
                     ("dng_dlam_fit_per_nm", fit["S"]),
                     ("neff2_per_nm2", fit["neff2"]),
                     ("D_ps_nm_km_fit", fit["D_ps_nm_km"]),
                     ("d_focus_opt_mm", d_opt * 1e3),
                     ("dng_dlam_autofocus_per_nm", S_af),
                     ("D_ps_nm_km_autofocus", S_af / C * 1e15),
                     ("facetA_drift_um", za_drift),
                     ("facetB_drift_um", zb_drift),
                     ("windows_used", fit["n_used"]),
                     ("chip_max_uncorr_dB", d0[chip].max()),
                     ("chip_max_deemb_dB", d1[chip].max()),
                     ("chip_median_uncorr_dB", np.median(d0[chip])),
                     ("chip_median_deemb_dB", np.median(d1[chip])),
                     ("mirror_median_uncorr_dB", np.median(d0[mirr])),
                     ("mirror_median_deemb_dB", np.median(d1[mirr]))):
            f.write("%s,%.6g\n" % (k, v))

    # ------------------- Bauteilmarken, mit dem GEMESSENEN ng umgerechnet
    marks = []
    if a.marks_csv and os.path.exists(a.marks_csv):
        with open(a.marks_csv) as f:
            for row in _csv.DictReader(f):
                s_um = float(row["s_from_facet_um"])
                marks.append((s_um * fit["A"] / NG * 1e-3, row["cell"]))
        print("\n%d Bauteilmarken, mit dem gemessenen ng = %.4f umgerechnet "
              "statt der im CSV unterstellten %.2f (Faktor %.4f)"
              % (len(marks), fit["A"], NG_MARKS_CSV, fit["A"] / NG_MARKS_CSV))
        print("  Bauteilbereich wandert von %.3f..%.3f mm nach %.3f..%.3f mm"
              % (min(m for m, _ in marks) * NG_MARKS_CSV / fit["A"],
                 max(m for m, _ in marks) * NG_MARKS_CSV / fit["A"],
                 min(m for m, _ in marks), max(m for m, _ in marks)))

    # ----------------------------------------------------------- Plots
    plot_ng(r, good, fit, lam0,
            os.path.join(out_dir, "%s_ng_vs_lambda.png" % a.tag),
            "Gruppenindex aus Gl. (3), Fenster %.0f nm / Schritt %.1f nm  --  "
            "%s" % (a.window_nm, a.step_nm, os.path.basename(a.scan)))

    dbf = 20 * np.log10(Rfull / ref + 1e-18)
    panels = []
    for wlo, whi in zoom_ladder(zfull[-1], za - 0.001, za + 0.0148):
        if (whi - wlo) > 0.060:
            tr = [("unkorrigiert", zfull, dbf, "0.45", 0.55)]
        else:
            tr = [("unkorrigiert", zz, d0, "0.55", 0.7),
                  ("de-embedded (Fokus %.2f mm)" % (d_opt * 1e3), zz, d1,
                   "tab:red", 0.85)]
            if ctl:
                tr.append(("Nullkontrolle, de-embedded",
                           ctl["zz"] - ctl["za"] + za, dc, "tab:blue", 0.55))
        panels.append(((wlo, whi), tr))
    plot_ladder(panels, za, marks,
                os.path.join(out_dir, "%s_zoom4.png" % a.tag),
                "De-Embedding nach Zhao et al. 2017  --  %s\n"
                "ng(%.0f nm) = %.5f,  dng/dlam = %.3e /nm,  D = %.0f "
                "ps/(nm km),  Fokus %.2f mm"
                % (os.path.basename(a.scan), lam0, fit["A"], fit["S"],
                   fit["D_ps_nm_km"], d_opt * 1e3), a.ymin)

    centers = np.arange(rec["lam"][0] + a.window_nm / 2,
                        rec["lam"][-1] - a.window_nm / 2 + 1e-9, 2.0)
    zs = np.linspace(za - 0.0015, za + a.chip_mm * 1e-3, 3000)
    M = np.empty((len(centers), len(zs)))
    for i, lc in enumerate(centers):
        s = ((rec["lam"] >= lc - a.window_nm / 2)
             & (rec["lam"] <= lc + a.window_nm / 2))
        i0 = int(np.argmax(s))
        i1 = int(len(s) - np.argmax(s[::-1]))
        _, M[i] = zspec(rec["y"][i0:i1], rec["dnu"], zs[0], zs[-1], len(zs))
    plot_map(centers, M, zs, za,
             os.path.join(out_dir, "%s_lambda_z_map.png" % a.tag),
             "Reflektogramm je Wellenlaengenfenster (Schritt 1 des Papers)"
             "  --  %s\nFacette B wandert mit lambda, Facette A nicht: "
             "das IST die Wellenleiterdispersion"
             % os.path.basename(a.scan))
    print("\nfertig.")


if __name__ == "__main__":
    main()
