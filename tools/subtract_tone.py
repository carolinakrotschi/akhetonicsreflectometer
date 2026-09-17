#!/usr/bin/env python3
"""Den staerksten Sinus aus dem Spektrum nehmen und neu transformieren.

Hintergrund: die Vermutung, dass der Untergrund im Reflektogramm gar kein
additives Rauschen ist, sondern die *Unvollkommenheit der einzelnen
Sinusse* -- ein Reflektor, dessen Schwebung nicht sauber monofrequent ist,
schmiert ueber die ganze Laengenachse. Dann muesste der Boden fallen,
sobald man den dominanten Ton aus dem Frequenz-Domain-Signal y(nu)
herausnimmt und erst DANACH fenstert und FFTt.

Das Skript nimmt genau diesen Ton weg, auf zwei Arten:

  ls    reiner Ton: gewichtete Kleinste-Quadrate-Anpassung von
        a*cos(2*pi*k*n/m) + b*sin(...) bei einer auf ~1e-4 Bins genau
        gesuchten Frequenz. Entfernt den IDEALEN Sinus. Was uebrig
        bleibt, ist alles, was dieser Reflektor NICHT als sauberen Ton
        beitraegt -- also genau die vermutete Verschmierung.

  band  Ton samt Hof: das Signal wird auf den Ton heruntergemischt, auf
        +-`--bw-mm` tiefpassgefiltert und wieder hochgemischt. Entfernt
        den Peak MIT seiner naeheren Umgebung.

`--n-tones N` wiederholt das (CLEAN-artig): staerksten Peak suchen,
wegnehmen, naechsten suchen.

WICHTIG -- was das Ergebnis beweisen kann und was nicht:

  Die FFT ist linear. Was man von y(nu) abzieht, aendert das
  Reflektogramm nur um die Transformierte des Abgezogenen. Ein
  `band`-Schnitt kann den Boden deshalb NUR innerhalb des Schnitts
  senken; faellt er ausserhalb trotzdem, ist etwas falsch gerechnet.
  Aussagekraeftig ist `ls`: ein idealer Ton hat ausserhalb der
  Fenster-Nebenzipfel keine Energie, ein UNVOLLKOMMENER dagegen schon.
  Bleibt der Boden nach `ls` unveraendert, traegt dieser Reflektor seine
  Energie nicht als breiten Teppich, und der Untergrund kommt woanders
  her.

Die Normierung bleibt in ALLEN Kurven die des Originals (staerkster Peak
vor dem Abzug = 0 dB). Sonst wuerde das Wegnehmen des groessten Peaks die
ganze Spur um dessen Hoehe anheben und der Vergleich waere hinfaellig
(dieselbe Falle wie bei der VOA-Serie, s. plot_voa_series.py).

Aufruf
    python tools/subtract_tone.py raw_data/hhi_1/2026-09-11-12-16hhi1fiber30.json \
        --tau-aux-ns 20.3809 --label "HHI-1 fiber 30 = b2 (11.09.)" \
        --out-dir results/2026-09-17 --cache <scratch>/f30_y.npz
    # zweiter Lauf, ohne die 157-MB-JSON noch einmal zu lesen:
    ... --cache <scratch>/f30_y.npz --method band --bw-mm 0.5
"""

import argparse
import os
import sys

import numpy as np
from scipy.signal.windows import blackmanharris, hann, kaiser

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import process_reflectogram_aux as pra  # noqa: E402

C = 299_792_458.0
NG = 1.468

# Fenster, in denen der Rauschboden getrennt berichtet wird. Getrennt
# deshalb, weil der Boden NICHT flach ist: 0-700 mm ist Faser, 2600-3000 mm
# steigt wieder an, dazwischen liegt das wirklich ruhige Stueck.
FLOOR_WINDOWS = [(0, 200), (200, 700), (700, 1200), (1200, 1600),
                 (1600, 1660), (1740, 2000), (2000, 2400), (2400, 2600),
                 (2600, 3000)]


# --------------------------------------------------------------- laden
def prepare(scan, tau_aux, trim, cache=None):
    """y(nu) auf gleichfoermiger Frequenzachse -- exakt die Kette aus
    process_reflectogram_aux.process() bis vor Fenster und FFT."""
    if cache and os.path.exists(cache):
        d = np.load(cache)
        print("cache: %s" % cache)
        return d["y"], float(d["dnu"]), float(d["span_nu"])

    ch1, ch2, ch3, ch4, _ = pra.load(scan)
    aux, _ = pra.balanced(ch2, ch4)
    meas, _ = pra.balanced(ch1, ch3)
    y, dnu, span_nu, diag = pra.resample_on_aux(meas, aux, tau_aux, trim)
    m = len(y)
    t = np.linspace(-1, 1, m)
    y = y - np.polyval(np.polyfit(t, y, 5), t)
    print("%s: %d Punkte, %.0f Aux-Fringes, %.1f Punkte/Fringe"
          % (os.path.basename(scan), m, diag["fringes"],
             diag["pts_per_fringe"]))
    if cache:
        np.savez_compressed(cache, y=y, dnu=dnu, span_nu=span_nu)
    return y, dnu, span_nu


# ------------------------------------------------- Ton finden und fitten
def refine_bin(y, w, k0, span=1.0, rounds=3, npts=64):
    """Frequenz des Tons in (gebrochenen) FFT-Bins.

    Zoom-FFT statt Parabel durch drei Punkte: bei -70 dB Boden entscheidet
    ein Hundertstel Bin darueber, ob der Abzug 30 oder 50 dB unterdrueckt.
    Drei Runden a 64 Punkte kommen auf ~1e-4 Bin.
    """
    from scipy.signal import ZoomFFT
    m = len(y)
    x = y * w
    lo, hi = k0 - span, k0 + span
    for _ in range(rounds):
        tr = ZoomFFT(m, [lo, hi], npts, fs=m, endpoint=True)
        mag = np.abs(tr(x))
        j = int(np.argmax(mag))
        step = (hi - lo) / (npts - 1)
        lo, hi = lo + max(j - 1, 0) * step, lo + min(j + 1, npts - 1) * step
    return 0.5 * (lo + hi)


def fit_tone(y, w, k):
    """Gewichtete KQ-Anpassung a*cos + b*sin bei Bin k. Gewicht = Fenster:
    gefittet wird der Ton so, wie ihn die spaetere FFT auch sieht."""
    m = len(y)
    ph = 2 * np.pi * k * np.arange(m) / m
    c, s = np.cos(ph), np.sin(ph)
    wc, ws = w * c, w * s
    A = np.array([[c @ wc, s @ wc], [c @ ws, s @ ws]])
    b = np.array([y @ wc, y @ ws])
    a, bb = np.linalg.solve(A, b)
    return a * c + bb * s


def fit_band(y, k, half_bins):
    """Ton MIT Hof: heruntermischen, tiefpass, hochmischen."""
    m = len(y)
    car = np.exp(-2j * np.pi * k * np.arange(m) / m)
    s = pra.analytic(y) * car
    S = np.fft.fft(s)
    B = max(1, int(round(half_bins)))
    S[B + 1:m - B] = 0.0
    return np.real(np.fft.ifft(S) * np.conj(car))


# ----------------------------------------------------------- Spektrum
def spectrum(y, win, dnu):
    m = len(y)
    R = np.abs(np.fft.rfft(y * win))
    z = np.arange(len(R)) * C / (2 * NG * dnu * m)
    return z, R


def to_db(R, ref):
    return 20 * np.log10(np.asarray(R) / ref + 1e-30)


def floor_table(z, db_list, names, zmax):
    """Markdown-Tabelle Rauschboden je Fenster, eine Spalte je Spur."""
    zm = z * 1e3
    rows = ["| Fenster [mm] | " + " | ".join(names) + " | Delta |",
            "|" + "---|" * (len(names) + 2)]
    for lo, hi in FLOOR_WINDOWS:
        if lo >= zmax * 1e3:
            continue
        sel = (zm >= lo) & (zm < min(hi, zmax * 1e3))
        if sel.sum() < 100:
            continue
        vals = [pra.noise_floor(d[sel]) for d in db_list]
        rows.append("| %d-%d | " % (lo, hi)
                    + " | ".join("%.2f" % v for v in vals)
                    + " | %+.2f |" % (vals[-1] - vals[0]))
    sel = zm < zmax * 1e3
    vals = [pra.noise_floor(d[sel]) for d in db_list]
    rows.append("| **gesamt 0-%.0f** | " % (zmax * 1e3)
                + " | ".join("**%.2f**" % v for v in vals)
                + " | **%+.2f** |" % (vals[-1] - vals[0]))
    return "\n".join(rows)


def width_table(z, db, zpk_mm, zmax):
    """Wie viel Energie des Peaks liegt innerhalb +-w um ihn herum?

    Die Frage "ist der Ton ein Ton?" laesst sich so ohne Modell
    beantworten: ein idealer Sinus steckt zu 99 % in einer Handvoll Bins,
    ein verschmierter braucht Millimeter.
    """
    zm = z * 1e3
    p = 10.0 ** (db / 10.0)
    tot = p[zm < zmax * 1e3].sum()
    inner = p[np.abs(zm - zpk_mm) < 50].sum()
    rows = ["| +-w um den Peak | Anteil an der Peak-Energie (+-50 mm) | "
            "Anteil an 0-%.0f mm |" % (zmax * 1e3), "|---|---|---|"]
    for w in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0):
        s = p[np.abs(zm - zpk_mm) < w].sum()
        rows.append("| %.2f mm | %6.2f %% | %6.2f %% |"
                    % (w, 100 * s / inner, 100 * s / tot))
    return "\n".join(rows)


def skirt_table(z, db, zpk_mm, zmax):
    """Rauschboden als Funktion des Abstands vom Peak -- das eigentliche
    Mass fuer 'wie weit schmiert dieser Reflektor?'."""
    zm = z * 1e3
    pk = db[np.abs(zm - zpk_mm) < 0.4].max()
    rows = ["| Abstand vom Peak [mm] | Boden [dB] | rel. zum Peak [dB] | Bins |",
            "|---|---|---|---|"]
    edges = [0.5, 1, 2, 5, 10, 20, 40, 80, 160, 320, 640, 1200]
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = ((np.abs(zm - zpk_mm) >= lo) & (np.abs(zm - zpk_mm) < hi)
               & (z <= zmax))
        if sel.sum() < 30:
            continue
        f = pra.noise_floor(db[sel])
        rows.append("| %g-%g | %.2f | %+.1f | %d |" % (lo, hi, f, f - pk,
                                                       sel.sum()))
    return "\n".join(rows)


def selftest(y, win, dnu, z, ref, bin_mm, z_probe_mm, amp):
    """Kontrolle: einen PERFEKTEN Ton bekannter Hoehe dazutun und wieder
    abziehen. Erst diese Zahl sagt, was eine schlechte Unterdrueckung des
    echten Tons bedeutet -- ohne sie koennte es auch am Verfahren liegen."""
    m = len(y)
    k_syn = z_probe_mm / bin_mm + 0.3137      # bewusst nicht auf dem Bingitter
    y2 = y + amp * np.cos(2 * np.pi * k_syn * np.arange(m) / m + 0.7)
    _, Rs = spectrum(y2, win, dnu)
    k0 = int(np.argmax(np.where(np.abs(z * 1e3 - z_probe_mm) < 5, Rs, 0.0)))
    kf = refine_bin(y2, win, float(k0))
    _, Rr = spectrum(y2 - fit_tone(y2, win, kf), win, dnu)
    before, after = to_db(Rs[k0], ref), to_db(Rr[k0], ref)
    print("Selbsttest: perfekter Ton bei %.4f mm (gefunden %.4f mm), "
          "%.2f dB -> %.2f dB  => %.1f dB unterdrueckt"
          % (k_syn * bin_mm, kf * bin_mm, before, after, before - after))
    return before - after


def ideal_peak(y, win, dnu, k_peak, amp):
    """Ein PERFEKTER Reflektor gleicher Hoehe am selben Ort, allein.

    Der fehlende Massstab. Erst im Vergleich mit dieser Kurve laesst sich
    sagen, was am Fuss eines Peaks Instrument ist (Fenster, Zahlenrauschen)
    und was echtes Licht. Gerechnet wird NUR der Ton, ohne Messdaten --
    was hier zu sehen ist, kann das Geraet gar nicht vermeiden.
    """
    m = len(y)
    syn = amp * np.cos(2 * np.pi * k_peak * np.arange(m) / m + 0.7)
    return spectrum(syn, win, dnu)[1]


def excess_table(z, db, db_ideal, zpk_mm, zmax):
    """Wie weit steht der Fuss des Peaks ueber dem eines idealen?"""
    zm = z * 1e3
    rows = ["| Abstand vom Peak [mm] | gemessen [dB] | idealer Reflektor [dB] "
            "| Ueberschuss |", "|---|---|---|---|"]
    edges = [0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40, 80]
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = ((np.abs(zm - zpk_mm) >= lo) & (np.abs(zm - zpk_mm) < hi)
               & (z <= zmax))
        if sel.sum() < 5:
            continue
        a = pra.noise_floor(db[sel])
        b = pra.noise_floor(db_ideal[sel])
        rows.append("| %g-%g | %.1f | %.1f | **%.1f dB** |"
                    % (lo, hi, a, b, a - b))
    return "\n".join(rows)


def peak_width(R, z, zc_mm, bin_mm):
    """-3-dB-Breite eines Peaks in um."""
    zm = z * 1e3
    idx = np.where(np.abs(zm - zc_mm) < 0.35)[0]
    i = idx[int(np.argmax(R[idx]))]
    half = R[i] / np.sqrt(2.0)
    l = r = i
    while l > 0 and R[l] > half:
        l -= 1
    while r < len(R) - 1 and R[r] > half:
        r += 1
    return (r - l) * bin_mm * 1000.0


def envelope(z, db, nout):
    """Blockweises Maximum auf ~nout Punkte, wie plot_chip_comparison.py."""
    k = max(len(z) // nout, 1)
    m = (len(z) // k) * k
    zz = z[:m].reshape(-1, k)
    dd = db[:m].reshape(-1, k)
    i = np.argmax(dd, axis=1)
    return zz[np.arange(len(i)), i], dd[np.arange(len(i)), i]


# --------------------------------------------------------------- Plot
def make_plot(z, db0, db1, zpk_mm, label, sub, out, spans, dbi=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    zm = z * 1e3
    fig, axs = plt.subplots(len(spans), 1, figsize=(13, 3.1 * len(spans)))
    axs = np.atleast_1d(axs)
    for ax, (lo, hi, title) in zip(axs, spans):
        sel = (zm >= lo) & (zm <= hi)
        a, b = zm[sel], db0[sel]
        c, d = zm[sel], db1[sel]
        if sel.sum() > 6000:
            a, b = envelope(a, b, 3000)
            c, d = envelope(c, d, 3000)
        ax.plot(a, b, lw=0.6, color="0.62", label="Original")
        ax.plot(c, d, lw=0.6, color="tab:blue", label="danach")
        if dbi is not None:
            e, f = zm[sel], dbi[sel]
            if sel.sum() > 6000:
                e, f = envelope(e, f, 3000)
            ax.plot(e, f, lw=0.9, color="crimson", ls="--",
                    label="idealer Reflektor gleicher Hoehe")
        f0 = pra.noise_floor(db0[sel])
        f1 = pra.noise_floor(db1[sel])
        ax.axhline(f0, color="0.55", ls=":", lw=1.0)
        ax.axhline(f1, color="tab:blue", ls=":", lw=1.0)
        ax.set_title("%s   (Boden %.1f -> %.1f dB)" % (title, f0, f1),
                     fontsize=9.5)
        ax.set_xlim(lo, hi)
        ax.set_ylabel("Amplitude [dB]")
        ax.grid(alpha=0.3)
        ax.axvline(zpk_mm, color="darkorange", lw=0.8, alpha=0.6)
    axs[0].legend(loc="upper right", fontsize=8)
    axs[-1].set_xlabel("absolute Distanz vom Instrument [mm]")
    fig.suptitle("%s -- %s" % (label, sub), fontsize=11)
    fig.text(0.01, 0.005,
             "Normierung in BEIDEN Kurven: staerkster Peak des Originals "
             "= 0 dB. Gepunktet: RMS-Rauschboden im jeweiligen Fenster. "
             "Orange: Lage des entfernten Tons.",
             fontsize=8, color="0.3")
    fig.tight_layout(rect=(0, 0.015, 1, 0.975))
    fig.savefig(out, dpi=150)
    print("geschrieben: %s" % out)


# --------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scan")
    ap.add_argument("--tau-aux-ns", type=float, default=20.3809)
    ap.add_argument("--trim", type=float, default=0.01)
    ap.add_argument("--cache", default=None,
                    help="npz mit y/dnu/span_nu; wird angelegt bzw. genutzt")
    ap.add_argument("--window", default="kaiser",
                    choices=["hann", "blackmanharris", "kaiser"])
    ap.add_argument("--kaiser-beta", type=float, default=12.0)
    ap.add_argument("--method", default="ls", choices=["ls", "band"])
    ap.add_argument("--bw-mm", type=float, default=0.5,
                    help="halbe Breite des Bandschnitts bei --method band")
    ap.add_argument("--n-tones", type=int, default=1)
    ap.add_argument("--ideal-overlay", action="store_true",
                    help="idealen Reflektor gleicher Hoehe am selben Ort "
                         "mitzeichnen -- der Massstab dafuer, was am Fuss "
                         "des Peaks Instrument ist und was echtes Licht")
    ap.add_argument("--selftest-mm", type=float, default=None,
                    help="Kontrolle: perfekten Ton dieser Hoehe an dieser "
                         "Stelle [mm] einspeisen und wieder abziehen")
    ap.add_argument("--zmax", type=float, default=3.0)
    ap.add_argument("--label", default=None)
    ap.add_argument("--out-dir", default="results")
    ap.add_argument("--name", default=None, help="Dateiname ohne Endung")
    a = ap.parse_args()

    label = a.label or os.path.basename(a.scan)
    y, dnu, span_nu = prepare(a.scan, a.tau_aux_ns * 1e-9, a.trim, a.cache)
    m = len(y)
    win = {"hann": hann(m), "blackmanharris": blackmanharris(m),
           "kaiser": kaiser(m, a.kaiser_beta)}[a.window]
    dz_bin = C / (2 * NG * span_nu)
    bin_mm = C / (2 * NG * dnu * m) * 1e3
    print("Aufloesung %.2f um, Plotgitter %.4f mm/Bin" % (dz_bin * 1e6, bin_mm))

    z, R0 = spectrum(y, win, dnu)
    ref = R0.max()
    db0 = to_db(R0, ref)
    lim = z <= a.zmax

    if a.selftest_mm is not None:
        selftest(y, win, dnu, z, ref, bin_mm, a.selftest_mm,
                 2 * R0.max() / win.sum())

    y_res = y.copy()
    removed = []
    for it in range(a.n_tones):
        _, Rc = spectrum(y_res, win, dnu)
        k0 = int(np.argmax(np.where(z <= a.zmax, Rc, 0.0)))
        kf = refine_bin(y_res, win, float(k0))
        z_pk = kf * C / (2 * NG * dnu * m)
        model = (fit_tone(y_res, win, kf) if a.method == "ls"
                 else fit_band(y_res, kf, a.bw_mm / bin_mm))
        y_res = y_res - model
        _, Rn = spectrum(y_res, win, dnu)
        before, after = to_db(Rc[k0], ref), to_db(Rn[k0], ref)
        removed.append((z_pk * 1e3, before, before - after))
        print("  Ton %2d: z = %10.4f mm (Bin %.4f), war %6.2f dB, "
              "Rest am Ort %6.2f dB  -> %.1f dB unterdrueckt"
              % (it + 1, z_pk * 1e3, kf, before, after, before - after))

    _, R1 = spectrum(y_res, win, dnu)
    db1 = to_db(R1, ref)

    sub = ("idealer Ton abgezogen (LS)" if a.method == "ls"
           else "Band +-%.2f mm herausgeschnitten" % a.bw_mm)
    if a.n_tones > 1:
        sub = "%d staerkste Toene, je %s" % (a.n_tones, sub)

    zpk = removed[0][0]
    dbi = None
    if a.ideal_overlay:
        k_pk = zpk * 1e-3 * (2 * NG * dnu * m) / C
        Ri = ideal_peak(y, win, dnu, k_pk, 2 * R0.max() / win.sum())
        dbi = to_db(Ri, ref)
        print("\nBreite -3 dB: gemessener Peak %.1f um, idealer Reflektor "
              "%.1f um (Zelle %.1f um)"
              % (peak_width(R0, z, zpk, bin_mm),
                 peak_width(Ri, z, zpk, bin_mm), dz_bin * 1e6))
        print("\n### Fuss des Peaks gegen einen idealen Reflektor\n")
        print(excess_table(z[lim], db0[lim], dbi[lim], zpk, a.zmax))

    print("\n### Rauschboden je Fenster -- %s\n" % sub)
    print(floor_table(z[lim], [db0[lim], db1[lim]],
                      ["Original", "danach"], a.zmax))
    print("\n### Wie breit ist der staerkste Ton wirklich? (Original)\n")
    print(width_table(z[lim], db0[lim], zpk, a.zmax))
    print("\n### Rockprofil um den staerksten Ton (Original)\n")
    print(skirt_table(z[lim], db0[lim], zpk, a.zmax))

    spans = [(0.0, a.zmax * 1e3,
              "ganze Achse 0-%.0f mm (Huellkurve ueber alle Bins)"
              % (a.zmax * 1e3)),
             (zpk - 170, zpk + 170, "+-170 mm um den entfernten Ton"),
             (zpk - 18, zpk + 18, "Chip-/Bauteilbereich, +-18 mm"),
             (zpk - 1.5, zpk + 1.5, "der Ton selbst, +-1.5 mm")]

    os.makedirs(a.out_dir, exist_ok=True)
    name = a.name or ("%s_ohne_%dton_%s"
                      % (os.path.basename(a.scan).rsplit(".", 1)[0],
                         a.n_tones, a.method))
    out_png = os.path.join(a.out_dir, name + ".png")
    make_plot(z[lim], db0[lim], db1[lim], zpk, label, sub, out_png, spans,
              dbi=None if dbi is None else dbi[lim])
    out_csv = os.path.join(a.out_dir, name + ".csv")
    np.savetxt(out_csv, np.column_stack([z[lim], db0[lim], db1[lim]]),
               delimiter=",", header="distance_m,amplitude_dB,residual_dB",
               comments="")
    print("geschrieben: %s" % out_csv)


if __name__ == "__main__":
    main()
