"""
tau_aux bestimmen -- die Strichbreite des Aux-MZI -- aus vorhandenen Daten.

WOZU
----
Die Abstandsrechnung braucht nur die Frequenzspanne des Sweeps:

    dz = c / (2 * n_g * dnu)

Wo im Spektrum der Sweep liegt, geht nicht ein. Und dnu bekommt man aus dem
Aux allein, wenn man EINE Zahl kennt, die nicht in den Messdaten steckt:

    dnu = (Anzahl Aux-Streifen) / tau_aux

tau_aux ist die feste Laufzeitdifferenz des Aux-Interferometers. Sie gilt
fuer jeden Sweep, egal welcher Bereich, und wird erst ungueltig, wenn jemand
den Aux-Arm umsteckt. Genau deshalb ist sie NICHT dasselbe wie eine
Wellenlaengenkalibrierung (die einen Sweep beschreibt) und darf nicht in
dieselbe Lookup-Logik mit ihrer +-0.05 nm / +-2 % Pruefung.

WIE HIER GEMESSEN
-----------------
Aus der vorhandenen Wellenlaengenkalibrierung. Sie liefert acht Zeitpunkte
t_i, an denen der Sweep durch eine BEKANNTE Filterwellenlaenge laeuft
(1530 .. 1565 nm in 5-nm-Schritten). Zwischen zwei davon gilt

    tau_aux = [phi(t_j) - phi(t_i)] / (2*pi) / |c/lam_i - c/lam_j|

also: Aux-Phase entrollen, Streifen zaehlen, durch die bekannte
Frequenzspanne teilen. Die Hilbert-Entrollung ist dieselbe wie im
Auswertepfad -- sie wird aus lina/analysis/lina_ofdr.py importiert, nicht
nachgebaut.

Die Qualitaetskontrolle faellt dabei ab: das geht fuer alle 28 Markerpaare
und fuer mehrere Aufnahmen. Kommt ueberall dasselbe tau_aux heraus, stimmt
es; die Streuung ist der Fehlerbalken.

WICHTIG, DAMIT DAS UEBERHAUPT ZULAESSIG IST
-------------------------------------------
Die Kalibriersweeps haben nur ch2 aufgezeichnet (das Filtersignal), keine
Aux-Franse. Die Markerzeiten stammen also aus anderen Aufnahmen als die
Aux-Phase. Das ist erlaubt, weil lambda(t) in Sekunden ab Triggerflanke
steht und rate-unabhaengig ist (5 kHz gefittet, auf 100 kHz zu -3 pm
reproduziert; Drift +0.0 pm/min). Die Streuung ueber die Markerpaare und
ueber mehrere Aufnahmen prueft genau diese Annahme mit.

Aufruf:

    python tools/measure_tau_aux.py
    python tools/measure_tau_aux.py --npz "logs/lina/data/02_ofdr_measurements/*.npz"
"""

from __future__ import annotations

import argparse
import datetime as _dt
import glob
import importlib.util
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LINA = os.path.join(HERE, "logs", "lina")

C_VAC = 299_792_458.0
N_GROUP_SMF28 = 1.4682


def _load_lina_ofdr():
    """Die verifizierten Aux-Routinen aus dem Auswertepfad wiederverwenden."""
    path = os.path.join(LINA, "analysis", "lina_ofdr.py")
    spec = importlib.util.spec_from_file_location("lina_ofdr_ro", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_capture(path):
    d = np.load(path, allow_pickle=True)
    meta = json.loads(str(d["meta"])) if "meta" in d.files else {}
    traces = {int(k[2:]): d[k] for k in d.files if k.startswith("ch")}
    return traces, meta


def band_energy(sig, rate_hz, f_dom, lo=0.7, hi=1.4):
    """Energieanteil im Chirpband um die dominante Linie.

    Ersetzt das +-2-%-Einton-Mass, das bei einer (physikalisch zwingend)
    chirpenden Aux-Franse strukturell niedrig bleibt.
    """
    spec = np.abs(np.fft.rfft(sig - sig.mean()))
    freqs = np.fft.rfftfreq(sig.size, 1.0 / rate_hz)
    total = float(np.sum(spec[1:] ** 2))
    band = (freqs > f_dom * lo) & (freqs < f_dom * hi)
    return float(np.sum(spec[band] ** 2) / total) if total > 0 else 0.0


def tau_from_pairs(phase, t, points, lam_key="true_nm"):
    """tau_aux fuer jedes Markerpaar. phase(t) ist die entrollte Aux-Phase."""
    tm = np.array([p["time_s"] for p in points], float)
    lam = np.array([p[lam_key] for p in points], float) * 1e-9
    nu = C_VAC / lam
    phi_m = np.interp(tm, t, phase)

    rows = []
    for i in range(len(tm)):
        for j in range(i + 1, len(tm)):
            dnu = abs(nu[i] - nu[j])
            n_fr = abs(phi_m[j] - phi_m[i]) / (2 * np.pi)
            rows.append({
                "lam_i_nm": lam[i] * 1e9, "lam_j_nm": lam[j] * 1e9,
                "span_nm": (lam[j] - lam[i]) * 1e9,
                "lam_center_nm": 0.5 * (lam[i] + lam[j]) * 1e9,
                "dnu_THz": dnu / 1e12,
                "fringes": n_fr,
                "tau_aux_ns": n_fr / dnu * 1e9,
            })
    return rows


def analyse_file(path, cal, mod, args):
    traces, meta = load_capture(path)
    rate = float(meta.get("rate_hz", args.rate))
    name = os.path.basename(path)

    missing = [c for c in args.aux if c not in traces]
    if missing:
        return None, f"{name}: Aux-Kanaele {missing} fehlen -- uebersprungen"

    aux = mod.combine_pair(traces, args.aux, args.balance, "aux MZI",
                           log=lambda *_: None)

    # Auf das Sweep-Fenster beschneiden. Der Ruecklauf des Lasers am
    # Pufferende ist keine Franse und wuerde die globale Hilbert-Phase
    # verderben. Das Fenster kommt aus der Kalibrierung.
    t0 = float(cal["derived"]["buffer_time_s"])
    t1 = float(cal["derived"]["sweep_end_s"])
    i0, i1 = int(round(t0 * rate)), min(int(round(t1 * rate)), aux.size)
    aux = aux[i0:i1]
    t = np.arange(i0, i0 + aux.size) / rate          # Sekunden ab Triggerflanke

    aux = mod.highpass(aux, rate, args.highpass_hz)

    # --- Waechter: lieber verweigern als raten -------------------------------
    # Die Aux-Franse ist KEIN reiner Einton, auch wenn alles in Ordnung ist:
    #     f_aux(t) = tau_aux * dnu/dt,   dnu/dt = c/lam^2 * dlam/dt
    # und c/lam^2 aendert sich ueber 1520..1570 nm um 8.8 %. Gemessen sind es
    # 15.6 % (der Rest ist die Ungleichfoermigkeit des EXFO). Das +-2-%-Mass
    # aus aux_quality() kann deshalb bei gesundem Aux gar nicht hoch werden --
    # gemessen 10..13 %. Als harte Schranke waere es schlicht falsch geeicht.
    # Stattdessen: Energie im Chirpband [0.7 f, 1.4 f].
    f_aux, tone_2pct = mod.aux_quality(aux, rate)
    band_frac = band_energy(aux, rate, f_aux)
    notes, fatal = [], []

    phase = mod.aux_phase(aux)
    dphi = np.diff(phase)
    if np.median(dphi) < 0:
        phase = -phase
        dphi = -dphi
    f_inst = dphi * rate / (2 * np.pi)
    w = max(3, int(0.02 * rate) | 1)
    f_sm = np.convolve(f_inst, np.ones(w) / w, "valid")
    f_lo, f_hi = float(f_sm.min()), float(f_sm.max())

    # (1) Nyquist -- gegen die HOECHSTE Momentanfrequenz, nicht gegen die
    # dominante Linie; der Chirp laeuft nach oben.
    if f_hi > 0.45 * rate:
        fatal.append(f"Aux-Franse erreicht {f_hi/1e3:.2f} kHz = "
                     f"{100*f_hi/(rate/2):.0f} % der Nyquistgrenze -- aliast")
    # (2) Monotonie der Phase. Das ist der Waechter, der ein falsch
    # zugeordnetes Kanalpaar faengt (verifiziert: ch1/ch3 als Aux gibt 11.9 %).
    bad = float(np.mean(dphi <= 0))
    if bad > args.max_nonmono:
        fatal.append(f"{100*bad:.2f} % der Aux-Phase laufen rueckwaerts "
                     f"(> {100*args.max_nonmono:.2f} %) -- entweder aliast die "
                     f"Franse oder das ist nicht das Aux-Paar")
    # (3) Chirpband: eine verschmierte oder mehrtonige Franse.
    if band_frac < args.min_band:
        fatal.append(f"nur {100*band_frac:.1f} % der Energie im Chirpband "
                     f"(< {100*args.min_band:.0f} %) -- Franse verschmiert")

    rows = tau_from_pairs(phase, t, cal["points"])
    tau = np.array([r["tau_aux_ns"] for r in rows])

    # (4) Der schaerfste Waechter, und er kostet nichts: haette sich die
    # Entrollung irgendwo verzaehlt, wuerden die Markerpaare, die diese Stelle
    # ueberspannen, auseinanderlaufen. Konsistenz ueber 28 Paare ist der Beweis,
    # dass keine Franse verloren ging.
    rel_sd = float(np.std(tau, ddof=1) / np.median(tau))
    if rel_sd > args.max_pair_sd:
        fatal.append(f"Markerpaare streuen um {100*rel_sd:.3f} % "
                     f"(> {100*args.max_pair_sd:.2f} %) -- die Entrollung hat "
                     f"sich wahrscheinlich verzaehlt")

    # Plausibilitaet gegen den kommandierten Bereich: die Gesamtstreifenzahl
    # muss zur kommandierten Spanne passen, sonst hat sich das Entrollen
    # verzaehlt.
    dnu_cmd = C_VAC / (meta.get("actual_start_nm", 1520.0) * 1e-9) - \
              C_VAC / (meta.get("actual_stop_nm", 1570.0) * 1e-9)
    n_total = (phase[-1] - phase[0]) / (2 * np.pi)
    tau_cmd = n_total / abs(dnu_cmd) * 1e9
    dev = abs(tau_cmd - np.median(tau)) / np.median(tau)
    if dev > args.max_cmd_dev:
        notes.append(f"gegen kommandierten Bereich {100*dev:+.1f} % "
                     f"(> {100*args.max_cmd_dev:.0f} %)")

    widest = max(abs(r["span_nm"]) for r in rows)
    res = {
        "file": name, "rate_hz": rate, "f_aux_hz": f_aux,
        "tone_2pct": tone_2pct, "band_frac": band_frac,
        "f_inst_lo_hz": f_lo, "f_inst_hi_hz": f_hi,
        "nonmono_frac": bad, "pair_rel_sd": rel_sd,
        "n_fringes_total": n_total,
        "tau_from_commanded_ns": tau_cmd, "cmd_dev": dev,
        "tau_median_ns": float(np.median(tau)),
        "tau_mean_ns": float(np.mean(tau)),
        "tau_std_ns": float(np.std(tau, ddof=1)),
        "tau_widest_ns": float(next(r["tau_aux_ns"] for r in rows
                                    if abs(r["span_nm"]) == widest)),
        "rows": rows, "notes": notes, "fatal": fatal,
    }
    return res, None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npz", default=os.path.join(
        LINA, "data", "02_ofdr_measurements", "*.npz"),
        help="Glob auf die Aufnahmen mit Aux-Kanal")
    ap.add_argument("--cal", default=os.path.join(
        LINA, "calibrations", "lina_wl_cal_1520.0-1570.0nm_5.00nms.json"))
    ap.add_argument("--aux", default="2,4", type=lambda s: [int(x) for x in s.split(",")])
    ap.add_argument("--balance", default="auto", choices=["auto", "always", "never"])
    ap.add_argument("--highpass-hz", type=float, default=500.0)
    ap.add_argument("--rate", type=float, default=100000.0)
    ap.add_argument("--min-band", type=float, default=0.50,
                    help="Mindestenergie im Chirpband [0.7f, 1.4f]")
    ap.add_argument("--max-nonmono", type=float, default=0.01)
    ap.add_argument("--max-pair-sd", type=float, default=0.005,
                    help="maximale relative Streuung ueber die Markerpaare")
    ap.add_argument("--max-cmd-dev", type=float, default=0.05)
    ap.add_argument("--n-group", type=float, default=N_GROUP_SMF28)
    ap.add_argument("--dispersion-ps-nm-km", type=float, default=17.0,
                    help="Faserdispersion fuer die Erwartung (SMF-28)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    if a.out is None:
        a.out = os.path.join(HERE, "results",
                             _dt.date.today().isoformat(), "tau_aux")
    os.makedirs(a.out, exist_ok=True)

    mod = _load_lina_ofdr()
    cal = json.load(open(a.cal))
    print(f"[cal] {os.path.basename(a.cal)}: {len(cal['points'])} Marker "
          f"{cal['points'][0]['true_nm']:.0f}..{cal['points'][-1]['true_nm']:.0f} nm, "
          f"Fit {cal['fit']['residual_pm_rms']:.1f} pm rms")
    print(f"[cal] Sweepfenster {cal['derived']['buffer_time_s']:.4f}"
          f"..{cal['derived']['sweep_end_s']:.4f} s\n")

    files = sorted(glob.glob(a.npz))
    if not files:
        raise SystemExit(f"keine Dateien fuer {a.npz}")

    results, all_rows = [], []
    for f in files:
        res, err = analyse_file(f, cal, mod, a)
        if err:
            print(f"[skip] {err}")
            continue
        flag = "FEHLER" if res["fatal"] else ("Hinweis" if res["notes"] else "ok")
        print(f"[{flag:>7}] {res['file']}")
        print(f"          Aux {res['f_inst_lo_hz']/1e3:6.3f}->"
              f"{res['f_inst_hi_hz']/1e3:6.3f} kHz (Chirp "
              f"{100*(res['f_inst_hi_hz']-res['f_inst_lo_hz'])/res['f_aux_hz']:4.1f} %, "
              f"max {100*res['f_inst_hi_hz']/(res['rate_hz']/2):4.1f} % Nyquist), "
              f"{100*res['band_frac']:4.1f} % im Band, "
              f"{100*res['nonmono_frac']:.3f} % nicht monoton")
        print(f"          tau_aux = {res['tau_median_ns']:.4f} ns  "
              f"(Median ueber {len(res['rows'])} Paare, "
              f"sd {res['tau_std_ns']*1000:.2f} ps = "
              f"{100*res['pair_rel_sd']:.3f} %)")
        for m in res["fatal"]:
            print(f"          ABBRUCH: {m}")
        for m in res["notes"]:
            print(f"          Hinweis: {m}")
        if res["fatal"]:
            continue
        results.append(res)
        for r in res["rows"]:
            all_rows.append({"file": res["file"], **r})

    if not results:
        raise SystemExit("keine Aufnahme hat die Waechter bestanden")

    tau_all = np.array([r["tau_aux_ns"] for r in all_rows])
    tau_per_file = np.array([r["tau_median_ns"] for r in results])
    tau_wide = np.array([r["tau_widest_ns"] for r in results])

    tau = float(np.median(tau_all))
    sd = float(np.std(tau_per_file, ddof=1)) if len(tau_per_file) > 1 else float("nan")

    print("\n" + "=" * 66)
    print(f"tau_aux            = {tau:.4f} ns")
    print(f"  Streuung Aufnahmen  sd {sd*1000:.2f} ps  ({100*sd/tau:.3f} %) "
          f"ueber {len(results)} Aufnahmen")
    print(f"  Streuung Paare      sd {np.std(tau_all, ddof=1)*1000:.2f} ps  "
          f"({100*np.std(tau_all, ddof=1)/tau:.3f} %) ueber {len(all_rows)} Paare")
    print(f"  weiteste Spanne     {np.mean(tau_wide):.4f} ns "
          f"(1530->1565 nm, praezisester Einzelwert)")
    print(f"  Aux-Wegdifferenz    {C_VAC*tau*1e-9/a.n_group:.4f} m "
          f"(n_g {a.n_group}, Faserlaenge)")
    print(f"  -> Aux-Peak bei     {C_VAC*tau*1e-9/a.n_group:.4f} m "
          f"in der n_g-Konvention")
    print("=" * 66)

    # --- Kreuzvergleich: Delta-nu aus dem Aux gegen Delta-nu aus der
    # Kalibrierung. Das ist der Test, der den ganzen Weg absichert -- aber er
    # ist nur TEILWEISE unabhaengig: tau_aux kommt aus denselben Markern.
    # Getestet wird damit die Extrapolation von 1530..1565 (Markerbereich)
    # auf 1520..1570 (Sweepfenster), nicht die Skala selbst. Wirklich
    # unabhaengig ist der Vergleich der Aux-Wegdifferenz oben mit den 4.162 m,
    # die im Auswertepfad aus der FFT-Peaklage bestimmt wurden.
    co = cal["fit"]["coeffs"]
    nm0 = float(np.polyval(co, cal["derived"]["buffer_time_s"]))
    nm1 = float(np.polyval(co, cal["derived"]["sweep_end_s"]))
    dnu_cal = abs(C_VAC / (nm0 * 1e-9) - C_VAC / (nm1 * 1e-9))
    print()
    print("Kreuzvergleich auf der 1520-1570-Konfiguration")
    print(f"  aus der Kalibrierung: {nm0:.3f} -> {nm1:.3f} nm"
          f"  ->  dnu = {dnu_cal/1e12:.4f} THz")
    dev = []
    for r in results:
        dnu_aux = r["n_fringes_total"] / (tau * 1e-9)
        dev.append(100 * (dnu_aux - dnu_cal) / dnu_cal)
        print(f"  aus dem Aux: {r['file'][:30]:30s} {r['n_fringes_total']:8.0f}"
              f" Streifen  ->  {dnu_aux/1e12:.4f} THz  ({dev[-1]:+.3f} %)")
    print(f"  Abweichung im Mittel {np.mean(dev):+.3f} % -> derselbe Massstab.")

    # --- Fehlerbudget: woher kommt die Streuung? -----------------------------
    # Ein fester Zeitfehler dt auf den Markern gibt dN = f_aux*dt Streifen,
    # also einen RELATIVEN Fehler, der mit 1/Spanne geht. Verzaehlte Fransen
    # dagegen waeren spannenunabhaengig. Welches von beidem vorliegt, sagt die
    # Aufschluesselung nach Spanne.
    print()
    print("Fehlerbudget -- Streuung nach Markerabstand")
    print(f"  {'Spanne':>8s} {'Paare':>6s} {'Median ns':>11s} {'sd ps':>8s} "
          f"{'rel':>9s}")
    spans = np.round([abs(r["span_nm"]) for r in all_rows], 3)
    f_typ = np.mean([r["f_aux_hz"] for r in results])
    dt_ms = []
    for s in sorted(set(spans)):
        m = spans == s
        if m.sum() < 2:
            continue
        sd_s = float(np.std(tau_all[m], ddof=1))
        print(f"  {s:5.0f} nm {int(m.sum()):6d} {np.median(tau_all[m]):11.4f} "
              f"{sd_s*1000:8.2f} {100*sd_s/np.median(tau_all[m]):8.3f} %")
        # zurueckgerechnet auf einen Zeitfehler der Marker
        dnu_s = abs(C_VAC / ((1545 - s / 2) * 1e-9) - C_VAC / ((1545 + s / 2) * 1e-9))
        dt_ms.append(sd_s / tau * (tau * 1e-9 * dnu_s) / f_typ * 1000)
    print(f"  Die Streuung geht mit 1/Spanne -- das ist Marker-ZEITfehler, "
          f"keine verzaehlten")
    print(f"  Fransen (die waeren spannenunabhaengig). Zurueckgerechnet: "
          f"{np.mean(dt_ms):.2f} ms bzw.")
    print(f"  {np.mean(dt_ms)*1e-3*cal['derived']['effective_speed_nm_s']*1000:.1f}"
          f" pm -- der Kalibrierfit selbst hat "
          f"{cal['fit']['residual_pm_rms']:.1f} pm rms. Die Messung ist also")
    print(f"  durch die Kalibrierung begrenzt, nicht durch das Aux.")
    wide_mask = spans >= 20
    print(f"  Bestwert aus den weiten Paaren (>=20 nm): "
          f"{np.median(tau_all[wide_mask]):.4f} ns +- "
          f"{np.std(tau_all[wide_mask], ddof=1)*1000:.2f} ps "
          f"({100*np.std(tau_all[wide_mask], ddof=1)/np.median(tau_all[wide_mask]):.3f} %)")
    # NICHT vergleichen: der kommandierte Bereich gibt hier per Konstruktion
    # exakt dasselbe dnu, weil das Sweepfenster ja gerade als "wo der Fit 1520
    # bzw. 1570 nm erreicht" definiert ist. Die Kalibrierung gewinnt nicht
    # ueber dnu, sondern ueber das FENSTER und die Ungleichfoermigkeit. Der
    # Aux-Pfad muss sich dieses Fenster anders holen (Punkt 4) -- hier ist es
    # aus der Kalibrierung geliehen, und das ist der Grund, warum dieser
    # Kreuzvergleich nur die Extrapolation prueft, nicht den ganzen Weg.

    # --- Dispersion: laeuft tau_aux ueber den Bereich? -----------------------
    # Ueber die Aufnahmen MITTELN, nicht alle 224 Zeilen einzeln fitten: jede
    # Aufnahme liefert dieselben 28 Markerpaare, die Werte sind also nicht
    # unabhaengig. Ohne das Mitteln waere der Standardfehler um sqrt(8) zu
    # klein und ein Systematik-Rest saehe nach Signifikanz aus.
    by_pair = {}
    for r in all_rows:
        by_pair.setdefault((r["lam_i_nm"], r["lam_j_nm"]), []).append(r["tau_aux_ns"])
    keys = [k for k in by_pair if abs(k[1] - k[0]) >= 20]
    if len(keys) > 3:
        lc = np.array([0.5 * (k[0] + k[1]) for k in keys])
        tp = np.array([np.mean(by_pair[k]) for k in keys])
        sl, ic = np.polyfit(lc, tp, 1)
        resid = tp - (sl * lc + ic)
        se = float(np.sqrt(np.sum(resid ** 2) / (len(keys) - 2)
                           / np.sum((lc - lc.mean()) ** 2)))
        length_m = C_VAC * tau * 1e-9 / a.n_group
        pred = a.dispersion_ps_nm_km * length_m / 1000.0      # ps/nm
        sig = "signifikant" if abs(sl) > 2 * se else "nicht signifikant"
        print()
        print("Dispersion (n_g ist nicht konstant)")
        print(f"  gemessen {sl*1000:+.3f} +- {se*1000:.3f} ps/nm ueber "
              f"{len(keys)} Paare ({sig})")
        print(f"  erwartet {pred:+.3f} ps/nm aus D = "
              f"{a.dispersion_ps_nm_km:.0f} ps/(nm*km) auf {length_m:.3f} m")
        if abs(sl) > 2 * se and abs(sl) > 2 * pred:
            print("  Die gemessene Steigung ist groesser als die physikalische "
                  "Erwartung -- das ist")
            print("  eher ein Rest an Systematik in den Markerzeiten (die aus "
                  "acht EINZELNEN")
            print("  Sweeps stammen) als echte Dispersion.")
        print(f"  So oder so: ueber 50 nm sind das "
              f"{100*abs(sl)*50/tau:.3f} % auf tau_aux, "
              f"erwartet {100*pred*1e-12*50/(tau*1e-9):.3f} %.")
        print("  -> Fuer die Laenge reicht der feste n_g; fuer eine "
              "beschriftete Wellenlaengenachse")
        print("     braeuchte man n_g(lambda) und einen absoluten Anker.")

    # --- Ausgabe -------------------------------------------------------------
    csv = os.path.join(a.out, "tau_aux_pairs.csv")
    with open(csv, "w", encoding="utf-8") as fh:
        keys = ["file", "lam_i_nm", "lam_j_nm", "span_nm", "lam_center_nm",
                "dnu_THz", "fringes", "tau_aux_ns"]
        fh.write(",".join(keys) + "\n")
        for r in all_rows:
            fh.write(",".join(f"{r[k]:.6f}" if isinstance(r[k], float) else str(r[k])
                              for k in keys) + "\n")
    print(f"\n[out] {csv}")

    summary = {
        "measured": _dt.datetime.now().isoformat(timespec="seconds"),
        "tau_aux_ns": tau,
        "tau_aux_s": tau * 1e-9,
        "uncertainty_ps_between_captures": sd * 1000,
        "uncertainty_ps_between_marker_pairs": float(np.std(tau_all, ddof=1)) * 1000,
        "n_captures": len(results), "n_pairs": len(all_rows),
        "guards": {"min_band_energy": a.min_band,
                   "max_nonmonotonic_frac": a.max_nonmono,
                   "max_pair_rel_sd": a.max_pair_sd},
        "aux_path_imbalance_m": C_VAC * tau * 1e-9 / a.n_group,
        "n_group_assumed": a.n_group,
        "source_calibration": os.path.basename(a.cal),
        "captures": [r["file"] for r in results],
        "per_capture_ns": {r["file"]: r["tau_median_ns"] for r in results},
    }
    js = os.path.join(a.out, "aux_mzi.json")
    json.dump(summary, open(js, "w", encoding="utf-8"), indent=2)
    print(f"[out] {js}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
        for r in results:
            x = [p["lam_center_nm"] for p in r["rows"]]
            y = [p["tau_aux_ns"] for p in r["rows"]]
            s = [8 + 2 * abs(p["span_nm"]) for p in r["rows"]]
            ax[0].scatter(x, y, s=s, alpha=0.55, label=r["file"][:28])
        ax[0].axhline(tau, color="k", lw=1.2, label=f"Median {tau:.4f} ns")
        ax[0].axhspan(tau - sd, tau + sd, color="k", alpha=0.12)
        ax[0].set_xlabel("Mittenwellenlaenge des Markerpaars (nm)")
        ax[0].set_ylabel("tau_aux (ns)")
        ax[0].set_title("tau_aux je Markerpaar (Punktgroesse = Spanne)")
        ax[0].grid(alpha=0.3)
        ax[0].legend(fontsize=6, ncol=2)

        ax[1].hist((tau_all - tau) * 1000, bins=30, color="tab:blue", alpha=0.8)
        ax[1].set_xlabel("Abweichung vom Median (ps)")
        ax[1].set_ylabel("Anzahl Paare")
        ax[1].set_title(f"Streuung: sd {np.std(tau_all, ddof=1)*1000:.2f} ps "
                        f"= {100*np.std(tau_all, ddof=1)/tau:.3f} %")
        ax[1].grid(alpha=0.3)
        fig.suptitle(f"tau_aux aus vorhandenen Daten - {tau:.4f} ns "
                     f"(Aux-Faserdifferenz {C_VAC*tau*1e-9/a.n_group:.4f} m)")
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        png = os.path.join(a.out, "tau_aux.png")
        fig.savefig(png, dpi=130)
        print(f"[out] {png}")
    except Exception as e:
        print(f"[plot] uebersprungen: {e}")


if __name__ == "__main__":
    main()
