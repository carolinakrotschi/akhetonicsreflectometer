#!/usr/bin/env python3
"""
Sweep-Parameter-Rechner fuer den Freilauf-Modus (schwarzer CoreDAQ).

Beantwortet: "wenn ich N Punkte pro Kanal habe, mit welcher Geschwindigkeit
muss ich sweepen, welchen Span bekomme ich, welche Reichweite, welche
Auflaesung -- und passt mein Aux-MZI da noch rein?"

Zwei von {N, Sweep-Geschwindigkeit, Span} legen alles fest, weil der
Taktschritt fest ist (1 us):

    T        = N * step
    Span     = v * T
    dlambda  = v * step
    dnu      = c * dlambda / lam^2
    z_max    = c / (4 * n_g * dnu)                <- Reichweite (Nyquist)
    dz       = c / (2 * n_g * dnu_span)           <- Ortsaufloesung
    z_max/dz = N/2                                <- immer, ohne Ausnahme

Aux-MZI-Pruefung (das wird gerne vergessen):
  Das Aux-MZI ist selbst ein Signal auf dem Detektor. Es erscheint bei
  z_aux = dL/2 und muss daher INNERHALB der Nyquist-Reichweite liegen.
  Gleichzeitig soll seine Verzoegerung mindestens so gross sein wie die des
  entferntesten Reflektors, den du messen willst:

      2 * z_mess  <=  dL  <  2 * z_max

Beispiele:
    python ofdr_sweep_calc.py --points 1000000 --span 60 --dl 4
    python ofdr_sweep_calc.py --points 1000000 --speed 120 --dl 4
    python ofdr_sweep_calc.py --points 250000 --span 120 --dl 4
    python ofdr_sweep_calc.py --table
"""

import argparse

C = 299_792_458.0
NG = 1.468
LAM = 1565e-9          # Mitte des Sweeps, fuer die Umrechnungen
STEP_US = 1.0          # Taktschritt des schwarzen CoreDAQ


def compute(points, speed_nms=None, span_nm=None, step_us=STEP_US):
    """Zwei von (speed, span) -- eines darf None sein."""
    T = points * step_us * 1e-6                       # s
    if speed_nms is None and span_nm is None:
        raise SystemExit("--speed oder --span angeben")
    if speed_nms is None:
        speed_nms = span_nm / T
    if span_nm is None:
        span_nm = speed_nms * T

    dlam_pm = speed_nms * step_us * 1e-6 * 1e3        # nm/s * s -> nm -> pm
    dnu = C * dlam_pm * 1e-12 / LAM**2                # Hz
    z_max = C / (4 * NG * dnu)                        # m
    dnu_span = C * span_nm * 1e-9 / LAM**2            # Hz
    dz = C / (2 * NG * dnu_span)                      # m
    dnudt = C / LAM**2 * speed_nms * 1e-9             # Hz/s
    f_nyq = 1.0 / (2 * step_us * 1e-6)                # Hz
    return dict(points=points, T=T, speed=speed_nms, span=span_nm,
                dlam_pm=dlam_pm, dnu=dnu, z_max=z_max, dz=dz,
                cells=z_max / dz, dnudt=dnudt, f_nyq=f_nyq)


def aux_check(r, dl_m, z_mess_m):
    tau_aux = NG * dl_m / C
    z_aux = dl_m / 2
    f_aux = tau_aux * r["dnudt"]
    pts_per_fringe = (1.0 / (STEP_US * 1e-6)) / f_aux
    tau_mess = 2 * NG * z_mess_m / C
    return dict(tau_aux=tau_aux, z_aux=z_aux, f_aux=f_aux,
                pts_per_fringe=pts_per_fringe, tau_mess=tau_mess,
                dl_min=2 * z_mess_m, dl_max=2 * r["z_max"])


def report(r, dl_m=None, z_mess_m=2.0):
    print(f"Punkte/Kanal        {r['points']:,}")
    print(f"Taktschritt         {STEP_US} us  ->  Aufnahmedauer {r['T']:.3f} s")
    print(f"Sweep               {r['speed']:.1f} nm/s ueber {r['span']:.1f} nm")
    print(f"Punktabstand        {r['dlam_pm']:.4f} pm  =  {r['dnu']/1e6:.2f} MHz")
    print()
    print(f"  REICHWEITE  z_max = {r['z_max']:.3f} m")
    print(f"  AUFLOESUNG  dz    = {r['dz']*1e6:.2f} um")
    print(f"  Zellen            = {r['cells']:,.0f}   (= N/2 = {r['points']//2:,})")
    print()
    print(f"Abtast-Nyquist im Zeitbereich: {r['f_nyq']/1e3:.0f} kHz")
    print(f"  Schwebung eines Reflektors bei z_max: "
          f"{2*NG*r['z_max']/C*r['dnudt']/1e3:.0f} kHz  (muss = Nyquist sein)")

    if dl_m is not None:
        a = aux_check(r, dl_m, z_mess_m)
        print()
        print(f"--- Aux-MZI mit dL = {dl_m:.2f} m ---")
        print(f"tau_aux             {a['tau_aux']*1e9:.2f} ns")
        print(f"erscheint bei       z = {a['z_aux']:.3f} m "
              f"(muss < z_max = {r['z_max']:.3f} m sein)")
        print(f"Schwebungsfrequenz  {a['f_aux']/1e3:.1f} kHz  "
              f"= {a['pts_per_fringe']:.1f} Punkte pro Fringe")
        print(f"Zielbereich fuer dL bei z_mess = {z_mess_m:.1f} m: "
              f"{a['dl_min']:.2f} m ... {a['dl_max']:.2f} m")
        ok = True
        if a["z_aux"] >= r["z_max"]:
            print("  FEHLER: Aux liegt jenseits Nyquist -- es aliasiert selbst. "
                  "dL kleiner machen oder langsamer sweepen.")
            ok = False
        if a["pts_per_fringe"] < 2:
            print("  FEHLER: unter 2 Punkte pro Aux-Fringe -- Phase nicht "
                  "rekonstruierbar.")
            ok = False
        elif a["pts_per_fringe"] < 4:
            print("  WARNUNG: unter 4 Punkte pro Aux-Fringe -- Unwrapping wird "
                  "empfindlich. Langsamer sweepen oder dL kleiner.")
        if dl_m < a["dl_min"]:
            print(f"  WARNUNG: dL kleiner als 2*z_mess ({a['dl_min']:.2f} m) -- "
                  "Rauschkorrektur fuer weit entfernte Peaks wird schlechter.")
        if ok and a["pts_per_fringe"] >= 4 and dl_m >= a["dl_min"]:
            print("  -> alle Bedingungen erfuellt.")


def default_table():
    rows = [(1_000_000, 120, None), (1_000_000, 60, None),
            (1_000_000, 30, None), (500_000, 120, None),
            (250_000, None, 120), (250_000, 120, None), (250_000, 60, None)]
    print(f"{'N/Kanal':>10} {'v(nm/s)':>8} {'T(s)':>6} {'Span(nm)':>9} "
          f"{'dlam(pm)':>9} {'z_max(m)':>9} {'dz(um)':>8} {'f_aux@dL=4m':>12}")
    for n, v, s in rows:
        r = compute(n, v, s)
        a = aux_check(r, 4.0, 2.0)
        if a["z_aux"] >= r["z_max"] or a["pts_per_fringe"] < 2:
            flag = "  <-- Aux passt NICHT"
        elif a["pts_per_fringe"] < 4:
            flag = "  <-- Aux grenzwertig"
        else:
            flag = ""
        print(f"{n:>10,} {r['speed']:>8.0f} {r['T']:>6.2f} {r['span']:>9.1f} "
              f"{r['dlam_pm']:>9.3f} {r['z_max']:>9.2f} {r['dz']*1e6:>8.1f} "
              f"{a['f_aux']/1e3:>9.0f} kHz{flag}")
    print("\n(f_aux = Schwebungsfrequenz eines Aux-MZI mit 4 m Armdifferenz;")
    print(" Abtast-Nyquist ist 500 kHz, komfortabel sind < 250 kHz.)")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--points", type=int, default=1_000_000,
                   help="Punkte pro Kanal im Puffer")
    p.add_argument("--speed", type=float, default=None, help="nm/s")
    p.add_argument("--span", type=float, default=None, help="nm")
    p.add_argument("--dl", type=float, default=None,
                   help="Armdifferenz des Aux-MZI in m")
    p.add_argument("--z-mess", type=float, default=2.0,
                   help="entferntester Reflektor, den du messen willst (m)")
    p.add_argument("--step-us", type=float, default=STEP_US)
    p.add_argument("--table", action="store_true", help="Uebersichtstabelle")
    a = p.parse_args()
    if a.table:
        default_table()
        return
    if a.speed is None and a.span is None:
        a.speed = 60.0
    r = compute(a.points, a.speed, a.span, a.step_us)
    report(r, a.dl, a.z_mess)


if __name__ == "__main__":
    main()
