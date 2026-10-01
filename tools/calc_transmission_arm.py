#!/usr/bin/env python3
"""
Design calculator for adding a TRANSMISSION arm to the OFDR setup.

The question this answers: if I stop measuring the light that comes back
from the DUT and instead collect the light that comes out of its far end,
where does the peak land, does it still fit under the 150 kHz analog
bandwidth of the CoreDAQ front end, how long may the return fiber be, and
how much signal do I gain?

Conventions (identical to process_reflectogram_aux.py):
  - z is REFLECTION convention: z = c*tau/(2*n_g), i.e. geometric
    fiber-equivalent length. A transmission path of optical length L
    therefore shows up at z = L/2, NOT at L.
  - resolution dz = c/(2*n_g*dnu_span), also reflection convention. In
    transmission one cell corresponds to 2*dz of real path.

The three constraints a transmission arm has to satisfy:

  1. Nyquist     z_T  <  z_max = c/(4*n_g*dnu)
  2. Bandwidth   f_beat(z_T) < 150 kHz   (CoreDAQ analog front end)
                 f_beat = 2*n_g*z/lam^2 * dlam/dt
  3. Separation  if reflection and transmission are merged into one scan,
                 |z_T - z_R| must be more than a few resolution cells.

Examples:
    python calc_transmission_arm.py                      # current chip setup
    python calc_transmission_arm.py --speed 30 --return-m 3.0
    python calc_transmission_arm.py --span 1505 1625 --speed 30
"""

import argparse

C = 299_792_458.0
NG = 1.468          # SMF-28 group index at 1550 nm
NG_CHIP = 1.90      # SiN waveguide group index (logs/2026-09-15.md)
STEP_US = 1.0       # CoreDAQ free-run clock step
F_BW = 150e3        # CoreDAQ analog bandwidth (logs/2026-09-17.md)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--span", nargs=2, type=float, default=[1520.0, 1570.0],
                   metavar=("LAM_START", "LAM_STOP"),
                   help="sweep start/stop in nm (default: standard span)")
    p.add_argument("--points", type=int, default=988_536,
                   help="points per channel (default: as measured)")
    p.add_argument("--speed", type=float, default=None,
                   help="sweep speed nm/s; default = points*1us fills the span")
    p.add_argument("--z-refl", type=float, default=1.582,
                   help="current reflection peak (chip facet A) in m")
    p.add_argument("--chip-mm", type=float, default=6.47,
                   help="one-way waveguide path through the chip in mm")
    p.add_argument("--return-m", type=float, default=None,
                   help="total return path (array fiber + patch cord) in m; "
                        "if omitted the tool reports the allowed budget")
    p.add_argument("--array-m", type=float, default=1.024,
                   help="length of one array fiber in m "
                        "(logs/2026-09-15.md: 1023-1026 mm)")
    p.add_argument("--circulator-m", type=float, default=0.0,
                   help="fiber length REMOVED from the forward path if you "
                        "take the circulator out (port1+port2 pigtails). "
                        "0 = circulator stays in.")
    p.add_argument("--dl-aux", type=float, default=4.214,
                   help="aux MZI arm imbalance in m")
    p.add_argument("--il-chip-db", type=float, default=30.0,
                   help="one-way insertion loss of the DUT in dB")
    p.add_argument("--r-target-db", type=float, default=45.0,
                   help="reflectivity of the structure you wanted to see, "
                        "as a positive dB number below incident")
    a = p.parse_args()

    lam0, lam1 = a.span
    span_nm = lam1 - lam0
    lam_c = 0.5 * (lam0 + lam1) * 1e-9
    duration = a.points * STEP_US * 1e-6
    speed = a.speed if a.speed is not None else span_nm / duration

    nu0, nu1 = C / (lam0 * 1e-9), C / (lam1 * 1e-9)
    span_nu = abs(nu0 - nu1)
    dnu = span_nu / a.points

    dz = C / (2 * NG * span_nu)
    z_max = C / (4 * NG * dnu)

    # beat frequency slope: Hz per metre of z (reflection convention)
    k_bw = 2 * NG / lam_c ** 2 * (speed * 1e-9)
    z_bw = F_BW / k_bw

    print("=" * 68)
    print(f"SWEEP   {lam0:.0f}-{lam1:.0f} nm ({span_nm:.0f} nm), "
          f"{a.points:,} pts, {speed:.1f} nm/s, {duration:.3f} s")
    print(f"        span_nu {span_nu/1e12:.3f} THz   dnu {dnu/1e6:.3f} MHz")
    print("=" * 68)

    print("\n1. RESOLUTION -- unchanged by transmission, set by the span alone")
    print(f"   reflection convention (what the scripts print) {dz*1e6:7.2f} um")
    print(f"   real path per cell in TRANSMISSION (x2)        {2*dz*1e6:7.2f} um")
    print(f"   -> transmission is 2x COARSER in physical length, not finer.")

    print("\n2. RANGE")
    print(f"   Nyquist    z_max   {z_max:6.3f} m")
    print(f"   Bandwidth  z_bw    {z_bw:6.3f} m   "
          f"(f_beat = {k_bw/1e3:.1f} kHz per metre of z)")
    lim, who = (z_max, "Nyquist") if z_max < z_bw else (z_bw, "150 kHz bandwidth")
    print(f"   -> BINDING LIMIT  {lim:6.3f} m  ({who})")

    f_aux = k_bw * a.dl_aux / 2
    print(f"   aux MZI sits at z = {a.dl_aux/2:.3f} m -> f_beat {f_aux/1e3:6.1f} kHz"
          f"  {'OK' if f_aux < F_BW else '*** OVER THE LIMIT ***'}")

    f_refl = k_bw * a.z_refl
    print(f"   reflection peak at  {a.z_refl:.3f} m -> f_beat {f_refl/1e3:6.1f} kHz"
          f"  {'OK' if f_refl < F_BW else '*** OVER THE LIMIT ***'}")

    print("\n3. WHERE THE TRANSMISSION PEAK LANDS")
    print("   NOTE: an internal anchor (529.43 / 573.98 mm) is a REFLECTION.")
    print("   It never reaches the detector in transmission, with or without")
    print("   the circulator. Replace it with a through-reference: input fiber")
    print("   straight into output fiber, no chip.")
    chip_eq = a.chip_mm * 1e-3 * NG_CHIP / NG     # chip path, fiber-equivalent
    z_fwd = a.z_refl - a.circulator_m
    if a.circulator_m:
        print(f"   CIRCULATOR REMOVED: forward path shortened by "
              f"{a.circulator_m:.3f} m")
        print(f"     -> saves ~0.8 dB and lowers f_beat, but the laser loses")
        print(f"        its back-reflection protection (does the EXFO have an")
        print(f"        isolator?).")
    print(f"   forward path to facet A                      {z_fwd:7.3f} m")
    print(f"   chip waveguide {a.chip_mm:.2f} mm at n_g {NG_CHIP}   "
          f"-> fiber-equiv  {chip_eq*1e3:7.2f} mm")
    budget = 2 * lim - z_fwd - chip_eq
    print(f"   => allowed RETURN fiber (array fibre + patch): "
          f"<= {budget:6.3f} m")

    if a.return_m is not None:
        z_t = abs(z_fwd + chip_eq + a.return_m) / 2
        f_t = k_bw * z_t
        sep = abs(z_t - a.z_refl)
        print(f"\n   with --return-m {a.return_m:.3f} m:")
        print(f"     transmission peak at z = {z_t:.3f} m"
              f"   (read it as {2*z_t:.3f} m of real path)")
        print(f"     f_beat = {f_t/1e3:.1f} kHz"
              f"   {'OK' if f_t < F_BW else '*** OVER THE LIMIT ***'}")
        print(f"     {'OK' if z_t < z_max else '*** BEYOND NYQUIST ***'}"
              f"  (z_max {z_max:.3f} m)")
        print(f"     separation from the reflection peak: {sep*1e3:.1f} mm"
              f"  = {sep/dz:,.0f} cells"
              f"   {'OK for one merged scan' if sep > 20*dz else '*** TOO CLOSE ***'}")
    else:
        print(f"   array fiber alone = {a.array_m:.3f} m; "
              f"rows below add a patch cord on top:")
        for Lp in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0):
            L = a.array_m + Lp
            z_t = abs(z_fwd + chip_eq + L) / 2
            f_t = k_bw * z_t
            bad = []
            if f_t >= F_BW:
                bad.append("BANDWIDTH")
            if z_t >= z_max:
                bad.append("NYQUIST")
            if z_t < 0.2:
                bad.append("TOO CLOSE TO DC")
            ok = "ok" if not bad else " + ".join(bad)
            print(f"     patch {Lp:4.1f} m (return {L:5.3f} m) -> "
                  f"z_T {z_t:6.3f} m, f_beat {f_t/1e3:6.1f} kHz   [{ok}]")

    print("\n4. SIGNAL BUDGET")
    refl_db = 2 * a.il_chip_db + a.r_target_db
    trans_db = a.il_chip_db
    print(f"   reflection:  2 x {a.il_chip_db:.0f} dB insertion loss "
          f"+ {a.r_target_db:.0f} dB target reflectivity = {refl_db:5.0f} dB")
    print(f"   transmission: 1 x {a.il_chip_db:.0f} dB insertion loss, "
          f"no reflectivity needed      = {trans_db:5.0f} dB")
    print(f"   -> GAIN {refl_db - trans_db:.0f} dB")
    print(f"      plus ~1.6 dB saved from the circulator (two passes),")
    print(f"      minus 3.0 dB if you merge both arms with a 50:50 coupler.")
    print(f"   dynamic range of the instrument: ~68.7 dB "
          f"(logs/2026-09-15.md) -> reflection needs {refl_db:.0f} dB, "
          f"{'IMPOSSIBLE' if refl_db > 68.7 else 'fits'}")
    print()


if __name__ == "__main__":
    main()
