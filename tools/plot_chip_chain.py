"""Erklaerplot: was am HHI-Chip gemessen ist und was Modell ist.

4 Reihen:
  1  Schema der ganzen Kette (gezeichnet, keine Daten)
  2  das gemessene Reflektogramm dazu, x in mm
  3  Zoom auf den Chip: NUR die zwei gemessenen Peaks, darueber das Modell
  4  Detail: wie weit Facette B von der Erwartung abweicht

Aufruf:  python tools/plot_chip_chain.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

NF, NG = 1.468, 3.48033      # Gruppenindex Faser / Chip
STR = NG / NF                # 1 um on-chip -> 2.3704 um im Reflektogramm
K = 1.0 / STR
SSC = 1200.0                 # um, Laenge der SSC-Zelle im GDS (ANNAHME)
RES = 16.59                  # um, eine Aufloesungszelle

S = [dict(name="top", ch="b27/b28", Ldes=640.1266, col="navy",
          csv="results/2026-09-11/measured_loop_top_reflectogram.csv",
          ref=573.9284, fa=1670.2266, fb=1677.4262, far=2773.6249),
     dict(name="bottom", ch="b0/b1", Ldes=556.5494, col="darkorange",
          csv="results/2026-09-11/measured_loop_bottom_reflectogram.csv",
          ref=573.9747, fa=1668.3752, fb=1675.3923, far=2770.5062)]
for v in S:
    v["fiber1"] = v["fa"] - v["ref"]
    v["fiber2"] = v["far"] - v["fb"]
    v["dz"] = v["fb"] - v["fa"]                  # mm, gemessener Chip-Pfad
    v["s_m"] = v["dz"] * K * 1e3                 # um on-chip, gemessen
    v["L_m"] = v["s_m"] - 2 * SSC                # um, Schleife (Modell)
    v["s_e"] = 2 * SSC + v["Ldes"]               # um on-chip, erwartet
    v["fb_e"] = v["fa"] + v["s_e"] * 1e-3 * STR  # mm, erwartete Facette B
    v["err"] = v["L_m"] - v["Ldes"]
    v["err_z"] = (v["fb"] - v["fb_e"]) * 1e3     # um faseraequivalent


def load(p):
    d = np.genfromtxt(p, delimiter=",", names=True)
    return d["distance_m"] * 1e3, d["amplitude_dB"]      # mm, dB


fig = plt.figure(figsize=(14, 17))
gs = fig.add_gridspec(4, 1, height_ratios=[.62, .95, 1.35, .95], hspace=.55)
v0 = S[0]

# ---------------------------------------------------------------- Reihe 1
ax = fig.add_subplot(gs[0])
ax.axis("off")
ax.set_xlim(450, 2900)
ax.set_ylim(-1.15, 1.5)
segs = [(v0["ref"], v0["fa"], "#9ecae1",
         "FIBER 1\n%.1f mm\n(the ~1 m patchcord)" % v0["fiber1"]),
        (v0["fa"], v0["fb"], "#e15759", "CHIP\n%.3f mm" % v0["dz"]),
        (v0["fb"], v0["far"], "#9ecae1",
         "FIBER 2\n%.1f mm\n(the ~1 m patchcord)" % v0["fiber2"])]
for x0, x1, c, lab in segs:
    ax.add_patch(Rectangle((x0, -.28), max(x1 - x0, 8), .56, fc=c, ec="k", lw=.8))
    ax.annotate(lab, ((x0 + x1) / 2, .42), fontsize=9.5, ha="center",
                va="bottom", fontweight="bold")
for x, lab in [(v0["ref"], "internal\nreflection\n%.2f mm" % v0["ref"]),
               ((v0["fa"] + v0["fb"]) / 2,
                "facet A %.2f mm  and  facet B %.2f mm\n(only 7.2 mm apart -- see row 3)"
                % (v0["fa"], v0["fb"])),
               (v0["far"], "open fiber end\n%.2f mm" % v0["far"])]:
    ax.plot([x], [-.30], marker="v", ms=9, color="crimson", clip_on=False)
    ax.annotate(lab, (x, -.44), fontsize=8.5, ha="center", va="top", color="crimson")
ax.annotate("light path:  reflectometer -> fiber 1 -> through the chip -> fiber 2 -> "
            "reflects off the open fiber end -> all the way back",
            (1675, 1.2), fontsize=10, ha="center", style="italic")
ax.set_title("1)  WHAT IS PHYSICALLY THERE  (top loop b27/b28 -- the chip is only 7 mm "
             "of a 2.2 m path)", fontsize=12, loc="left")

# ---------------------------------------------------------------- Reihe 2
ax = fig.add_subplot(gs[1])
z, db = load(v0["csv"])
m = (z > 450) & (z < 2900)
ax.plot(z[m], db[m], lw=.45, color=v0["col"])
for x0, x1, c, lab in segs:
    ax.axvspan(x0, x1, color=c, alpha=.35, zorder=0)
for x, lab in [(v0["ref"], "internal reflection"),
               ((v0["fa"]+v0["fb"])/2, "facets A + B  =  the chip"),
               (v0["far"], "open fiber end")]:
    ax.axvline(x, color="crimson", ls=":", lw=.9)
    ax.annotate(lab, (x, 8), fontsize=8.5, ha="center", va="bottom", color="crimson")
ax.set_xlim(450, 2900)
ax.set_ylim(-70, 26)
ax.set_xlabel("distance  [mm]")
ax.set_ylabel("amplitude  [dB]")
ax.grid(alpha=.3)
ax.set_title("2)  THE SAME THING AS MEASURED -- each peak is one reflector, the coloured "
             "stretches are the two fibres", fontsize=12, loc="left")

# ---------------------------------------------------------------- Reihe 3
ax = fig.add_subplot(gs[2])
for v in S:
    z, db = load(v["csv"])
    m = (z > v["fa"] - 1.4) & (z < v["fa"] + 9.5)
    ax.plot(z[m] - v["fa"], db[m], lw=1.1, color=v["col"],
            label="%s loop (%s)" % (v["name"], v["ch"]), zorder=3)
    ax.axvline(v["dz"], color=v["col"], ls="-", lw=1.6, zorder=2)
ax.axvline(0, color="crimson", ls="-", lw=1.3, zorder=2)
ax.annotate("MEASURED PEAK\nfacet A  (x=0 by definition)", (0, -45), fontsize=9,
            ha="center", va="bottom", color="crimson",
            bbox=dict(fc="w", ec="crimson", lw=.8, alpha=.92))
for v, dy in zip(S, (-45, -36)):
    ax.annotate("MEASURED PEAK\nfacet B at +%.4f mm" % v["dz"], (v["dz"], dy),
                fontsize=9, ha="center", va="bottom", color=v["col"],
                bbox=dict(fc="w", ec=v["col"], lw=.8, alpha=.92))
ax.annotate("these two peaks are EVERYTHING the instrument knows about the chip",
            (4.0, -50.5), fontsize=10.5, ha="center", style="italic")

for v, ytop in zip(S, (13.2, 5.4)):
    h = 3.0
    parts = [(0, SSC, "#d9d9d9", "SSC  1200 um"),
             (SSC, SSC + v["L_m"], "gold",
              "THE LOOP  %.1f um   (design %.1f)" % (v["L_m"], v["Ldes"])),
             (SSC + v["L_m"], v["s_m"], "#d9d9d9", "SSC  1200 um")]
    for s0, s1, c, lab in parts:
        ax.add_patch(Rectangle((s0 * 1e-3 * STR, ytop), (s1 - s0) * 1e-3 * STR, h,
                               fc=c, ec="k", lw=.7, alpha=.85, zorder=4))
        ax.annotate(lab, ((s0 + s1) / 2 * 1e-3 * STR, ytop + h / 2), fontsize=8.5,
                    ha="center", va="center", zorder=5,
                    fontweight="bold" if c == "gold" else "normal")
    ax.annotate("MODEL, %s loop\nNOT measured" % v["name"], (-3.3, ytop + h / 2),
                fontsize=8.5, ha="left", va="center", style="italic", color="dimgray")
ax.set_xlim(-3.4, 9.5)
ax.set_ylim(-54, 19)
ax.set_xlabel("distance behind facet A  [mm]")
ax.set_ylabel("amplitude  [dB]")
ax.grid(alpha=.25)
ax.legend(fontsize=9, loc="lower right", framealpha=.95)
ax.set_title("3)  ZOOM ON THE CHIP.  Measured: 2 peaks.  The bars on top are the MODEL -- "
             "the loop edges are NOT in the data.\n"
             "     Why: the loop passes smoothly into the SSC, so its ends do not reflect. "
             "And dips in the curve are window\n"
             "     sidelobes plus noise, not structure -- so the band cannot be pinned to a dip.",
             fontsize=11.5, loc="left")

# ---------------------------------------------------------------- Reihe 4
ax = fig.add_subplot(gs[3])
for v in S:
    z, db = load(v["csv"])
    m = np.abs(z - v["fb"]) < .075
    x = (z[m] - v["fb_e"]) * 1e3
    ax.plot(x, db[m] - db[m].max(), lw=1.8, color=v["col"],
            label="%s loop: peak is %+.1f um from expected" % (v["name"], v["err_z"]))
    ax.plot([v["err_z"]], [0], marker="o", ms=9, color=v["col"], zorder=5)
ax.axvline(0, color="green", ls="--", lw=2)
ax.annotate("WHERE THE DESIGN SAYS facet B should be\n"
            "(GDS length converted with n_g = 3.48033)",
            (-40, -13.2), fontsize=9.5, ha="center", color="green",
            bbox=dict(fc="w", ec="green", lw=.9, alpha=.92))
ax.add_patch(Rectangle((-RES / 2, -1.7), RES, 1.7, fc="k", alpha=.13, ec="k", lw=.7))
ax.annotate("grey box = one resolution cell = %.1f um\n"
            "(the instrument cannot resolve finer, so\nboth peaks are effectively spot on)" % RES,
            (40, -3.4), fontsize=9, ha="center")
for v, dy in zip(S, (-8.6, -11.0)):
    ax.annotate("%s loop:  %+.1f um here  =  %+.2f um on the chip"
                % (v["name"], v["err_z"], v["err"]), (-40, dy),
                fontsize=9.5, ha="center", color=v["col"], fontweight="bold")
ax.set_xlim(-70, 70)
ax.set_ylim(-16, 2.5)
ax.set_xlabel("distance from the expected facet-B position  [um]      "
              "(negative = chip measured shorter than the design)")
ax.set_ylabel("amplitude relative to own peak  [dB]")
ax.grid(alpha=.3)
ax.legend(fontsize=9.5, loc="lower right")
ax.set_title("4)  HOW CLOSE IS IT?  Both peaks land within half a resolution cell of the "
             "design, on opposite sides --\n     which is what you get when one n_g is fitted "
             "to both loops instead of tuned to one.", fontsize=11.5, loc="left")

fig.savefig("results/2026-09-11/measured_loops_explained.png", dpi=120, bbox_inches="tight")
print("results/2026-09-11/measured_loops_explained.png")
for v in S:
    print("%-7s fiber1=%.3f mm  chip=%.4f mm (%.1f um on-chip)  fiber2=%.3f mm  "
          "loop=%.2f um (design %.2f, err %+.2f)  facetB err %+.1f um = %.2f cells"
          % (v["name"], v["fiber1"], v["dz"], v["s_m"], v["fiber2"],
             v["L_m"], v["Ldes"], v["err"], v["err_z"], v["err_z"] / RES))
