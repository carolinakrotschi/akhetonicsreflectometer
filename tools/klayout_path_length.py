"""Exakte Mittellinienlaenge eines Waveguide-Routings in KLayout.

VERWENDUNG (KLayout 0.30.x)
   1. Layout oeffnen, mit dem Select-Tool das schwarze Routing selektieren
      (klicken, weitere mit Shift+Klick dazu, oder Rahmen ziehen).
      Die gruenen Rechtecke / Edge-Coupler NICHT mitselektieren.
   2. Macros -> Macro Development -> neues Python-Makro -> diese Datei
      oeffnen -> Run (F5).
   3. Ausgabe im Konsolenfenster: Mittellinienlaenge + Breite in um.

STANDALONE (Flaeche/Umfang schon bekannt, z.B. aus einem DRC-Skript):
   python tools/klayout_path_length.py --area 1234.5 --perimeter 2345.6

METHODE
   Ein Band konstanter Breite w mit Mittellinienlaenge L und stumpf
   abgeschnittenen Enden erfuellt exakt:

       A = L * w                (Flaeche)
       P = 2*L + 2*N*w          (Umfang, N = Anzahl getrennter Baender)

   Das gilt auch fuer Boegen: Aussenbogen (R+w/2)*theta plus Innenbogen
   (R-w/2)*theta = 2*R*theta = 2*L_mittellinie -- der w/2-Ueberhang kuerzt
   sich exakt weg. Die Breite muss also nicht bekannt sein, L und w folgen
   gemeinsam aus A und P:

       L und N*w = Wurzeln von  t^2 - (P/2)*t + N*A = 0
       (groessere Wurzel = L)

   Vor der Rechnung werden die selektierten Formen zu einer Region gemergt,
   damit aneinanderstossende Teilstuecke (Gerade + Bogen + S-Bogen ...)
   nicht falsche Stirnflaechen einrechnen; N ist dann die Anzahl der
   gemergten, zusammenhaengenden Baender.

   Restfehler nur aus der Polygon-Diskretisierung der Boegen in GDS
   (Sekantenzug statt Kreis -> minimal zu kurz). Bei PCell-Boegen mit
   feiner Rundung << 0.1 %. Wer es analytisch exakt braucht: die
   PCell-Parameter (radius/angle) werden mit ausgegeben, dann ist
   L_bogen = R * theta.
"""

import argparse
import math


def length_from_area_perimeter(area, perimeter, n_shapes=1):
    """Mittellinienlaenge L und Breite w aus Flaeche und Umfang.

    Loest t^2 - (P/2)*t + N*A = 0; groessere Wurzel = L, kleinere = N*w.
    """
    half_p = perimeter / 2.0
    disc = half_p * half_p - 4.0 * n_shapes * area
    if disc < 0:
        raise ValueError(
            "Diskriminante < 0 (%.6g): die Auswahl ist kein Band konstanter "
            "Breite (Coupler/Pad mitselektiert?)." % disc
        )
    root = math.sqrt(disc)
    return (half_p + root) / 2.0, (half_p - root) / 2.0 / n_shapes


def _run_in_klayout():
    import pya

    lv = pya.LayoutView.current()
    if lv is None:
        raise RuntimeError("Kein Layout-Fenster aktiv.")
    ly = lv.active_cellview().layout()
    dbu = ly.dbu

    region = pya.Region()
    n_sel = 0
    path_len = 0.0
    path_widths = []
    pcells = []

    for obj in lv.each_object_selected():
        if obj.is_cell_inst():
            # PCell-Parameter zur analytischen Kontrolle mitloggen.
            try:
                cell = ly.cell(obj.inst().cell_index)
                if cell.is_pcell_variant():
                    pcells.append((cell.name,
                                   cell.pcell_parameters_by_name()))
            except Exception:
                pass
            continue

        sh = obj.shape
        n_sel += 1

        # Transformation der Instanzierungskette in die Top-Zelle.
        try:
            trans = obj.trans
        except Exception:
            trans = pya.ICplxTrans()

        if sh.is_path():
            # Echtes Path-Objekt: Mittellinie direkt aus den Stuetzpunkten.
            p = sh.path.transformed(trans)
            pts = list(p.each_point())
            path_len += sum(
                (pts[i + 1] - pts[i]).length() for i in range(len(pts) - 1)
            ) * dbu
            path_widths.append(p.width * dbu)
            continue

        poly = None
        for attr in ("polygon", "simple_polygon"):
            try:
                poly = getattr(sh, attr)
            except Exception:
                poly = None
            if poly is not None:
                break
        if poly is None:
            n_sel -= 1
            print("  uebersprungen (kein Polygon/Pfad): %s" % sh)
            continue
        region.insert(pya.Polygon(poly).transformed(trans))

    if n_sel == 0 and not pcells:
        print("Nichts selektiert -- bitte das schwarze Routing auswaehlen.")
        return

    print("Selektierte Formen: %d" % n_sel)

    if path_len > 0:
        print("--- Path-Objekt(e): Mittellinie direkt summiert ---")
        print("Laenge  = %.6f um" % path_len)
        for w in path_widths:
            print("Breite  = %.6f um" % w)

    if not region.is_empty():
        region.merge()
        n_bands = region.count()
        area = region.area() * dbu * dbu
        perim = region.perimeter() * dbu
        print("--- Polygone: aus Flaeche + Umfang (nach Merge) ---")
        print("zusammenhaengende Baender N = %d" % n_bands)
        print("Flaeche = %.6f um^2" % area)
        print("Umfang  = %.6f um" % perim)
        try:
            L, w = length_from_area_perimeter(area, perim, n_bands)
            print("LAENGE  = %.6f um   <-- Mittellinie" % L)
            print("Breite  = %.6f um" % w)
            print("Kontrolle A/w = %.6f um" % (area / w))
        except ValueError as exc:
            print("Warnung: %s" % exc)
            print("Mit bekannter Breite w gilt: L = A/w")

    if pcells:
        print("--- PCell-Parameter der selektierten Instanzen ---")
        for name, params in pcells:
            print("  %s" % name)
            for k, v in sorted(params.items()):
                print("      %-16s = %s" % (k, v))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--area", type=float, required=True, help="Flaeche um^2")
    ap.add_argument("--perimeter", type=float, required=True,
                    help="Umfang um")
    ap.add_argument("--shapes", type=int, default=1,
                    help="Anzahl getrennter Baender (default 1)")
    a = ap.parse_args()
    L, w = length_from_area_perimeter(a.area, a.perimeter, a.shapes)
    print("Mittellinienlaenge = %.6f um" % L)
    print("Breite             = %.6f um" % w)


if __name__ == "__main__":
    try:
        import pya  # noqa: F401
        _run_in_klayout()
    except ImportError:
        main()
