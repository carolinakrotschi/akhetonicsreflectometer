"""Fiber-Box: 6 gestapelte Thorlabs-BFCT-Trays, dahinter 1 Turm a 3 Connector-Halter (je 3 Plaetze),
Seitenfach links fuer 4 Fiber-Rollen, Schiebedeckel, 6 Wand-Kupplungen zum Einstecken.

Koordinaten der Box: X = Breite (links->rechts), Y = Tiefe (vorne->hinten),
Z = Hoehe. Ursprung = aeussere Ecke vorne-links-unten. Alle Masse in mm.

Ausgabe (alles in step/):
  fiberbox.step             Box offen, Deckel daneben, Halter + Riegel montiert
  fiberbox_closed.step      dasselbe mit eingeschobenem Deckel
  fiberbox_with_trays.step  wie fiberbox.step plus 6 BFCT-Trays (nur Referenz, nicht drucken)
  box_print.step            Druck: Box (Schrift/Logo vertieft, einfarbig; danach ausmalen)
  lid_print.step, connector_holder_print.step   Druckteile einzeln, in Druckorientierung
  test_print.step           Schnell-Probedruck: Wand-Kupplung + Turm-Halter (je 1 Platz)

Ausfuehren:  .venv\\Scripts\\python.exe fiberbox.py
"""
from pathlib import Path

import cadquery as cq

HERE = Path(__file__).parent
OUT = HERE / "step"
TRAY_STEP = HERE.parent / "general documents" / "fiberholder.step"

# ---------------------------------------------------------------- Parameter
WALL = 2.0            # Wand (5 Perimeter bei 0.4er Duese)
FLOOR = 2.0
CLR = 1.0             # Spiel Tray <-> Wand je Seite

# Thorlabs BFCT (aus fiberholder.step gemessen; dort Y = hoch, Z = Tiefe)
TRAY_W = 147.92       # x
TRAY_D = 88.38        # z: -51.76 (Clip-Lasche) .. +36.62
TRAY_ZMAX = 36.62
TRAY_H = 12.83        # Einzelhoehe inkl. Schnapphaken
TRAY_PITCH = 11.2     # Stapelmass laut Thorlabs-Zeichnung (0.44 in), Haken rasten ein
N_TRAYS = 6
TRAY_HOLES = [(-38.1, -18.29), (38.1, 18.29)]   # Bohrung Ø7.74, diagonal genutzt
PIN_D = 7.2
PIN_H = (N_TRAYS - 1) * TRAY_PITCH + TRAY_H + 3.0

# Fiber-Rollen (Annahme Ø70 x 10, hochkant -- NACHMESSEN)
COIL_SLOT = 13.0
N_COILS = 4
COIL_DIV_T = 1.6
COIL_DIV_H = 45.0
COIL_D = 70.0
DIVIDER_T = 2.0       # Trennwand Seitenfach <-> Hauptraum
DIVIDER_H = 55.0      # darueber Durchgang fuer die Fibers
PASS_Z0 = 15.0        # waagrechtes Kabel-Fenster in allen Waenden des Rollenfachs,
PASS_Z1 = 30.0        # Unter-/Oberkante ueber dem Boden
PASS_END = 10.0       # Abstand zur Front-/Rueckwand
PASS_SEG = 35.0       # max. Fensterlaenge zwischen zwei Stegen (Bruecke beim Druck)
PASS_POST = 3.0       # Stegbreite

# FC-Kupplung: Flansch quadratisch 15 x 15 (+-0.5, Zeichnung), 4 mm dick (gemessen),
# Haelse Ø8.1 x 5.5 mm je Seite (gemessen). 13.2 war im Probedruck zu schmal; 17 auf Wunsch.
FC_FLANGE = 15.0
FC_POCKET = 17.0      # Schlitzbreite (Carolina)
FLANGE_T = 4.0

# Wand-Kupplung: Flansch sitzt IN der Wand (Schlitz, von oben eingesteckt), aussen nur
# eine duenne Lippe -> der aeussere Hals schaut ganz heraus. Lippe und Innenwand haben
# eine Oeffnung NECK_W, damit Mutter/Kappe bis an den Flansch kommen.
NECK_L = 5.5          # Halslaenge je Seite ab Flansch (gemessen)
MOUNT_LIP = 1.2       # aeussere Lippe vor dem Flansch
MOUNT_LEDGE = 2.0     # Material unter dem Flansch
MOUNT_HALF = FC_POCKET / 2 + 2.0   # Verdickung halbe Breite

# Connector-Halter (Modul, stapelbar auf 2 Stiften), 1 Turm, 3 Kupplungen pro Ebene
N_HOLDERS = 3         # Ebenen pro Turm
N_HOLDER_STACKS = 1
HOLDER_PITCH = 20.0
HOLDER_PLATE_T = 2.0
SLOT_T = FLANGE_T + 0.4
SLOT_WALL = 1.2       # duenn: Ueberwurfmutter muss weit genug auf den 5.5-mm-Hals
MOUNT_T = MOUNT_LIP + SLOT_T + SLOT_WALL       # Wanddicke an der Kupplung (ab Aussenseite)
NECK_W = 10.5         # U-Kerbe: Platz fuer gruene Schutzkappe / Mutter (Hals Ø8.1)
HOLDER_W = 24.0       # X
HOLDER_D = 78.0       # Y
HOLDER_SLEEVE_YS = (-20.0, 0.0, 20.0)
HOLDER_PINS = ((0.0, -34.0), (0.0, 34.0))      # an den Enden, nicht im Weg der Stecker
HOLDER_PIN_D = 5.0
HOLDER_PIN_CLR = 0.6
HOLDER_BOSS_D = 9.0
HOLDER_ZONE_D = 80.0
HOLDER_BACK_GAP = 14.5  # Luft hinter dem hinteren Turm: Wand-Kupplungen haengen nicht ueber den Halter-Steckern

# Schiebedeckel
LID_T = 2.0
LID_SLOT_H = 2.4      # 0.4 Spiel vertikal
LID_SIDE_CLR = 0.3
TOP_LIP = 1.2         # Material ueber der Nut
BAND_H = 8.0          # verdickter Rand oben an den Seitenwaenden
BAND_IN = 2.6         # so weit ragt der Rand nach innen (= Nuttiefe)

# Beschriftung / Logo als Inlay (eigener Koerper -> eigenes Filament im AMS)
INLAY_DEPTH = 0.8     # vertiefte Schrift, einfarbig gedruckt, zum Ausmalen
FONT = "Arial"
PORT_TEXT_H = 5.0
TITLE = "LINA OBR"
TITLE_H = 20.0
LOGO_SVG = HERE / "logo_akhetonics.svg"
LOGO_D = 30.0
FRONT_LABELS = ["LASER", "DUT"]
REAR_LABELS = ["CH1", "CH2", "CH3", "CH4"]
COLOR_BOX = cq.Color(0.80, 0.80, 0.82, 1)      # hellgrau

# ---------------------------------------------------------------- Ableitungen
coil_w = N_COILS * COIL_SLOT + (N_COILS - 1) * COIL_DIV_T
tray_x0 = WALL + coil_w + DIVIDER_T
tray_zone_w = TRAY_W + 2 * CLR
in_d = TRAY_D + 2 * CLR + HOLDER_ZONE_D + HOLDER_BACK_GAP

OUT_W = tray_x0 + tray_zone_w + WALL
OUT_D = WALL + in_d + WALL

tray_cx = tray_x0 + tray_zone_w / 2
tray_y_front = WALL + CLR
_zone_y0 = WALL + TRAY_D + 2 * CLR
HOLDER_CYS = [_zone_y0 + HOLDER_ZONE_D * (2 * k + 1) / (2 * N_HOLDER_STACKS)
              for k in range(N_HOLDER_STACKS)]

# Wand-Kupplungen oberhalb der Stapel; Hoehe mindestens fuer die stehenden Rollen
# Wand-Kupplungen so hoch wie moeglich (kurzer Einsteck-Schlitz von oben), aber
# hinten ueber dem Halter-Turm, vorne mit der 45-Grad-Verdickung ueber den Rollen-Stegen
_cz_min = max(FLOOR + N_HOLDERS * HOLDER_PITCH + 1.0 + MOUNT_LEDGE + FC_POCKET / 2,
              FLOOR + COIL_DIV_H + 1.0 + (MOUNT_T - WALL) + MOUNT_LEDGE + FC_POCKET / 2)
slot_z0 = max(_cz_min + FC_POCKET / 2 + 1.0,           # Nutunterkante
              FLOOR + COIL_D + 3.0)
panel_cz = slot_z0 - 1.0 - FC_POCKET / 2
OUT_H = slot_z0 + LID_SLOT_H + TOP_LIP

# vorne vor dem Rollenfach: ueber den Trays ist kein Platz fuer die 45-Grad-Schraege
FRONT_X = [WALL + MOUNT_HALF + 0.5, WALL + 3 * MOUNT_HALF + 2.5]   # Laser, DUT
REAR_X = [tray_cx - 51, tray_cx - 17, tray_cx + 17, tray_cx + 51]


def tray_to_box(x, z):
    """BFCT-Koordinaten (x, z) -> Box (X, Y). Clip-Lasche (z<0) zeigt nach hinten."""
    return tray_cx + x, tray_y_front + (TRAY_ZMAX - z)


# ---------------------------------------------------------------- Hilfen
def box(x0, y0, z0, dx, dy, dz):
    return cq.Workplane("XY").box(dx, dy, dz, centered=False).translate((x0, y0, z0))


def cyl_y(x, y0, z, d, length):
    return cq.Workplane("XY").add(
        cq.Solid.makeCylinder(d / 2, length, cq.Vector(x, y0, z), cq.Vector(0, 1, 0)))


def _wall_slab(face_y, d, u0, u1, w0, w1, z0, z1, x):
    """Quader in Wand-Koordinaten: u entlang der Wand (um x), w = Abstand von der
    Wand-Innenseite face_y in Richtung d (+1 Front, -1 Rueckwand)."""
    ya, yb = sorted((face_y + d * w0, face_y + d * w1))
    return box(x + u0, ya, z0, u1 - u0, yb - ya, z1 - z0)


def mount_geom(x, face_out, d, cz, z_top):
    """Wand-Kupplung: Verdickung + Flansch-Schlitz in der Wand. face_out = Aussenflaeche,
    d = Richtung nach innen (+1 Front, -1 Rueckwand). Schlitz und Hals-Oeffnung sind
    nach oben offen bis z_top. Gibt (add, cut) zurueck."""
    zf0 = cz - FC_POCKET / 2
    zb = zf0 - MOUNT_LEDGE
    extra = MOUNT_T - WALL                       # so weit ragt die Verdickung nach innen
    prof = [(face_out + d * 0, zb - extra), (face_out + d * MOUNT_T, zb),
            (face_out + d * MOUNT_T, z_top), (face_out + d * 0, z_top)]
    if d < 0:
        prof = prof[::-1]
    add = (cq.Workplane("YZ").polyline(prof).close().extrude(2 * MOUNT_HALF)
           .translate((x - MOUNT_HALF, 0, 0)))
    big = 50
    cut = _wall_slab(face_out, d, -FC_POCKET / 2, FC_POCKET / 2, MOUNT_LIP, MOUNT_LIP + SLOT_T,
                     zf0, z_top + big, x)
    cut = cut.union(_wall_slab(face_out, d, -NECK_W / 2, NECK_W / 2, -1, MOUNT_T + 1,
                               cz, z_top + big, x))
    y0 = min(face_out - d * 1, face_out + d * (MOUNT_T + 1))
    cut = cut.union(cyl_y(x, y0, cz, NECK_W, MOUNT_T + 2))
    return add, cut


def wall_mount(x, front):
    face, d = (0.0, +1) if front else (OUT_D, -1)
    return mount_geom(x, face, d, panel_cz, slot_z0)


def make_test_print():
    """Schnell-Probedruck mit minimalem Material: ein Turm-Halter-Platz und ein Stueck
    Wand mit einer Wand-Kupplung, beide in Originalgeometrie."""
    t = WALL
    w = 2 * MOUNT_HALF + 4.0
    cz = FLOOR + MOUNT_LEDGE + (MOUNT_T - WALL) + 1.0 + FC_POCKET / 2
    h = cz + FC_POCKET / 2 + 1.0
    wall = box(0, 0, 0, w, MOUNT_T + 4.0, FLOOR).union(box(0, 0, 0, w, t, h))
    add, cut = mount_geom(w / 2, 0.0, +1, cz, h)
    wall = wall.union(add).cut(cut)
    block_t = SLOT_T + 2 * SLOT_WALL
    tower = make_holder(sleeve_ys=(0.0,), pins=(), w=block_t + 4.0, d=FC_POCKET + 4.0)
    tower = tower.translate((w + 10 + (block_t + 4.0) / 2, (FC_POCKET + 4.0) / 2, 0))
    return [("test_wall_connector", wall), ("test_tower_holder", tower)]


# ---------------------------------------------------------------- Box
def make_box():
    b = box(0, 0, 0, OUT_W, OUT_D, OUT_H)
    b = b.cut(box(WALL, WALL, FLOOR, OUT_W - 2 * WALL, OUT_D - 2 * WALL, OUT_H))

    # Rand oben links/rechts mit 45-Grad-Fase darunter (stuetzfrei druckbar)
    for side in (0, 1):
        x_in = WALL if side == 0 else OUT_W - WALL
        s = 1 if side == 0 else -1
        prof = (cq.Workplane("XZ")
                .polyline([(x_in, OUT_H - BAND_H - BAND_IN),
                           (x_in + s * BAND_IN, OUT_H - BAND_H),
                           (x_in + s * BAND_IN, OUT_H),
                           (x_in, OUT_H)]).close()
                .extrude(-(OUT_D - WALL)))       # XZ extrudiert in -Y -> Vorzeichen
        b = b.union(prof)
        # Nut fuer den Deckel
        nx0 = x_in if side == 0 else x_in - BAND_IN
        b = b.cut(box(nx0 - (0.01 if side == 0 else 0), -1, slot_z0,
                      BAND_IN + 0.01, OUT_D - WALL + 1, LID_SLOT_H))

    # Frontwand oben bis Nutunterkante absenken -> Deckel schiebt von vorne ein
    b = b.cut(box(WALL, -1, slot_z0, OUT_W - 2 * WALL, WALL + 1.01, OUT_H))
    # Seitenwand-Material vor der Nut bleibt; Deckel-Frontblock fuellt den Ausschnitt.

    # Seitenfach: Trennwand + Rollen-Stege
    b = b.union(box(WALL + coil_w, WALL, FLOOR, DIVIDER_T, in_d, DIVIDER_H))
    for i in range(1, N_COILS):
        x = WALL + i * COIL_SLOT + (i - 1) * COIL_DIV_T
        b = b.union(box(x, WALL, FLOOR, COIL_DIV_T, in_d, COIL_DIV_H))

    # Waagrechtes Kabel-Fenster in allen Stegen + Trennwand, mit kurzen Stuetzstegen
    walls = [(WALL + i * COIL_SLOT + (i - 1) * COIL_DIV_T, COIL_DIV_T) for i in range(1, N_COILS)]
    walls.append((WALL + coil_w, DIVIDER_T))
    y0, y1 = WALL + PASS_END, WALL + in_d - PASS_END
    n = int((y1 - y0 + PASS_POST) // (PASS_SEG + PASS_POST)) + 1
    seg = (y1 - y0 - (n - 1) * PASS_POST) / n
    for wx, wt in walls:
        for j in range(n):
            ys = y0 + j * (seg + PASS_POST)
            b = b.cut(box(wx - 0.01, ys, FLOOR + PASS_Z0, wt + 0.02, seg, PASS_Z1 - PASS_Z0))

    # Tray-Stifte mit Fase oben
    for hx, hz in TRAY_HOLES:
        X, Y = tray_to_box(hx, hz)
        pin = (cq.Workplane("XY").circle(PIN_D / 2).extrude(PIN_H)
               .faces(">Z").edges().chamfer(1.0).translate((X, Y, FLOOR)))
        b = b.union(pin)

    # Halter-Stifte, 2 pro Turm
    for cy in HOLDER_CYS:
        for sx, sy in HOLDER_PINS:
            pin = (cq.Workplane("XY").circle(HOLDER_PIN_D / 2)
                   .extrude(N_HOLDERS * HOLDER_PITCH - 3)
                   .faces(">Z").edges().chamfer(0.8)
                   .translate((tray_cx + sx, cy + sy, FLOOR)))
            b = b.union(pin)

    # Wand-Kupplungen: 2 vorne, 4 hinten
    mounts = [wall_mount(x, True) for x in FRONT_X] + [wall_mount(x, False) for x in REAR_X]
    for add, _ in mounts:
        b = b.union(add)
    for _, cut in mounts:
        b = b.cut(cut)

    for _, inl in make_inlays():
        b = b.cut(inl)
    return b


# ---------------------------------------------------------------- Inlays
def _text_front(txt, x, z, size, y_face=0.0, outward=-1):
    """Schrift in der Front- (outward=-1) bzw. Rueckwand (+1), lesbar von aussen."""
    plane = cq.Plane(origin=(x, y_face - outward * INLAY_DEPTH, z),
                     xDir=(1, 0, 0) if outward < 0 else (-1, 0, 0),
                     normal=(0, outward, 0))
    return cq.Workplane(plane).text(txt, size, INLAY_DEPTH, font=FONT, kind="bold",
                                    halign="center", valign="center", combine=False)


def _logo(cx, cz, d):
    """Akhetonics-Logo (SVG, even-odd) als Inlay in der Frontwand, Mitte (cx, cz)."""
    import svgelements as se
    loops = []
    for el in se.SVG.parse(str(LOGO_SVG)).elements():
        if isinstance(el, se.Path):
            for sp in el.as_subpaths():
                sp = se.Path(sp)
                k = max(24, int(sp.length() / 4))
                loops.append([(sp.point(i / k).x, -sp.point(i / k).y) for i in range(k)])
    xs = [p[0] for l in loops for p in l]; ys = [p[1] for l in loops for p in l]
    sc = d / (max(xs) - min(xs))
    mx, my = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2
    loops = [[((u - mx) * sc, (v - my) * sc) for u, v in l] for l in loops]
    wires = [cq.Wire.makePolygon([cq.Vector(u, v, 0) for u, v in l], close=True) for l in loops]
    faces = [cq.Face.makeFromWires(w) for w in wires]

    def pip(pt, poly):  # Punkt-in-Polygon (Ray casting)
        x, y = pt; c = False
        for (x1, y1), (x2, y2) in zip(poly, poly[1:] + poly[:1]):
            if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1):
                c = not c
        return c

    def inside(i, j):  # liegt Loop i in Loop j?
        return i != j and faces[i].Area() < faces[j].Area() and pip(loops[i][0], loops[j])
    depth = [sum(inside(i, j) for j in range(len(loops))) for i in range(len(loops))]
    solid = None
    for i in range(len(loops)):
        if depth[i] % 2:
            continue
        holes = [wires[j] for j in range(len(loops))
                 if depth[j] == depth[i] + 1 and inside(j, i)]
        f = cq.Face.makeFromWires(wires[i], holes)
        sl = cq.Solid.extrudeLinear(f, cq.Vector(0, 0, -INLAY_DEPTH))
        solid = sl if solid is None else solid.fuse(sl)
    # XY (v nach oben) -> Frontwand: v -> Z, Extrusion -z -> +Y (in die Wand)
    solid = solid.rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), 90).translate(cq.Vector(cx, 0, cz))
    return cq.Workplane("XY").add(solid)


def make_inlays():
    """Liste (name, Workplane). Die Box bekommt dieselben Koerper als Vertiefung."""
    z_lab = panel_cz - FC_POCKET / 2 - 2.0 - PORT_TEXT_H / 2   # unter den Kupplungen
    out = [(f"label_front_{t}", _text_front(t, x, z_lab, PORT_TEXT_H))
           for t, x in zip(FRONT_LABELS, FRONT_X)]
    out += [(f"label_rear_{t}", _text_front(t, x, z_lab, PORT_TEXT_H, OUT_D, +1))
            for t, x in zip(REAR_LABELS, reversed(REAR_X))]   # von hinten gelesen: CH1 links
    out.append(("label_title", _text_front(TITLE, tray_cx, 38.0, TITLE_H)))
    out.append(("logo", _logo(WALL + coil_w / 2, 28.0, LOGO_D)))
    out.append(("label_akhetonics", _text_front("AKHETONICS", OUT_W / 2, 30.0, 16.0, OUT_D, +1)))
    return out


# ---------------------------------------------------------------- Deckel
def make_lid():
    lid_w = OUT_W - 2 * WALL - 2 * LID_SIDE_CLR            # liegt auf dem Rand, Lippe drueber
    lid_d = OUT_D - WALL - LID_SIDE_CLR                    # bis Rueckwand-Innenseite
    lid = box((OUT_W - lid_w) / 2, 0, slot_z0 + (LID_SLOT_H - LID_T) / 2,
              lid_w, lid_d, LID_T)
    # Frontblock: fuellt den Ausschnitt der Frontwand ueber dem Deckel und darunter
    fb_w = OUT_W - 2 * WALL - 2 * LID_SIDE_CLR
    lid = lid.union(box(WALL + LID_SIDE_CLR, 0, slot_z0 + 0.2, fb_w, WALL,
                        OUT_H - slot_z0 - 0.2))
    # Griffrillen vorne oben
    top = slot_z0 + (LID_SLOT_H - LID_T) / 2 + LID_T
    for i in range(4):
        lid = lid.cut(box(OUT_W / 2 - 30, 6 + i * 3, top - 0.6, 60, 1.2, 1))
    return lid


# ---------------------------------------------------------------- Halter
def make_holder(sleeve_ys=HOLDER_SLEEVE_YS, pins=HOLDER_PINS, w=HOLDER_W, d=HOLDER_D):
    """Steck-Halter fuer 3 FC-Kupplungen (wie der gedruckte Halter im Labor):
    Flansch wird von oben in einen Schlitz gesteckt, der Hals liegt in einer U-Kerbe.
    Der naechste Halter im Stapel liegt oben auf und haelt die Flansche fest.
    Lokal: Mitte bei X=Y=0, Unterseite Z=0, Kupplungsachse = X."""
    block_t = SLOT_T + 2 * SLOT_WALL
    block_h = HOLDER_PITCH - HOLDER_PLATE_T
    h = cq.Workplane("XY").box(w, d, HOLDER_PLATE_T, centered=(True, True, False))
    h = h.union(cq.Workplane("XY").box(block_t, d, block_h, centered=(True, True, False))
                .translate((0, 0, HOLDER_PLATE_T)))
    for sx, sy in pins:
        h = h.union(cq.Workplane("XY").circle(HOLDER_BOSS_D / 2).extrude(HOLDER_PITCH)
                    .translate((sx, sy, 0)))
        h = h.cut(cq.Workplane("XY").circle((HOLDER_PIN_D + HOLDER_PIN_CLR) / 2)
                  .extrude(HOLDER_PITCH + 2).translate((sx, sy, -1)))
    flange_z0 = HOLDER_PLATE_T + 0.5
    cz = flange_z0 + FC_FLANGE / 2
    for sy in sleeve_ys:
        # Flansch-Schlitz, oben offen
        h = h.cut(cq.Workplane("XY").box(SLOT_T, FC_POCKET, HOLDER_PITCH, centered=(True, True, False))
                  .translate((0, sy, flange_z0)))
        # U-Kerbe fuer den Hals: Halbkreis unten, nach oben offen
        neck = (cq.Workplane("YZ").center(sy, cz).circle(NECK_W / 2).extrude(block_t + 2)
                .union(cq.Workplane("YZ").center(sy, cz + HOLDER_PITCH / 2)
                       .rect(NECK_W, HOLDER_PITCH).extrude(block_t + 2))
                .translate((-(block_t + 2) / 2, 0, 0)))
        h = h.cut(neck)
    return h


# ---------------------------------------------------------------- Trays (Referenz)
def load_tray():
    t = cq.importers.importStep(str(TRAY_STEP))
    # STEP: Y hoch, Z Tiefe  ->  Box: Z hoch, Y Tiefe (z -> -Y)
    return t.rotate((0, 0, 0), (1, 0, 0), 90)


def build_assembly(lid_pos, with_trays):
    """lid_pos: "closed" (eingeschoben), "beside" (neben der Box, Box offen) oder None."""
    bx, lid, holder = make_box(), make_lid(), make_holder()
    asm = cq.Assembly(name="fiberbox")
    asm.add(bx, name="box", color=COLOR_BOX)
    if lid_pos == "closed":
        asm.add(lid, name="lid", color=COLOR_BOX)
    elif lid_pos == "beside":
        asm.add(lid, name="lid", color=COLOR_BOX,
                loc=cq.Location((OUT_W + 30, 0, -slot_z0)))
    for k, cy in enumerate(HOLDER_CYS):
        for i in range(N_HOLDERS):
            asm.add(holder, name=f"holder_{k+1}_{i+1}",
                    loc=cq.Location((tray_cx, cy, FLOOR + i * HOLDER_PITCH)),
                    color=COLOR_BOX)
    if with_trays:
        tray = load_tray()
        for i in range(N_TRAYS):
            asm.add(tray, name=f"bfct_{i+1}",
                    loc=cq.Location((tray_cx, tray_y_front + TRAY_ZMAX, FLOOR + i * TRAY_PITCH)),
                    color=cq.Color(0.95, 0.95, 0.95, 1))
    return asm


def main():
    # Haupt-Datei: Box offen, Deckel daneben -> man sieht von oben hinein
    OUT.mkdir(exist_ok=True)
    parts = OUT
    build_assembly("beside", False).save(str(OUT / "fiberbox.step"))
    build_assembly("beside", True).save(str(OUT / "fiberbox_with_trays.step"))
    build_assembly("closed", False).save(str(OUT / "fiberbox_closed.step"))

    bx, lid, holder = make_box(), make_lid(), make_holder()
    cq.exporters.export(bx, str(parts / "box_print.step"))
    cq.exporters.export(lid.translate((0, 0, -slot_z0)), str(parts / "lid_print.step"))
    cq.exporters.export(holder, str(parts / "connector_holder_print.step"))
    test = cq.Assembly(name="test_print")
    for name, w in make_test_print():
        test.add(w, name=name, color=COLOR_BOX)
    test.save(str(parts / "test_print.step"))

    print(f"Box aussen: {OUT_W:.1f} x {OUT_D:.1f} x {OUT_H:.1f} mm")
    print(f"Seitenfach innen: {coil_w:.1f} x {in_d:.1f} mm, Trays-Raum {tray_zone_w:.1f} x {TRAY_D+2*CLR:.1f}")
    print(f"Panel-Connector-Mitte Z = {panel_cz:.1f}, Nut ab Z = {slot_z0:.1f}")


if __name__ == "__main__":
    main()
