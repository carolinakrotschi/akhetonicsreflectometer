# OFDR-Reflektometer — Erklärung zum Handover

**Für:** Carolina · **Stand:** 18.08.2026
**Bezieht sich auf:** `HANDOVER.md`, `ofdr_process.py`, `ofdr_diagnose.py` und die
E-Mail deines Supervisors zum neuen MZI

Dieses Dokument erklärt (A) das Handover — jede Formel, jeden Begriff, von Grund
auf — und (B) den Vorschlag deines Supervisors, übersetzt und durchgerechnet,
mit den Zahlen, die du morgen am Tisch brauchst.

---

# TEIL A — Das Handover verstehen

## A.1 Was das Gerät überhaupt tut

Ein **OFDR** (Optical Frequency-Domain Reflectometer) beantwortet die Frage:
*"Wo entlang dieser Faser sitzt welche Reflexion, und wie stark ist sie?"*
Also im Prinzip eine Landkarte der Faser: hier ein Stecker, dort ein Spleiß,
am Ende die Faserendfläche.

Ein klassisches OTDR macht das mit einem kurzen Puls und einer Stoppuhr
(Zeitbereich). Auflösung ~cm bis m. Ein OFDR macht es anders und viel
feiner: es schickt **kein** Puls, sondern einen Laser, der seine Wellenlänge
kontinuierlich durchfährt (*sweep*), und misst **Interferenz**. Auflösung:
Mikrometer.

Der Aufbau in einem Satz: Licht wird geteilt. Ein Teil geht über einen
**Zirkulator** in die Testfaser und kommt reflektiert zurück. Der andere Teil
läuft als unangetastete Referenz (**Lokaloszillator, LO**) nebenher. Beide
werden wieder zusammengeführt und interferieren auf einem Detektor.

## A.2 Die eine Grundgleichung — hergeleitet

Im Handover steht:

```
I(ν) = Σᵢ 2·√(P_LO · Pᵢ) · cos(2π · τᵢ · ν + φᵢ)
```

Woher kommt das?

**Schritt 1 — Zwei Felder überlagern sich.** Auf dem Detektor treffen das
LO-Feld E_LO und ein reflektiertes Feld Eᵢ zusammen. Ein Detektor misst
Leistung, also Betragsquadrat:

```
P = |E_LO + Eᵢ|² = |E_LO|² + |Eᵢ|² + 2·Re(E_LO* · Eᵢ)
                   └──DC──┘  └klein┘  └── der interessante Teil ──┘
```

Der dritte Term ist der **Interferenzterm**. Er ist proportional zu
`2·√(P_LO·Pᵢ)` — das ist der Grund, warum ein OFDR so empfindlich ist: ein
extrem schwaches Rücksignal Pᵢ wird durch die Wurzel und durch den starken LO
*verstärkt*. Das nennt man **kohärente Detektion** oder heterodynen Gewinn.

**Schritt 2 — Woher der Kosinus kommt.** Das reflektierte Licht ist später
zurück als der LO, um die Laufzeit **τᵢ** (griech. tau, *round-trip delay*).
Eine Verzögerung um τ bedeutet bei optischer Frequenz ν eine Phasendifferenz
von `2π·τ·ν`. Daher `cos(2π·τᵢ·ν + φᵢ)`.

**Schritt 3 — Der Trick.** Jetzt fährt der Laser ν durch. Weil τᵢ **fest**
ist, aber ν sich ändert, läuft die Phase durch — der Detektor sieht eine
**Schwebung** (*fringes*, Interferenzstreifen). Und jetzt der Kern:

> Trägt man das Detektorsignal gegen die **optische Frequenz ν** auf, ist es
> eine Summe von Sinusschwingungen. Die "Frequenz" jeder Sinusschwingung
> (gemessen in *pro Hz optischer Frequenz*, also in Sekunden) **ist** genau
> die Laufzeit τᵢ des jeweiligen Reflektors.

Also: **FFT über ν → jeder Reflektor wird ein Peak, und seine Position auf der
FFT-Achse ist seine Laufzeit, also seine Entfernung.** Das ist das gesamte
Prinzip. Der Rest der Skripte ist Sorgfalt darum herum.

Ein Bild dazu: ein naher Reflektor gibt eine langsame Schwebung (wenige
Streifen über den ganzen Sweep), ein weit entfernter eine schnelle. Wie beim
Radar-FMCW oder bei OCT in der Medizin — dieselbe Mathematik.

## A.3 Die Tabelle der Formeln, Zeile für Zeile

### Auflösung: `Δz = c / (2 · n_g · Δν_span)` → **7,0 µm**

*Δν_span* = die gesamte durchfahrene optische Frequenzspanne. Von 1505 nm bis
1625 nm sind das **14,7 THz**.

**Warum diese Formel?** Fourier-Grundregel: wer ein Signal über ein Fenster der
Breite B beobachtet, kann Frequenzen bis auf 1/B unterscheiden. Hier ist das
Fenster die Frequenzspanne Δν_span, also unterscheidbare Laufzeiten
`Δτ = 1/Δν_span`. Umgerechnet in Entfernung: `Δz = c·Δτ/(2·n_g)`.

**Merksatz:** *Die Auflösung hängt NUR von der Breite des Sweeps ab.* Mehr nm
durchfahren = feiner auflösen. Nichts anderes hilft.

### Nyquist-Reichweite: `z_max = c / (4 · n_g · δν_step)` → **0,42 m**

*δν_step* = der Abstand **zwischen zwei Messpunkten** in optischer Frequenz
(bei 1 pm Schritt: ca. 122 MHz).

**Warum?** Nyquist-Abtasttheorem: eine Schwingung braucht mindestens 2 Punkte
pro Periode, sonst wird sie falsch rekonstruiert. Der Reflektor bei τ erzeugt
eine Schwingung mit Periode `1/τ` in ν. Bedingung: `1/τ > 2·δν_step`, also
`τ_max = 1/(2·δν_step)` und in Entfernung `z_max = c/(4·n_g·δν_step)`.

**Merksatz:** *Die Reichweite hängt NUR davon ab, wie eng deine Messpunkte
liegen.* Kleinerer Schritt = größere Reichweite.

### Die wichtigste Konsequenz aus beidem

```
Reichweite / Auflösung = N/2       (N = Anzahl Messpunkte)
```

Aktuell: 120.000 Punkte → 60.000 "Zellen" → 0,42 m ÷ 7 µm ≈ 60.000. ✔

**Das ist das ganze Budget.** Du kannst Reichweite gegen Auflösung tauschen,
aber das Produkt ist durch die Punktzahl fixiert. Willst du *beides* besser:
mehr Punkte. Punkt. Und **Sweep-Geschwindigkeit kürzt sich überall raus** —
sie entscheidet nur, ob dein Digitizer schnell genug ist bzw. ob dein Speicher
reicht.

### Punktabstand in ν: `δν = c·δλ/λ²` → **132 → 114 MHz (18 % Chirp!)**

Das ist nur die Umrechnung Wellenlänge ↔ Frequenz: `ν = c/λ`, differenziert
`|dν| = c·dλ/λ²`.

**Der Haken:** Der Laser triggert alle **1 pm in λ** — aber weil λ² im Nenner
steht, ist 1 pm bei 1505 nm ein Frequenzschritt von 132 MHz und bei 1625 nm
nur noch 114 MHz. Deine Punkte sind also gleichmäßig in λ, aber **ungleichmäßig
in ν** — 18 % Unterschied über den Sweep. Das nennt das Handover *Chirp*.

Eine FFT setzt gleichmäßige Abstände voraus. Wenn du die Rohdaten einfach
FFTst, verschmiert jeder Peak. Deshalb ist Schritt 3 der Pipeline
(λ→ν-Umrechnung mit Neuabtastung) **Pflicht, nicht Kosmetik**.

### Streifenperiode: `λ²/(2·n_g·z)` → **0,78 pm bei z = 1,05 m**

Umkehrung der gleichen Rechnung: wie viele Pikometer muss der Laser fahren,
damit ein Reflektor bei z eine volle Interferenzperiode durchläuft.

Bei z = 1,05 m ist die Periode **0,78 pm** — aber du misst nur alle 1 pm.
**Weniger als 2 Punkte pro Periode = unterabgetastet.** Genau daher kommt das
große Problem in §5 des Handovers (siehe A.6).

### Gruppenindex `n_g = 1,468`

Licht in Glasfaser ist langsamer als im Vakuum: `v = c/n_g`. Für Laufzeiten
zählt der **Gruppen**index (Geschwindigkeit des Signals/Wellenpakets), nicht
der Phasenindex (~1,4682 vs. ~1,468). Für SMF-28 bei 1550 nm: 1,468.

Praktisch: **1 m Faser ≈ 4,9 ns einfacher Weg, also 9,8 ns hin und zurück.**

## A.4 Die Distanz-Konvention — die Verwirrungsquelle

Alle Skripte rechnen `z = c·τ / (2·n_g)`. Der Faktor 2 heißt: *"τ ist ein
Hin- und Rückweg"* — Reflexionskonvention. z ist die einfache Strecke bis zum
Reflektor.

**Aber:** Ein **Transmissions-MZI** (Licht läuft nur einmal durch, keine
Reflexion) mit Armdifferenz ΔL erzeugt τ = n_g·ΔL/c — nur *ein* Weg. Das
Skript teilt trotzdem durch 2 und zeigt daher `z = ΔL/2`.

Der dominante Peak bei **46,4 mm** ist also entweder
- eine MZI-Armdifferenz von **92,8 mm** (Transmissions-Lesart), **oder**
- ein echter Reflektor **46,4 mm** hinter dem LO-Abgleichpunkt.

Das Handover sagt ausdrücklich: **das ist noch nicht geklärt, ein Maßband auf
dem Tisch klärt es.** Das ist deine Aufgabe (1). Nicht überspringen — die
gesamte Interpretation der Topologie hängt daran.

## A.5 Die Pipeline (`ofdr_process.py`) Schritt für Schritt

**1. Laden + Aufräumen.** Der letzte Messpunkt ist der *retrace*: der Laser
springt nach dem Sweep zurück unter die Startwellenlänge. Wird per **Wert**
verworfen (nicht per Position), weil zusätzlich ~75 Trigger am Anfang und ~33
am Ende fehlen (Anlaufphase des Lasers).

> **Regel aus dem Handover: nie zwei Scans über den Index vergleichen, immer
> über den Wellenlängenwert.** Sonst vergleichst du versetzte Datenpunkte.

**2. Balancierte Subtraktion in Software.** `P = Ch1 − g(λ)·Ch2`.

Der letzte Koppler hat **zwei komplementäre Ausgänge**: wo einer hell ist, ist
der andere dunkel (Energieerhaltung). Beide Ausgänge tragen dasselbe
Interferenzsignal, aber gegenphasig — und denselben DC-Sockel und dasselbe
Intensitätsrauschen **gleichphasig**. Also:

- Differenz Ch1 − Ch2 → Interferenz **verdoppelt**, DC & gemeinsames Rauschen **weg**
- Summe Ch1 + Ch2 → Interferenz weg, reine Leistungsüberwachung geschenkt

`g` ist ein Korrekturfaktor, weil die beiden Photodioden unterschiedlich
empfindlich sind und das über 120 nm driftet (0,84 → 1,14). Wird
abschnittsweise per Median gefittet.

*(Das ist die Software-Version. Echte balancierte Photodetektoren machen das
in Hardware und sind viel besser — siehe Handover §6, Punkt 1.)*

**3. λ → ν umrechnen und neu abtasten.** Wegen des 18 %-Chirps. Pflicht.

**4. "Slow-bow"-Korrektur.** Der Laser fährt nicht exakt so, wie er behauptet.
Die tatsächliche Abweichung vom nominellen 1-pm-Raster beträgt **~59 pm
Spitze-Spitze** in einem langsamen Bogen über den Sweep (siehe rechte Hälfte
von `real_data_diagnosis.png` — der Bogen von −25 pm über +18 pm zurück auf
−40 pm).

Ohne Korrektur zerfasert der Hauptpeak über ±0,8 mm.

Wie korrigiert wird: man nimmt den dominanten Ton, holt sich per **Hilbert-
Transformation** (im Code: `analytic_signal`) seine momentane Phase, glättet
sie stark, und verzerrt damit die ν-Achse.

> **Hilbert / analytisches Signal in einem Satz:** ein reeller Kosinus wird in
> einen rotierenden Zeiger umgewandelt, dessen Winkel man direkt als Phase
> ablesen kann — man löscht dazu einfach die negativen Frequenzen im Spektrum.
> Das steht genau so in den drei Zeilen von `analytic_signal()`.

Ergebnis: Hauptpeak **21,0 µm** breit bei einer theoretischen Grenze von
18,2 µm → **15 % vom Ideal, sehr gut.**

> **Die Peak-Breite ist der Gesundheitsindikator des Systems.** 46,4-mm-Peak
> ≤ ~21 µm mit Kaiser-12 → alles in Ordnung. Größer → erst untersuchen, dann
> weitermessen. Das Skript warnt bei > 2×.

**Der `--mode`-Schalter (wichtig, kein Kosmetikdetail):**

| Modus | Was er tut | Wann benutzen |
|---|---|---|
| `diagnostic` (Standard) | Korrektur folgt nur dem *langsamen* Bogen | **Immer, wenn du etwas identifizieren willst.** Echte Merkmale bleiben scharf und an der richtigen Stelle |
| `cosmetic` | Korrektur folgt der Rohphase | Nur zum Vorzeigen. **Frisst schwache Merkmale (−8 dB gemessen) und verschiebt sie.** Niemals daraus diagnostizieren |
| `none` | keine Korrektur | Demo, wie schlimm es ohne wäre |

Das siehst du direkt in der linken Hälfte von `real_data_diagnosis.png`: blau
= rohes EXFO-Raster (Peak versinkt in Gras bei −10 dB), orange =
phasenkorrigiert (sauberer Peak, Untergrund bei −45 dB).

**5. Fenstern + FFT.** Ein endlich langer Datensatz erzeugt bei der FFT
**Nebenkeulen** (*sidelobes*) — falsche kleine Peaks neben jedem echten. Ein
Fenster (weiche Gewichtung zu den Rändern hin) drückt die.

- **Hann**: schmaler Hauptpeak, aber Nebenkeulen nur −31 dB → in einem OFDR
  werden die für Reflektoren gehalten. Gefährlich.
- **Kaiser β = 12** (Standard hier): Hauptpeak 2,6 Bins breit, Nebenkeulen
  **< −80 dB**. Richtige Wahl, wenn man schwache echte Reflexionen sucht.

**Bekannte Harmonische:** Peaks bei exakt 2×, 3× der Hauptverzögerung
(92,8 mm, 139,2 mm) sind **keine Reflektoren**, sondern Verzerrungsprodukte
von Detektor/Verarbeitung. Die balancierte Subtraktion löscht geradzahlige und
lässt ungeradzahlige durch, die Summe umgekehrt. Der Peak-Lister markiert sie
automatisch.

## A.6 Die Artefakt-Geschichte (§5) — und warum sie die wichtigste Seite ist

**Symptom:** ein diffuses Band von Peaks zwischen 6 und 32 cm, −38 bis −53 dB.
Zwischen zwei Scans blieb die **Einhüllende** identisch (Korrelation +0,979),
aber die **Feinstruktur** war jedes Mal völlig neu gewürfelt (+0,014).

**Vier Erklärungen wurden getestet, drei starben:**

1. *Laser-FM-Rauschen* → getötet vom **PM/AM-Test**. Logik: eine echte
   zusätzliche Lichtstrecke moduliert Amplitude und Phase des Haupttons
   gleich stark; reines Phasenrauschen moduliert fast nur die Phase. Gemessen:
   **PM/AM = 1,04** → es ist echtes Licht auf einem echten Weg (**additiv**).
2. *Mehrfachreflexion an Steckern* → getötet, weil feste Wege ihre Laufzeit auf
   0,1 % halten; dieses Band wanderte um 26 %.
3. *Polarisationseffekt im Zirkulator* → getötet von Scan C: Aufwickeln
   (= Dämpfung) schwächte das Band, reines Verdrehen der Polarisation nicht.
   Verlust-empfindlich + Polarisations-unempfindlich = keine Polarisation.
4. *Fabry-Pérot in der Endkappe* → getötet von Scan B: Kappe ab, Band
   unverändert.

**Die Wahrheit: es ist das Faserende bei ~1,05 m — aliasiert.**

**Aliasing anschaulich:** Wie das Wagenrad im Film, das rückwärts zu drehen
scheint. Die Kamera filmt zu langsam für die echte Drehung, also erfindet das
Auge eine falsche, langsamere. Genauso hier: das Faserende bei 1,05 m
erzeugt eine Streifenperiode von **0,78 pm**, du misst aber nur alle **1 pm**
— unterabgetastet. Also erscheint es an einer **falschen** Position:

```
Nyquist-Reichweite 0,42 m; Faltungsperiode 2 × 0,42 = 0,83 m
1,05 m − 0,83 m = 0,22 m  →  erscheint bei ~21 cm ✔ (gemessen: 21 cm)
```

`python ofdr_diagnose.py fold 1.046` rechnet dir das vor.

Und weil es unterabgetastet ist, kann die Neuabtastung es nicht kohärent
behandeln — das Abstimmrauschen des Lasers zerreißt es in ein Rauschband,
das bei jedem Sweep **neu ausgewürfelt** wird.

**Die drei Konsequenzen, die du dir merken musst:**

1. **Aliasing ist in Software nicht reparierbar.** Die Information geht beim
   Abtasten verloren, nicht danach. Es gibt kein optisches Anti-Alias-Filter.
   Nur zwei Auswege: (a) keine Reflektoren jenseits Nyquist zulassen
   (terminieren, kürzen, Index-Matching), oder (b) **feiner abtasten** — und
   genau das will dein Supervisor jetzt bauen.
2. Solange das Faserende dort steht, liegt in 6–32 cm ein **Rauschboden bei
   ~−40 dB**. Alles Schwächere dort ist unmessbar.
3. **Nicht wegkalibrierbar** durch einen Referenzscan — die Feinstruktur ist
   jedes Mal anders.

**Und die Arbeitsmethode, die du übernehmen sollst:** drei plausible Geschichten
waren falsch. Jede wurde von einer **Messung** getötet, nicht von einem besseren
Argument. `ofdr_diagnose.py` ist genau diese vier Tests in Kommandoform:

```bash
python ofdr_diagnose.py pm-am  scan.json          # additiv oder Phasenrauschen?
python ofdr_diagnose.py compare scanA.json scanB.json  # fest oder zufällig?
python ofdr_diagnose.py fold   1.046              # wohin faltet was?
python ofdr_diagnose.py tuning scan.json          # wie gut fährt der Laser heute?
```

## A.7 Glossar

| Begriff | Bedeutung |
|---|---|
| **MZI** | Mach-Zehnder-Interferometer: aufteilen, zwei verschieden lange Arme, wieder zusammenführen. Erzeugt einen sauberen sinusförmigen Streifen als Funktion von ν |
| **LO** (Local Oscillator) | Der Referenzstrahl, der nicht durch die Testfaser läuft und mit dem Rücksignal interferiert |
| **Zirkulator** | 3-Tor-Bauteil, Licht geht 1→2→3 aber nie rückwärts. Trennt hinlaufendes von rücklaufendem Licht |
| **eVOA** | Elektronisch verstellbarer Abschwächer. Hier zum Anpassen der LO-Leistung |
| **Balanced detection** | Beide komplementären Koppler-Ausgänge messen und subtrahieren |
| **Visibility** (0,70–0,88) | Streifenkontrast, `(max−min)/(max+min)`. 1 = perfekt gleich starke Felder. 0,745 heißt: die beiden Felder liegen ~7 dB auseinander |
| **SOP** | State of Polarization. Ändert sich in normaler SM-Faser bei jeder Berührung |
| **SM / PM Faser** | Single-Mode (Polarisation läuft frei) / Polarization-Maintaining (Polarisation bleibt fest, teurer) |
| **MPI** | Multi-Path Interference: unerwünschtes Licht über einen zweiten Weg |
| **Rayleigh-Rückstreuung** | Das extrem schwache, überall in der Faser gestreute Licht (~−100 dB/mm). Damit misst ein Luna OBR Dehnung/Temperatur — braucht echte balancierte Detektoren |
| **Fringe** | Ein Interferenzstreifen = eine volle Periode des Kosinus |
| **Chirp** | Hier: dass 1 pm nicht überall dieselbe Frequenzänderung bedeutet |
| **Slow bow** | Langsamer Bogen der Abstimmfehler des Lasers über den Sweep (~59 pm p-p) |
| **Retrace** | Rücklauf des Lasers nach dem Sweep |
| **Transform-limited** | So schmal wie die Physik es maximal erlaubt |
| **Aux-Interferometer** | Ein zweites, mitlaufendes MZI, das ausschließlich als **Frequenz-Lineal** dient. **← das, was dein Supervisor bauen lassen will** |
| **Dispersion / autofocus** | Verschiedene Wellenlängen laufen unterschiedlich schnell → weit entfernte Peaks verschmieren; korrigierbar mit quadratischem Phasenterm |
| **DGD** | Differential Group Delay: in PM-Faser laufen die zwei Polarisationsachsen unterschiedlich schnell (~1,5 ps/m). Deshalb nie SM→PM im Rückweg spleißen |

---

# TEIL B — Der Vorschlag deines Supervisors, übersetzt

## B.1 Was er sagt, in einem Absatz

Das Problem: hinter dem Zirkulator hängen **~2 m Testfaser** (1 m Patchcord +
1 m aus dem Bauteil). Deine Nyquist-Reichweite ist aber **0,42 m**. Also
aliasiert alles jenseits von 0,42 m in dein Messfenster hinein und macht genau
den Dreck, den §5 des Handovers beschreibt.

Die Lösung: **viel feiner abtasten als alle 1 pm.** Der Laser kann aber nur
alle 1 pm triggern. Also: **den Trigger komplett weglassen** und stattdessen
mit dem *internen Takt* des Leistungsmessers abtasten (alle 1 µs). Dann hast du
zwar viel mehr Punkte — aber **keine Ahnung mehr, bei welcher Wellenlänge
welcher Punkt aufgenommen wurde.** Und die brauchst du zwingend, weil die FFT
über ν läuft.

Deshalb das neue MZI: **es ist dein Lineal.** Es misst bei jedem Zeitpunkt mit,
wie weit der Laser gerade gefahren ist. Aus seiner Streifenphase rekonstruierst
du die ν-Achse.

Das ist exakt die Architektur, die in Handover §6 schon als Entscheidung
festgehalten ist ("free-running internal-clock sampling — the Luna
architecture"). Dein Supervisor sagt dasselbe, nur mit konkreten Bauteilen.

## B.2 ⚠ Terminologie-Falle: es gibt zwei "Auflösungen"

Er schreibt *"we'd need a much better resolution than 1pm"*. Das Handover
schreibt *"Resolution 7,0 µm"*. Das sind **zwei verschiedene Dinge** und du
verwechselst sie sonst garantiert:

| | Was gemeint ist | Formel | Wovon abhängig |
|---|---|---|---|
| **Abtastschritt** (er meint das) | Abstand zwischen zwei Messpunkten, 1 pm | δν = c·δλ/λ² | Trigger bzw. Taktrate |
| **Ortsauflösung** (Handover) | Wie nah zwei Reflektoren sein dürfen, 7 µm | Δz = c/(2n_g·Δν_span) | **nur** Sweep-Breite |

Der Abtastschritt bestimmt die **Reichweite**. Die Sweep-Breite bestimmt die
**Auflösung**. Er will die Reichweite verbessern, nicht die Ortsauflösung.

Wenn du das im Gespräch sauber trennst, verstehst du sofort, warum die
Rechnungen unten so aussehen, wie sie aussehen.

## B.3 Der Aufbau, den du bauen sollst

```
        EXFO Laser
             │
        ┌────┴────┐  90/10 Splitter
        │90%      │10%
        │         │
        │      ┌──┴──┐  50/50 Splitter
        │      │     │
        │   kurzer  langer Arm  (Armdifferenz ΔL — siehe B.5!)
        │      │     │
        │      └──┬──┘  50/50 Koppler (rekombinieren)
        │         ├────────► Ch3   ┐ AUX-MZI = dein Lineal
        │         └────────► Ch4   ┘ (komplementäre Ausgänge)
        │
   20 dB Koppler ──────► LO-Arm ──┐
        │                          │
    Zirkulator ──► TESTFASER       │
        │  ◄── Reflexionen ────────┤
        │                          │
        └──────────────► Koppler ──┘
                          ├────────► Ch1   ┐ MESS-Interferometer
                          └────────► Ch2   ┘ (wie bisher)
```

**Warum 4 Kanäle:** Ch1/Ch2 = die zwei komplementären Ausgänge deiner
bestehenden Messung. Ch3/Ch4 = die zwei komplementären Ausgänge des neuen
Aux-MZI. Bei beiden Paaren machst du die balancierte Subtraktion (siehe A.5,
Schritt 2) — beim Aux ist das besonders wichtig, weil du seine Phase sehr
sauber brauchst.

**Warum 90/10 und nicht 50/50:** Das Aux-MZI hat kaum Verlust (kein
Zirkulator, keine Testfaser), es braucht also viel weniger Licht. Die 90 %
gehen dahin, wo die Verluste sind.

**Wichtig aus Handover §6 (nicht übersehen):** Wenn du mit dem Aux
neu abtastest, benutzt du dessen **rohe, ungeglättete** Phase — nicht wie in
`--mode diagnostic` eine geglättete. Der ganze Diagnostic/Cosmetic-Konflikt
verschwindet dann, weil die Korrektur nicht mehr aus dem Messsignal selbst
kommt. **Aux-Phase niemals glätten** — sie muss die ~4,7-pm-Welligkeit
mittragen.

## B.4 Die Rechnung: Buffer vs. Sweep vs. Span

Der schwarze CoreDAQ füllt seinen Puffer mit **1 Punkt alle 1 µs**. Daraus
folgt zwingend:

```
Aufnahmezeit T = N_Punkte × 1 µs
Sweep-Span    = Sweep-Geschwindigkeit × T
δλ pro Punkt  = Sweep-Geschwindigkeit × 1 µs
```

Und daraus wieder die zwei Größen, die du willst. Hier die Tabelle — **die
solltest du dir ausdrucken**:

| Punkte/Kanal | Sweep | T | Span | δλ | **Reichweite z_max** | **Auflösung Δz** |
|---|---|---|---|---|---|---|
| 1.000.000 | 120 nm/s | 1,00 s | 120 nm | 0,12 pm | **3,48 m** | **7,0 µm** |
| 1.000.000 | 60 nm/s | 1,00 s | 60 nm | 0,06 pm | **6,95 m** | **13,9 µm** |
| 1.000.000 | 30 nm/s | 1,00 s | 30 nm | 0,03 pm | 13,90 m | 27,8 µm |
| 250.000 | 120 nm/s | 0,25 s | 30 nm | 0,12 pm | 3,48 m | 27,8 µm |
| 250.000 | 60 nm/s | 0,25 s | 15 nm | 0,06 pm | 6,95 m | 55,6 µm |
| *aktuell (Trigger)* | *10 nm/s* | *12 s* | *120 nm* | *1 pm* | *0,42 m* | *7,0 µm* |

**Lies die Tabelle so:**

- **Alle Zeilen mit 1M Punkten sind gut genug** für 2 m Testfaser. Die erste
  Zeile ist sogar identisch gut wie jetzt — nur mit **8× mehr Reichweite**.
- Der Sprung von 0,42 m auf 3,5 m löst dein Problem vollständig: das Faserende
  bei 1,05 m und alles bei 2 m liegt dann **innerhalb** der Reichweite und
  erscheint als scharfer, ehrlicher Peak statt als Alias-Dreck.
- **Bei 250k Punkten pro Kanal musst du Auflösung opfern** — deshalb ist die
  Frage deines Supervisors ("1M pro Kanal oder geteilt?") die **entscheidende**
  Frage vor dem Aufbau. Kläre die zuerst mit Giulio.

**Zum Puffer-Überlauf, den er erwähnt:** Er hat es leicht andersherum
formuliert, die Logik ist aber: der Puffer läuft nach `N × 1 µs` voll. Ist dein
Sweep **länger** als das, verlierst du das Ende (oder er überschreibt den
Anfang). Also **immer vorher rechnen**: `Sweep-Dauer = Span / Geschwindigkeit`
muss `< N × 1 µs` sein — plus etwas Reserve für Beschleunigung/Trigger-Latenz.

Deshalb steht in der Tabelle 60–120 nm/s statt der bisherigen 10 nm/s: mit
1 µs Takt und 1 s Puffer **musst** du schneller fahren, sonst passt der Span
nicht rein. Das ist übrigens akustisch unkritisch (die Warnung im Handover
betraf zu *langsames* Sweepen).

## B.5 Die wichtigste offene Zahl: wie lang muss der lange Arm sein?

Dein Supervisor sagt nur *"one arm should be longer than the other"* — ohne
Zahl. Die Zahl ist aber nicht beliebig, und sie ist **größer, als du denkst**.

**Zwei gegenläufige Bedingungen:**

**(1) Nach oben begrenzt** — die Aux-Streifen müssen bei 1 µs Takt sauber
abgetastet sein. Die Streifenfrequenz ist `f = τ_aux × dν/dt`. Bei 1 MHz
Abtastrate ist Nyquist 500 kHz, komfortabel sind ~7–10 Punkte pro Streifen.

**(2) Nach unten begrenzt** — Standardregel im OFDR: das Aux-Interferometer
sollte eine **mindestens so große Verzögerung haben wie der entfernteste
Reflektor, den du messen willst**, sonst korrigiert es das Abstimmrauschen für
weit entfernte Peaks nicht mehr vollständig (der Restfehler wächst grob mit
τ_mess/τ_aux).

**Die Zahlen:**

| | Wert |
|---|---|
| Testfaser 2 m → Laufzeit τ_mess | 19,6 ns |
| Nötige Aux-Armdifferenz (Regel: τ_aux ≥ τ_mess) | **ΔL ≈ 4 m Faser** |
| Aux-Streifenfrequenz bei 60 nm/s, ΔL = 4 m | 144 kHz → ~7 Punkte/Streifen ✔ |
| Aux-Streifenfrequenz bei 120 nm/s, ΔL = 4 m | 288 kHz → ~3,5 Punkte/Streifen (knapp) |
| Dein **bestehendes** MZI (92,8 mm) entspricht | z ≈ 46 mm — **~40× zu kurz** |

> **Das ist der Punkt, den du unbedingt mit deinem Supervisor klären musst:**
> "ein Arm länger" heißt hier **mehrere Meter Faser** (eine kleine Spule), nicht
> ein paar Zentimeter. Frag ihn, ob er das so meint, oder ob er auf
> Sub-Fringe-Interpolation setzt (dann darf das Aux kürzer sein, aber die
> Rauschunterdrückung wird schlechter).

**Empfohlener Startpunkt zum Diskutieren:** 1M Punkte/Kanal, Sweep 60 nm/s,
Span 60 nm, Aux-ΔL ≈ 4 m → Reichweite 6,95 m, Auflösung 13,9 µm, 7 Punkte pro
Aux-Streifen. Bequem in allen Dimensionen, und wenn es läuft, kannst du auf
120 nm/s / 120 nm hochgehen für die vollen 7 µm.

**Praktische Punkte zum Aux-MZI:**

- Der lange Arm sollte **thermisch ruhig** liegen (Spule in eine Schachtel,
  nicht frei über den Tisch). Seine Länge *ist* dein Maßstab — wenn sie driftet,
  driftet deine Distanzachse.
- In SM-Faser können die zwei Arme in der Polarisation auseinanderlaufen →
  Streifenkontrast bricht ein (**polarization fading**). Wenn das passiert:
  Polarisationsregler in einen Arm, oder mittelfristig PM-Faser (Handover §6,
  Punkt 3).
- **τ_aux einmal kalibrieren**: einen Sweep im *bisherigen* Triggermodus fahren,
  bei dem du die Wellenlängentabelle bekommst, und `gesamte entwickelte Phase /
  (2π × Frequenzspanne) = τ_aux` rechnen. Genau so wurden die 454,5 ps des
  bestehenden MZI bestimmt (Handover §6).

## B.6 Die Trigger-Umstellung

| | bisher | neu |
|---|---|---|
| EXFO | Trigger alle 1 pm | **Single-Trigger**: ein Puls am Sweep-Start |
| CoreDAQ | rot, 1 Sample pro Trigger, ~130k Puffer | **schwarz**: 1 Trigger → dann selbstständig alle 1 µs, 1M Puffer |
| ν-Achse kommt von | EXFO-Wellenlängentabelle | **Aux-MZI-Phase** |

Sein Argument, den schwarzen CoreDAQ zu nehmen statt den roten umzukonfigurieren,
ist pragmatisch: der schwarze **kann von Haus aus nur genau das**, was du
brauchst. Gerät tauschen ist einfacher als Firmware-Einstellungen kämpfen.
**Giulio kennt das Gerät** — frag ihn früh, nicht erst wenn es klemmt.

**Was du dabei verlierst und ersetzen musst:** ohne Trigger gibt es **keine
Wellenlängentabelle mehr**. Deine ν-Achse existiert nur noch, wenn das Aux
funktioniert. Deshalb: **Aux zuerst aufbauen und mit dem bisherigen Triggermodus
verifizieren**, erst dann auf Single-Trigger umstellen. Sonst debuggst du zwei
neue Dinge gleichzeitig.

Die absolute Startwellenlänge liest du weiterhin vom Laser ab — laut Handover
§6, Punkt 4 ist dessen absolute Genauigkeit ~100× besser als nötig.

## B.7 Deine Checkliste

**Zuerst — kostet nichts, klärt viel (aus Handover §7):**

1. **Maßband auf den Tisch.** Ist der 46,4-mm-Peak eine 92,8-mm-Armdifferenz
   oder ein Reflektor bei 46,4 mm? Ungeklärt, blockiert das Verständnis der
   Topologie.
2. **Faserende terminieren** (Index-Matching-Gel, enge Wicklung, oder ein
   < 35 cm Patchcord) und mit `ofdr_diagnose.py compare` vorher/nachher
   verifizieren. Vorhersage: das 6–32-cm-Band fällt ins Rauschen. Damit
   bestätigst du die Aliasing-Diagnose mit eigenen Händen und hast sofort
   saubere Daten.

**Vor dem Aufbau klären:**

3. **1M Punkte pro Kanal oder auf 4 Kanäle geteilt?** (→ Giulio). Entscheidet
   über die ganze Tabelle in B.4.
4. **Welche maximale Sweep-Geschwindigkeit kann der EXFO?** Du brauchst
   voraussichtlich 60–120 nm/s.
5. **Wie lang soll der lange Aux-Arm sein?** (→ Supervisor, mit den Zahlen aus
   B.5 in der Hand). Meine Rechnung sagt ~4 m.
6. Hat der schwarze CoreDAQ wirklich 4 Kanäle und feste Messbereiche?
   **Autoranging mitten im Sweep ist tödlich** (Phasensprung → Unwrap kaputt
   bis zum Ende des Scans, Handover §7).

**Beim Aufbau:**

7. Aux-MZI bauen, τ_aux **im alten Triggermodus** kalibrieren.
8. Erst dann auf Single-Trigger + Freilauf umstellen.
9. Vor jeder Pipeline-Änderung: an synthetischen Daten mit bekannter Wahrheit
   testen (`make_test_data.py`).
10. Nach jedem Scan: **Peak-Breite prüfen.** ≤ ~21 µm bei Kaiser-12, sonst
    stoppen und untersuchen.

---

## Wenn du eine Sache aus diesem Dokument mitnimmst

```
Reichweite / Auflösung = N/2
```

Die Auflösung kommt aus der **Sweep-Breite**. Die Reichweite kommt aus dem
**Punktabstand**. Alles, was dein Supervisor vorschlägt, dient genau einem
Zweck: mehr Punkte in denselben Sweep zu bekommen — und das Aux-MZI ist der
Preis, den man dafür zahlt, dass der Laser nicht mehr sagt, wo er gerade ist.
