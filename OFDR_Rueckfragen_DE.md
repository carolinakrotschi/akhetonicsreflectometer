# OFDR — Antworten auf deine Rückfragen

**Für:** Carolina · **Stand:** 18.08.2026 · Ergänzung zu `OFDR_Erklaerung_DE.md`

Mitgeliefert (im selben Ordner):

- `ofdr_sweep_calc.py` — rechnet dir jede Parameterkombination durch, inkl. Aux-Prüfung
- `make_aux_test_data.py` — erzeugt synthetische 4-Kanal-Daten mit bekannter Wahrheit
- `ofdr_aux.py` — die 4-Kanal-Pipeline (Antwort auf Frage 10, läuft schon)
- `vergleich_alt_neu.png` — alt vs. neu, an simulierten Daten

---

## 1. Warum kann man das `/2` im Skript nicht einfach weglassen?

Weil das `/2` keine Einstellung ist, sondern die **Definition der Achse**. Und
weil auf einer Achse nur eine Definition gelten kann, während auf dem Detektor
gleichzeitig beide Sorten Signal liegen.

Was die FFT dir gibt, ist immer nur **τ**, die Laufzeit. Die ist eindeutig, da
gibt es nichts zu entscheiden. Erst der Schritt von τ zu einer Entfernung
braucht eine Annahme, nämlich *"wie oft ist das Licht diese Strecke gelaufen?"*:

| Signalart | Weg | Umrechnung |
|---|---|---|
| **Reflexion** (Stecker, Faserende) | hin **und** zurück | z = c·τ/(2·n_g) |
| **Transmission** (MZI-Arm) | nur hin | ΔL = c·τ/n_g |

In einem Reflektometer sind die *meisten* Signale Reflexionen. Also ist `/2`
richtig. Würdest du es entfernen, würden **alle echten Reflektoren doppelt so
weit** angezeigt — du hättest ein Problem gegen ein größeres getauscht.

Das MZI-Signal ist die eine Ausnahme, und Ausnahmen kann man nicht durch
Ändern der Achse behandeln, sondern nur durch Wissen. Deshalb macht das Skript
genau das Richtige: es **druckt beide Lesarten** für den Hauptpeak aus
(`ofdr_process.py`, Zeile 222–224):

```
delay 454.50 ps | as transmission-MZI imbalance: 92.82 mm fiber
                | as reflection: 46.41 mm fiber one-way
```

Die Information ist also da. Was fehlt, ist die Messung mit dem Maßband, die
sagt, welche der beiden Zeilen stimmt.

**Und die gute Nachricht:** Sobald das Aux-MZI auf eigenen Kanälen (Ch3/Ch4)
sitzt, ist das Problem strukturell weg. Dann liegt das Transmissionssignal
gar nicht mehr auf den Messkanälen, und die Achse auf Ch1/Ch2 ist eindeutig
eine Reflexionsachse. Das ist ein weiteres, oft übersehenes Argument für den
Umbau.

---

## 2. Woher kommen die 122 MHz bei 1 pm Schritt?

Reine Umrechnung. Zwei Wege, beide gleich:

**Weg A — differenzieren.** ν = c/λ, also

```
|dν/dλ| = c/λ²

c/λ²  bei λ = 1565 nm  =  2,998·10⁸ / (1,565·10⁻⁶)²
                       =  2,998·10⁸ / 2,449·10⁻¹²
                       =  1,224·10²⁰  Hz pro Meter
```

Das ist der Umrechnungsfaktor. Jetzt mit dem Schritt multiplizieren:

```
δν = 1,224·10²⁰ Hz/m  ×  1 pm
   = 1,224·10²⁰  ×  1·10⁻¹² m
   = 1,224·10⁸ Hz
   = 122 MHz
```

**Weg B — über relative Änderungen (finde ich anschaulicher).** Bei ν = c/λ
ist die *relative* Änderung in λ und in ν gleich groß (nur mit Vorzeichen):

```
1 pm von 1565 nm         =  1·10⁻¹² / 1,565·10⁻⁶  =  6,39·10⁻⁷   (also 0,64 ppm)
optische Frequenz        =  c/1565 nm  =  191,6 THz
0,64 ppm von 191,6 THz   =  122 MHz  ✔
```

Merkhilfe: **1 pm ≈ 122 MHz bei 1550 nm.** Und weil λ² im Nenner steht, ist der
Faktor am blauen Ende größer: 132 MHz bei 1505 nm, 114 MHz bei 1625 nm. Genau
das sind die 18 % Chirp.

---

## 3. Was bedeutet "der Reflektor bei τ erzeugt eine Schwingung mit Periode 1/τ"?

Es bedeutet: **um wie viel muss der Laser seine Frequenz verstimmen, damit das
Signal einmal von hell nach dunkel nach hell durchläuft?**

Aus der Grundgleichung: die Phase ist `2π·τ·ν`. Eine volle Periode ist erreicht,
wenn sich diese Phase um 2π ändert:

```
2π · τ · Δν = 2π      →      Δν = 1/τ
```

Die Einheit passt: τ in Sekunden → 1/τ in Hertz. Es ist also eine Periode
**gemessen in Hertz optischer Frequenz** — nicht in Sekunden. Ungewohnt, aber
genau darum geht es: wir arbeiten im Frequenzbereich, ν ist unsere x-Achse.

**Konkret an deinem Setup:**

| Reflektor | τ | 1/τ | in Wellenlänge | Punkte pro Fringe bei 1 pm |
|---|---|---|---|---|
| MZI-Peak (46,4 mm) | 454,5 ps | 2,20 GHz | **18,0 pm** | **18,0** ✔ |
| Faserende 1,05 m | 10,28 ns | 97,3 MHz | **0,79 pm** | **0,79** ✘ |
| Bauteil bei 2,0 m | 19,59 ns | 51,1 MHz | **0,42 pm** | **0,42** ✘ |

(Das Handover nennt "17 pts/fringe" für das MZI — passt zu meinen 18, die
Differenz ist die Wahl von λ in der Umrechnung.)

**Lies die letzte Spalte:** Der Aufbau kann das MZI-Signal problemlos abtasten
(18 Punkte pro Streifen, üppig). Das Faserende bei 1,05 m bekommt **weniger als
einen** Punkt pro Streifen. Nyquist verlangt mindestens 2. Deshalb aliasiert es
— und deshalb ist die Grenze bei 0,42 m, wo genau 2 Punkte pro Streifen
herauskommen.

**Anschaulich:** weit weg = schnelle Streifen. Zu weit weg = Streifen schneller
als deine Abtastung = du siehst eine falsche, langsamere Schwingung. Wagenrad
im Film.

---

## 4. Was bedeuten die 59 pm Spitze-Spitze?

Du siehst den Bogen — die Frage ist, was daran schlimm ist. Antwort in drei
Schritten.

**Schritt 1: Was die Zahl beschreibt.** Der Laser behauptet, bei Punkt Nummer k
sei er bei `λ_start + k · 1 pm`. Die rechte Hälfte von `real_data_diagnosis.png`
zeigt, wo er **wirklich** ist: bis zu +18 pm darüber (in der Sweep-Mitte) und
bis zu −40 pm darunter (am Ende). Gesamter Ausschlag ≈ 59 pm.

Deine x-Achse ist also falsch. Nicht zufällig verrauscht, sondern glatt und
systematisch verbogen.

**Schritt 2: Warum 59 pm von 120.000 pm trotzdem viel ist.** Relativ sind das
0,05 % — klingt harmlos. Aber die FFT interessiert sich nicht für relative
Fehler in der Achse, sondern für **Phasenfehler**. Und der ist groß:

```
59 pm  ×  122 MHz/pm            =  7,2 GHz Frequenzfehler
Phasenfehler auf dem 454,5-ps-Ton:
  Δφ = 2π · τ · Δν
     = 2π · 454,5 ps · 7,2 GHz
     = 20,6 rad
     = 3,3 volle Perioden
```

Der Ton, der eine perfekt gleichmäßige Schwingung sein *sollte*, läuft also über
den Sweep um **3,3 ganze Perioden** aus dem Takt. Für die FFT ist das nicht mehr
eine Frequenz, sondern ein Frequenzgemisch → der Peak zerfasert.

**Schritt 3: Warum daraus ±0,8 mm werden.** Ein Peak sitzt dort, wo die
*Steigung* der Phase liegt (Phase pro Frequenz = Verzögerung). Weil der Fehler
ein glatter Bogen ist, hat er in der Mitte eine andere Steigung als am Rand:

```
Bogenamplitude A ≈ 10 rad, Sweep-Spanne S = 14,7 THz
maximale Steigungsabweichung  dτ ≈ A·π/S  ≈ 2,2 ps
                              dz  = c·dτ/(2·n_g)  ≈ 0,22 mm
```

In derselben Größenordnung wie die im Handover genannten ±0,8 mm (mein
Bogenmodell ist eine glatte Sinushälfte; die echte Kurve ist kantiger und
liefert entsprechend mehr).

**Vergleich zur Sollbreite: 21 µm.** Ohne Korrektur wird der Peak also um einen
Faktor ~30–40 breiter, als er sein müsste. Das ist der ganze Punkt.

**Und das ist der Grund, warum das Aux-MZI keine Verfeinerung, sondern eine
Notwendigkeit ist.** Aktuell wird der Bogen aus dem Messsignal selbst
geschätzt, was funktioniert, aber schwache echte Merkmale mit-verbiegt (daher
`diagnostic` vs. `cosmetic`). Das Aux misst ihn **direkt und unabhängig**.

---

## 5. Wofür braucht man `cosmetic`?

Kurz: **nicht zum Messen. Als Obergrenze und als Test, ob die Achse schuld
ist.** Drei Verwendungen:

**(a) Beweisen, dass die Hardware es kann.** Wenn du in `cosmetic` einen 21-µm-
Peak bekommst, dann *ist* der Reflektor physikalisch scharf und deine Optik
liefert die volle Kohärenz. Das ist eine Aussage über den Aufbau, nicht über
die Reflektorpositionen.

**(b) Der eigentliche diagnostische Wert: "ist die Frequenzachse das Problem?"**

| `diagnostic` | `cosmetic` | Diagnose |
|---|---|---|
| breit | **scharf** | Die Achse ist schuld. Bessere Korrektur → Aux-MZI |
| breit | **auch breit** | Die Achse ist *nicht* schuld. Suche in der Physik: Dispersion, Polarisation, oder der Reflektor ist wirklich ausgedehnt |

Das ist ein echter Diskriminator, in derselben Familie wie die vier Tests in
`ofdr_diagnose.py`. Er verhindert, dass du tagelang die Software optimierst,
während das Problem auf dem Tisch liegt — oder umgekehrt.

**(c) Bilder für den Supervisor.** "Wir erreichen 21 µm, 15 % über dem
Transformlimit" ist eine legitime Kennzahl, und sie kommt aus diesem Modus.

**Warum du damit nicht diagnostizieren darfst:** Die Korrektur folgt der
*Rohphase* des Hauptsignals. Alles, was zusätzlich im Signal steckt — also
genau die schwachen Reflektoren, die du suchst — verbiegt diese Phase mit und
wird von der Korrektur teilweise **weggerechnet**. Gemessen wurden −8 dB
Verlust an echten Merkmalen plus Positionsverschiebungen. Der Modus macht das
Bild schöner, indem er löscht, was dich interessiert.

**Und:** nach dem Umbau ist `cosmetic` überflüssig. Die Aux-Referenzierung gibt
dir gleichzeitig scharfe Peaks *und* korrekte Positionen, weil die Korrektur
nicht mehr aus dem Messsignal kommt. Der ganze Modus ist ein Symptom des
aktuellen Provisoriums.

---

## 6. Warum ist das Signal im PNG im negativen dB-Bereich?

Weil die Skala **relativ zum größten Peak** ist, nicht absolut. Im Code
(`ofdr_process.py`, Zeile 204):

```python
db = 20 * np.log10(R / R.max() + 1e-15)
```

Es wird durch `R.max()` geteilt. Der stärkste Peak ist damit per Definition
`20·log10(1) = 0 dB`, und alles andere ist kleiner als 1 → Logarithmus negativ.
Es gibt in diesem Diagramm gar keine positiven Werte, das ist strukturell so.

**Was die Zahlen bedeuten** (beachte: `20·log10`, weil `R` eine *Amplitude* ist,
nicht eine Leistung):

| dB | Amplitude relativ zum Hauptpeak |
|---|---|
| 0 | 1 : 1 (das ist der Hauptpeak selbst) |
| −21 | 1 : 11 |
| −45 | 1 : 178 |
| −53 | 1 : 447 |
| −80 | 1 : 10 000 |

Also: das Artefaktband bei −45 dB ist ~180× schwächer als der Hauptpeak. Die
Kaiser-Nebenkeulen bei −80 dB sind 10.000× schwächer — deshalb ist Kaiser die
richtige Wahl, wenn man schwache echte Reflexionen sucht.

**Warum man es so macht:** dich interessiert das **Verhältnis** (ist dieser
Reflektor stark oder schwach *im Vergleich* zum Hauptsignal), nicht die
absolute Leistung in mW. Für absolute Werte (Return Loss in dB) bräuchtest du
eine Kalibrierung mit einer bekannten Referenzreflexion — die habt ihr noch
nicht, steht aber auch nicht auf der Liste.

---

## 7. 250k Punkte pro Kanal mit 120 nm Span — geht das?

**Rechnerisch ja, praktisch nein.** Und der Grund ist interessant. Erst der
Rechenweg, dann die Werte.

### Der Rechenweg (immer derselbe, fünf Zeilen)

Fest ist der Takt: 1 Punkt alle 1 µs. Daraus folgt alles:

```
1)  T        = N × 1 µs                     Aufnahmedauer
2)  v        = Span / T                     nötige Sweep-Geschwindigkeit
3)  δλ       = v × 1 µs                     Punktabstand in Wellenlänge
4)  δν       = δλ × 122 MHz/pm              Punktabstand in Frequenz
5a) z_max    = c / (4·n_g·δν)               REICHWEITE
5b) Δz       = c / (2·n_g·Δν_span)          AUFLÖSUNG
```

### Eingesetzt für 250.000 Punkte und 120 nm

```
1)  T   = 250 000 × 1 µs                    = 0,25 s
2)  v   = 120 nm / 0,25 s                   = 480 nm/s
3)  δλ  = 480 nm/s × 1 µs                   = 0,48 pm
4)  δν  = 0,48 × 122,4 MHz                  = 58,8 MHz
5a) z_max = 3·10⁸ / (4 × 1,468 × 58,8·10⁶)  = 0,87 m
5b) Δν_span für 120 nm = 14,7 THz
    Δz  = 3·10⁸ / (2 × 1,468 × 14,7·10¹²)   = 7,0 µm
```

Kontrolle: 0,87 m ÷ 7,0 µm = 125.000 = N/2 ✔

### Und jetzt die drei Gründe, warum das trotzdem nicht geht

1. **0,87 m Reichweite reicht nicht** für 2 m Testfaser. Besser als die
   jetzigen 0,42 m, aber das Bauteil bei 2 m aliasiert weiter. Problem nicht
   gelöst, nur verkleinert.
2. **480 nm/s** — kann der EXFO das überhaupt? Musst du nachschauen. Viele
   abstimmbare Laser hören bei 100–200 nm/s auf.
3. **Der Killer: das Aux-MZI passt nicht mehr rein.** Bei 480 nm/s und ΔL = 4 m
   schlägt das Aux mit **1,15 MHz**. Deine Abtastrate ist 1 MHz, Nyquist also
   500 kHz. **Das Lineal selbst aliasiert.** Damit ist die Frequenzachse
   kaputt, und ohne Frequenzachse gibt es überhaupt kein Ergebnis.

Das Skript sagt dir das direkt:

```
$ python ofdr_sweep_calc.py --points 250000 --span 120 --dl 4
  REICHWEITE  z_max = 0.869 m
  AUFLOESUNG  dz    = 6.95 um
  --- Aux-MZI mit dL = 4.00 m ---
  erscheint bei       z = 2.000 m (muss < z_max = 0.869 m sein)
  Schwebungsfrequenz  1150.8 kHz  = 0.9 Punkte pro Fringe
    FEHLER: Aux liegt jenseits Nyquist -- es aliasiert selbst.
    FEHLER: unter 2 Punkte pro Aux-Fringe -- Phase nicht rekonstruierbar.
```

### Was mit 250k stattdessen geht

| Span | v | z_max | Δz | Aux ΔL=4 m |
|---|---|---|---|---|
| 120 nm | 480 nm/s | 0,87 m | 7,0 µm | ✘ aliasiert |
| 60 nm | 240 nm/s | 1,74 m | 13,9 µm | ✘ aliasiert |
| 30 nm | 120 nm/s | 3,48 m | 27,8 µm | grenzwertig (3,5 Pkt/Fringe) |
| **15 nm** | **60 nm/s** | **6,95 m** | **55,6 µm** | ✔ 7,0 Pkt/Fringe |

Vergleich mit 1M Punkten:

| Span | v | z_max | Δz | Aux ΔL=4 m |
|---|---|---|---|---|
| 120 nm | 120 nm/s | 3,48 m | **7,0 µm** | grenzwertig (3,5 Pkt/Fringe) |
| **60 nm** | **60 nm/s** | **6,95 m** | **13,9 µm** | ✔ 7,0 Pkt/Fringe |

**Fazit: Der Pufferunterschied ist ein Faktor 4 in der Auflösung** (55,6 µm vs.
13,9 µm bei gleichem Komfort). Deshalb ist "1M pro Kanal oder geteilt?" die
erste Frage an Giulio, noch vor allem Löten.

Spiel selbst damit:

```bash
python ofdr_sweep_calc.py --table
python ofdr_sweep_calc.py --points 1000000 --span 60 --dl 4
```

---

## 8. Wie kommt man auf 4 m Armdifferenz? Schritt für Schritt

### Untergrenze — Genauigkeit

**Schritt 1.** Was ist das Entfernteste, das du messen willst? → 2 m Testfaser.

**Schritt 2.** Dessen Laufzeit (hin **und** zurück, es ist eine Reflexion):

```
τ_mess = 2 · n_g · z / c = 2 · 1,468 · 2 m / 2,998·10⁸ = 19,59 ns
```

**Schritt 3.** Die Regel: **τ_aux ≥ τ_mess.**

Warum? Das Aux ist dein Lineal, und seine Teilstriche sind seine Streifen. Der
Abstand der Teilstriche in Frequenz ist `1/τ_aux` (Frage 3). Der zu messende
Reflektor schwingt mit `1/τ_mess`. Ein Lineal, dessen Teilstriche gröber sind
als das, was du messen willst, kann die Verzerrung dort nicht mehr auflösen —
der Restfehler wächst grob mit `τ_mess/τ_aux`.

**Schritt 4.** Umrechnen. Das Aux ist ein **Transmissions**-MZI, das Licht läuft
nur einmal durch den längeren Arm — **kein Faktor 2**:

```
τ_aux = n_g · ΔL / c        →      ΔL = c · τ_aux / n_g
                                      = 2,998·10⁸ · 19,59 ns / 1,468
                                      = 4,00 m
```

**Die Abkürzung, die man sich merken kann:**

```
τ_aux = τ_mess
n_g·ΔL/c = 2·n_g·z/c
        ΔL = 2 · z_max
```

> **ΔL = 2 × der weiteste Reflektor.** 2 m messen → 4 m Armdifferenz.
> 3 m messen → 6 m. So einfach ist es.

Äquivalent und praktischer zu prüfen: das Aux **erscheint im Reflektogramm bei
z = ΔL/2**. Diese Position muss mindestens bei deinem entferntesten Reflektor
liegen. ΔL = 4 m → Aux-Peak bei 2,0 m ✔

### Obergrenze — Abtastung

Das Aux ist selbst ein Signal auf dem Detektor. Es muss also **innerhalb der
Nyquist-Reichweite** liegen:

```
ΔL/2  <  z_max          also        ΔL < 2 · z_max
```

Und in Frequenz ausgedrückt — die Streifenrate des Aux:

```
dν/dt  = (c/λ²) · dλ/dt = 122,4 MHz/pm · Sweep-Geschwindigkeit
f_aux  = τ_aux · dν/dt
Punkte pro Fringe = 1 MHz / f_aux
```

| ΔL | Sweep | f_aux | Punkte/Fringe | Urteil |
|---|---|---|---|---|
| 4 m | 30 nm/s | 72 kHz | 13,9 | sehr komfortabel |
| **4 m** | **60 nm/s** | **144 kHz** | **7,0** | **gut** |
| 4 m | 120 nm/s | 288 kHz | 3,5 | über Nyquist, aber knapp |
| 4 m | 480 nm/s | 1151 kHz | 0,9 | kaputt |
| 8 m | 60 nm/s | 288 kHz | 3,5 | knapp |

### Also, zusammengesetzt

```
2 · z_mess   ≤   ΔL   <   2 · z_max
   4 m       ≤   ΔL   <   13,9 m        (bei 1M Punkten, 60 nm/s)
```

**ΔL = 4 m ist die Untergrenze und gleichzeitig eine gute Wahl** — genug Lineal
für 2 m Messstrecke, 7 Punkte pro Streifen, und viel Luft nach oben, falls die
Testfaser später länger wird.

**Wenn dein Supervisor konservativer sein will:** In der Literatur findet man
oft `τ_aux ≥ 2·τ_mess` (also ΔL = 8 m). Das gilt, wenn man **hart auf
Streifenflanken taktet** (ein Messpunkt pro Aux-Streifen) — dann ist das Aux
selbst der Abtaster und muss Nyquist für das Messsignal erfüllen. Wir tasten
aber mit 7 Punkten pro Streifen und interpolieren die Phase, damit ist ΔL = 4 m
ausreichend. **Kläre mit ihm, welche der beiden Varianten er meint** — es macht
den Unterschied zwischen einer 4-m- und einer 8-m-Spule.

---

## 9. Muss ich für die 4 Kanäle neue Software schreiben?

**Nein — und ich habe sie dir schon geschrieben und getestet.** `ofdr_aux.py`
liegt im Ordner und funktioniert. Aber lies, was sich ändert, damit du sie
verstehst und anpassen kannst.

### Was sich am Code ändert

| | |
|---|---|
| **Unverändert** | balancierte Subtraktion (jetzt 2× aufgerufen), Fensterung, FFT, Peak-Liste, Breitenprüfung, Plot |
| **Fällt weg** | `to_uniform_nu()` — es gibt keine Wellenlängentabelle mehr |
| **Fällt weg** | `phase_correct()` mit `diagnostic`/`cosmetic` — keine Selbstreferenzierung mehr, damit verschwindet der ganze Konflikt |
| **Neu** | `resample_on_aux()` — ~25 Zeilen |

Netto also **weniger** Code als jetzt. Es ist eine Änderung, kein Neuanfang.

### Was `resample_on_aux()` tut, in vier Zeilen Logik

```python
phi = np.unwrap(np.angle(analytic(aux)))      # Aux-Phase, NICHT geglättet
phi = np.maximum.accumulate(phi)              # Monotonie erzwingen
phi_u = np.linspace(phi[0], phi[-1], m)       # gleichmäßiges Phasenraster
y = PchipInterpolator(phi, meas)(phi_u)       # Messsignal darauf umtasten
```

Der Kern in einem Satz: die Aux-Phase ist `φ(t) = 2π·τ_aux·ν(t) + const`.
Gleiche Schritte in φ sind also **gleiche Schritte in ν** — völlig egal, wie
ungleichmäßig der Laser tatsächlich gefahren ist. Danach ist die Achse fertig
und die FFT erlaubt.

Und `δν = Δφ / (2π·τ_aux)` gibt dir Reichweite und Auflösung, ohne dass du
irgendwo eine Wellenlänge gebraucht hättest.

### Der Beweis, dass es funktioniert

Ich habe synthetische 4-Kanal-Daten mit **bekannter Wahrheit** erzeugt
(1 Mio. Punkte, 60 nm/s, Aux ΔL = 4 m, mit dem gemessenen 59-pm-Bogen und
Abstimmwelligkeit drin) und durch die Pipeline geschickt:

```
   z_wahr  0.0460 m (  0.0 dB) -> gefunden  0.0460 m, Fehler  +4.1 um (+0.30 Zellen)
   z_wahr  0.3000 m (-25.0 dB) -> gefunden  0.3000 m, Fehler  -4.8 um (-0.35 Zellen)
   z_wahr  1.0500 m (-30.0 dB) -> gefunden  1.0500 m, Fehler  +3.6 um (+0.26 Zellen)
   z_wahr  2.0000 m (-35.0 dB) -> gefunden  2.0000 m, Fehler  +4.2 um (+0.31 Zellen)
   groesster Fehler: 4.8 um (eine Zelle = 13.6 um)   -> BESTANDEN

Hauptpeak 46.0041 mm, -3 dB Breite 40.9 um (Fenstergrenze ~35.5 um)
```

Alle vier Reflektoren auf besser als eine halbe Auflösungszelle, Amplituden
innerhalb 0,5 dB, Peakbreite 15 % über dem Transformlimit — **dieselbe
Gesundheitskennzahl wie jetzt, aber mit 6,7 m statt 0,42 m Reichweite.**

`vergleich_alt_neu.png` zeigt beides übereinander: oben derselbe Reflektorsatz
im alten Triggermodus (1,05 m und 2,0 m falten nach 0,22 m und 0,33 m und
füllen den ganzen Bereich mit Dreck), unten aux-referenziert (vier saubere
Peaks). Das ist eine Simulation, keine Messung — aber es ist genau das
Artefakt aus §5 des Handovers, reproduziert und dann behoben.

### Was du selbst noch tun musst

1. **Das Dateiformat des schwarzen CoreDAQ anpassen.** `ofdr_aux.py` liest
   `.npz` und die bekannte JSON-Struktur mit `Ch1..Ch4 [mW]`. Was das schwarze
   Gerät wirklich exportiert, weißt du erst, wenn du es angeschlossen hast
   (→ Giulio). Es ist die Funktion `load()`, ganz oben, ~15 Zeilen.
2. **Prüfen, welche Kanäle Aux und welche Messung sind.** Steht als Annahme im
   Kopf der Datei: Ch1/Ch2 Messung, Ch3/Ch4 Aux.
3. **τ_aux kalibrieren** und mit `--tau-aux-ns` übergeben. Ohne das schätzt das
   Skript aus `--dl`, was für Positionen auf ~1 % genau ist — gut für erste
   Bilder, nicht für Zahlen im Bericht.

### Ein Nebenbefund, den du kennen solltest

Beim Simulieren ist mir aufgefallen: die im Handover genannte schnelle
Abstimmwelligkeit von **90–110 MHz rms bei ~4,5 pm Periode liegt gefährlich nah
an einer harten Grenze.** Wenn die Steigung dieser Welligkeit die nominelle
Sweep-Steigung übersteigt, fährt der Laser momentan **rückwärts** in ν. Dann ist
die Aux-Phase nicht mehr monoton, und **kein** Verfahren kann das noch
auflösen — auch das Aux nicht.

```
Grenze:  A_max = (c/λ²) · Periode / (2π)
              = 122,4 MHz/pm · 4,5 pm / (2π)
              ≈ 88 MHz
```

Das Handover selbst warnt, dass diese Zahl "durch den Hauptton gemessen" und
damit von additiven Störern kontaminiert ist — sie ist also wahrscheinlich
überschätzt. **Aber es ist eine Annahme, auf der der ganze Umbau ruht.** Das
Gute: das Aux-MZI misst diese Welligkeit sauber und direkt. Also: sobald das Aux
steht, **als erstes die Welligkeit nachmessen** — noch bevor du auf Freilauf
umstellst. `make_aux_test_data.py --ripple-mhz 100` zeigt dir, wie das
Fehlerbild aussieht, wenn es doch zu viel ist.

---

## 10. Was du jetzt genau machen sollst

### Phase 0 — heute/morgen, ohne einzuschrauben

**0.1** Maßband an den 46,4-mm-Peak. Ist es ein 92,8-mm-Armunterschied oder ein
Reflektor bei 46,4 mm? Zähle die Fasern und Steckerlängen nach.
→ *Ergebnis in ein Laborbuch, das ist eine seit dem Handover offene Frage.*

**0.2** Referenzscan im aktuellen Zustand aufnehmen und ablegen. Dann:

```bash
python ofdr_diagnose.py tuning  scan_heute.json
python ofdr_process.py          scan_heute.json --mode diagnostic
```

Notiere: Peakbreite (soll ≤ ~21 µm), Bogen in pm, schnelle Welligkeit in MHz.
→ *Das ist deine Grundlinie. Ohne sie kannst du später nicht sagen, ob etwas
besser geworden ist.*

**0.3** Faserende terminieren (Index-Matching-Gel auf die APC-Endfläche, oder
enge Wicklung über die letzten Zentimeter, oder ein Patchcord < 35 cm) und
neuen Scan aufnehmen. Dann:

```bash
python ofdr_diagnose.py compare scan_heute.json scan_terminiert.json
```

Erwartung: das Band bei 6–32 cm fällt ins Rauschen.
→ *Damit bestätigst du die Aliasing-Diagnose selbst, und du hast ab sofort
saubere Daten zum Arbeiten. Höchster Nutzen pro Aufwand von allem auf dieser
Liste.*

**0.4** Zum Aufwärmen: Simulation laufen lassen, damit du die Pipeline kennst,
bevor du Hardware anfasst.

```bash
python make_aux_test_data.py --points 400000 --speed 60 --dl 4 --out sim.npz
python ofdr_aux.py sim.npz --dl 4 --zmax 2.5
```

### Phase 1 — Fragen klären, bevor etwas gekauft/gelötet wird

**1.1 → Giulio:** Hat der schwarze CoreDAQ **1M Punkte pro Kanal oder 1M
geteilt durch 4**? *(Faktor 4 in der Auflösung. Wichtigste Zahl.)*

**1.2 → Giulio:** 4 Kanäle gleichzeitig? Feste Messbereiche einstellbar?
*(Autoranging mitten im Sweep zerstört den Scan.)* Welches Exportformat?

**1.3 → EXFO-Handbuch:** maximale Sweep-Geschwindigkeit? *(Du brauchst
voraussichtlich 60–120 nm/s.)* Und: wie schaltet man Single-Trigger ein?

**1.4 → Supervisor:** ΔL = 4 m oder 8 m? Zeig ihm Abschnitt 8 mit der Rechnung.
*(Er hat keine Zahl genannt, und sein bestehendes MZI mit 92,8 mm ist ~40×
zu kurz.)*

**1.5 → Supervisor:** Startkonfiguration bestätigen. Mein Vorschlag:

| | |
|---|---|
| Punkte | 1M pro Kanal |
| Sweep | 60 nm/s |
| Span | 60 nm |
| Aux ΔL | 4 m |
| ⇒ Reichweite | 6,95 m |
| ⇒ Auflösung | 13,9 µm |
| ⇒ Aux | 7 Punkte pro Streifen |

*Erst wenn das läuft, auf 120 nm/s und 120 nm Span für die vollen 7 µm.*

### Phase 2 — Aufbau, in dieser Reihenfolge

**2.1** 90/10-Splitter direkt nach dem Laser. 90 % in den bestehenden Aufbau.
*Vorher/nachher-Scan: die 10 % Verlust sollen nichts kaputt machen.*

**2.2** Aux-MZI bauen: 50/50 → zwei Arme mit ΔL ≈ 4 m → 50/50 zurück, beide
Ausgänge auf Ch3/Ch4. Lange Spule thermisch ruhig legen (in eine Schachtel,
nicht frei über den Tisch — ihre Länge *ist* dein Maßstab).

**2.3** Aux **im alten Triggermodus** prüfen: sauberer Sinus? Kontrast stabil
über den ganzen Sweep? *(Bricht er irgendwo ein → Polarisationsfading →
Polarisationsregler in einen Arm.)*

**2.4** τ_aux kalibrieren, im alten Triggermodus mit Wellenlängentabelle:

```
τ_aux = gesamte entwickelte Aux-Phase / (2π · Frequenzspanne)
```

**2.5** Die schnelle Welligkeit mit dem Aux nachmessen (Abschnitt 9,
Nebenbefund). Ist sie unter ~50–80 MHz rms → grünes Licht.

**2.6** Erst **jetzt** auf Single-Trigger + Freilauf umstellen. Ein neues Ding
nach dem anderen, sonst debuggst du zwei Unbekannte gleichzeitig.

**2.7** Puffer-Kontrolle bei jedem Scan: `Span / Geschwindigkeit < N × 1 µs`.
Wenn die Streifen im letzten Teil des Puffers abbrechen, war der Sweep zu lang.

**2.8** Erste echte Messung durch `ofdr_aux.py`. Prüfen:
Peakbreite ≲ 1,5× Fenstergrenze, Aux ≥ 4 Punkte/Streifen, und — das ist der
Moment der Wahrheit — **das Faserende bei ~1 m als scharfer Peak an seiner
wahren Position, nicht mehr als Band bei 21 cm.**

### Dauerregeln

- Nach jedem Scan: **Peakbreite prüfen.** Regression → stoppen und untersuchen.
- Vor jeder Pipeline-Änderung: **an synthetischen Daten mit bekannter Wahrheit
  testen.** Du hast jetzt den Generator dafür.
- Jede Hypothese über ein unerwartetes Merkmal bekommt eine **Messung**, keine
  Geschichte. Die vier Tests in `ofdr_diagnose.py` sind je ein Kommando. Das
  Handover hält fest: drei plausible Theorien waren falsch, jede wurde von einer
  Messung getötet.
