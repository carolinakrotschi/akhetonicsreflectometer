# Ablaufplan — welches Programm, wann, mit welchem Aufbau

**Für:** Carolina · 18.08.2026

---

## ⚠ Zuerst: eine Korrektur an meinem eigenen Plan

In der Liste, die du zitierst, standen die Schritte 2.3–2.5 als *"Aux im alten
Triggermodus prüfen und kalibrieren"*. **Das geht nicht.** Ich habe es
nachgerechnet:

```
Alter Triggermodus, 1 pm Schritt  ->  Reichweite 0,42 m
Ein Aux mit 4 m Armdifferenz      ->  erscheint bei 2,00 m
```

Das Aux läge weit außerhalb der Reichweite und würde **selbst aliasieren**. Im
Triggermodus ließe sich höchstens ein Aux mit **ΔL < 0,83 m** prüfen — für
2 m Messstrecke brauchst du aber 4 m.

**Die Lösung ist einfacher als der ursprüngliche Plan:** τ_aux lässt sich im
Freilauf durch **Fringe-Zählen** kalibrieren. Du brauchst dafür gar keine
Wellenlängentabelle, nur die Start- und Endwellenlänge, die du am Laser
einstellst:

```
τ_aux = Anzahl gezählter Fringes / Frequenzspanne des Sweeps
```

Genauigkeit: ~7 ppm. Völlig ausreichend. Dafür gibt es jetzt ein eigenes
Programm, `aux_check.py`.

Die Reihenfolge unten ist entsprechend korrigiert.

---

## Die sechs Programme auf einen Blick

| Programm | Messgerät | Trigger | Kanäle | wofür |
|---|---|---|---|---|
| `ofdr_process.py` | rot | alter Modus | Ch1, Ch2 | die bisherige Auswertung |
| `ofdr_diagnose.py` | rot | alter Modus | Ch1, Ch2 | die vier Diagnosetests |
| `ofdr_sweep_calc.py` | — | — | — | Planung am Schreibtisch |
| `make_aux_test_data.py` | — | — | — | Simulation zum Üben und Testen |
| **`aux_check.py`** | **schwarz** | **Freilauf** | **Ch3, Ch4** | **Abnahme des Lineals** |
| **`ofdr_aux.py`** | **schwarz** | **Freilauf** | **Ch1–Ch4** | **die neue Auswertung** |

---

# Schritt für Schritt

## S0 · Basislinie · *alter Aufbau, nichts verändert*

**Aufbau:** wie er jetzt ist. **Gerät:** roter CoreDAQ, Triggermodus (1 pm),
2 Kanäle, feste Messbereiche.

```bash
python ofdr_process.py  scan_S0.json
python ofdr_diagnose.py tuning scan_S0.json
```

**Notiere:**

| | erwartet |
|---|---|
| Hauptpeak-Position | ~46,4 mm |
| Hauptpeak-Breite | ≤ ~21 µm |
| langsamer Bogen | ~59 pm p-p |
| schnelles Zittern | ~90–110 MHz rms |

**Wozu:** Das ist dein Vorher-Bild. Ohne diese vier Zahlen kannst du später
nichts vergleichen. Die Datei aufheben.

---

## S1 · Faserende abtöten · *nur die Testfaser angefasst*

**Aufbau:** Index-Matching-Gel auf die Endfläche, oder die letzten Zentimeter
eng aufwickeln, oder ein Patchcord < 35 cm. **Gerät:** unverändert.

```bash
python ofdr_process.py  scan_S1.json
python ofdr_diagnose.py compare scan_S0.json scan_S1.json
```

**Erwartung:** Der Median im Bereich 6–32 cm fällt deutlich. Das Dreckband
verschwindet ins Rauschen.

**Weiter, wenn:** Das Band ist kleiner geworden. Damit hast du die
Aliasing-Diagnose selbst bestätigt — und ab jetzt saubere Daten.

*(Falls es NICHT kleiner wird: stopp. Dann ist die Erklärung im Handover doch
nicht die ganze Wahrheit, und das musst du vor dem Umbau wissen.)*

---

## S2 · Parameter festlegen · *Schreibtisch, kein Gerät*

```bash
python ofdr_sweep_calc.py --table
python ofdr_sweep_calc.py --points 1000000 --span 60 --dl 4
```

Trag ein, was Giulio dir zur Pufferfrage gesagt hat. Wenn 1 Mio. **pro Kanal**:

| | |
|---|---|
| Punkte | 1.000.000 |
| Sweep-Geschwindigkeit | 60 nm/s |
| Sweep von / bis | 1505 → 1565 nm |
| Aux-Armdifferenz | 4 m |
| → Reichweite | 6,95 m |
| → Auflösung | 13,9 µm |
| → Aux erscheint bei | 2,00 m |
| → Punkte pro Aux-Fringe | 7,0 |

Wenn 1 Mio. **geteilt durch 4**, dann rechne mit `--points 250000` neu — dann
musst du auf 15 nm Span runter und landest bei 55,6 µm Auflösung.

**Ausdrucken und an den Aufbau hängen.**

---

## S3 · Üben an der Simulation · *Schreibtisch, kein Gerät*

Mach das **bevor** du Hardware anfasst, damit du die Ausgaben kennst.

```bash
python make_aux_test_data.py --points 400000 --speed 60 --dl 4 \
       --lam0 1505 --out sim.npz

python aux_check.py sim.npz --lam-start 1505 --lam-stop 1529

python ofdr_aux.py sim.npz --tau-aux-ns 19.587 --zmax 2.5
```

*(Die 1529 sind Start + 400.000 × 1 µs × 60 nm/s = 1505 + 24 nm.)*

**Erwartung:** `aux_check` meldet Kontrast gut, Monotonie gut, τ_aux ≈ 19,59 ns,
ΔL ≈ 4,00 m. `ofdr_aux` findet vier Reflektoren und schreibt am Ende
`-> BESTANDEN`.

Zum Anschauen, wie ein **kaputter** Fall aussieht:

```bash
python make_aux_test_data.py --points 400000 --ripple-mhz 100 --out sim_bad.npz
python aux_check.py sim_bad.npz --lam-start 1505 --lam-stop 1529
```

Dort wird die Monotonieprüfung Alarm schlagen. So sieht es aus, wenn der Laser
zu stark zittert.

---

## S4 · 90/10-Splitter einbauen · *eine einzige Änderung*

**Aufbau:** 90/10 direkt hinter den Laser. 90 % in den bestehenden Aufbau,
10 % erstmal offen lassen (oder abgeschlossen). **Gerät:** noch roter CoreDAQ,
noch Triggermodus.

```bash
python ofdr_process.py scan_S4.json
```

**Weiter, wenn:**

- Hauptpeak steht an derselben Stelle wie in S0 (± ein paar µm)
- Peakbreite unverändert
- Alles ~0,5 dB dunkler — das ist der eingefügte Verlust, normal

**Wozu:** Du prüfst, dass der Splitter nichts kaputt gemacht hat, **mit
Werkzeug, das du schon kennst**. Eine neue Sache nach der anderen.

---

## S5 · Aux bauen und abnehmen · *jetzt wird umgeschaltet*

**Aufbau:** 10 %-Zweig → 50/50-Splitter → zwei Arme, einer 4 m länger →
50/50-Koppler → beide Ausgänge auf **Ch3 und Ch4**.
Die 4-m-Spule in eine Schachtel legen, nicht frei über den Tisch.

**Gerät:** schwarzer CoreDAQ. Single-Trigger am EXFO, Freilauf mit 1 µs.
4 Kanäle, **feste Messbereiche** (Autoranging mitten im Sweep zerstört den
Scan). Sweep 1505 → 1565 nm bei 60 nm/s.

```bash
python aux_check.py scan_S5.json --lam-start 1505 --lam-stop 1565
```

**Das Programm beantwortet vier Fragen. Weiter nur, wenn alle vier passen:**

| Prüfung | Kriterium | wenn nicht |
|---|---|---|
| **1 Kontrast** | Sichtbarkeit überall > 0,4 | Polarisationsfading → Polarisationsregler in einen Arm |
| **2 Monotonie** | 0 % rückwärts laufende Schritte | **STOPP.** Laser zittert zu stark. Ursache suchen, bevor du weitermachst |
| **3 Kalibrierung** | τ_aux ≈ 19,6 ns, ΔL ≈ 4 m, ≥ 4 Punkte/Fringe | ΔL stimmt nicht mit der verlegten Faser überein → nachmessen |
| **4 Laserqualität** | Zittern deutlich unter der genannten Grenze | siehe Prüfung 2 |

Am Ende druckt das Programm die Zeile, die du brauchst:

```
--> DIESE ZAHL WEITERVERWENDEN:  ofdr_aux.py ... --tau-aux-ns 19.5871
```

**Schreib sie ins Laborbuch.** Sie ist dein Maßstab.

**Pufferkontrolle bei diesem Scan:**

```
Sweepdauer  = 60 nm / 60 nm/s     = 1,00 s
Pufferdauer = 1.000.000 × 1 µs    = 1,00 s     -> passt genau, keine Reserve
```

Bau lieber etwas Reserve ein: Sweep 55 nm statt 60, oder 65 nm/s statt 60.
Wenn das Flackern im letzten Teil der Datei abbricht, war der Sweep zu lang.

> **Warum Prüfung 4 wichtiger ist, als sie aussieht:** Das Handover nennt
> 90–110 MHz Zittern, gemessen *durch das Messsignal hindurch* und deshalb von
> Störlicht verfälscht. Bei ~88 MHz mit 4,5 pm Periode fährt der Laser
> zwischendurch rückwärts, und dann kann **kein** Verfahren mehr helfen. Das
> Aux misst diese Zahl zum ersten Mal sauber. Das ist die Annahme, auf der der
> ganze Umbau ruht — deshalb wird sie geprüft, bevor irgendjemand einem
> Reflektogramm glaubt.

---

## S6 · Erste echte Messung · *derselbe Scan, jetzt die Messkanäle*

Du brauchst keinen neuen Scan — `scan_S5.json` enthält bereits alle vier Kanäle.

```bash
python ofdr_aux.py scan_S5.json --tau-aux-ns 19.5871 --zmax 3.0
```

**Drei Kontrollen, in dieser Reihenfolge:**

**Kontrolle 1 — steht das Aux, wo es soll?**
Im Reflektogramm muss ein Peak bei **z = ΔL/2 = 2,00 m** sein. Das ist das Aux
selbst, das über den 20-dB-Koppler mit einstreut. Wenn er woanders sitzt,
stimmt τ_aux nicht.

**Kontrolle 2 — der beste Test, den du hast:**
Der alte, bekannte Peak bei **46,4 mm** muss an **genau derselben Stelle**
stehen wie in S0 und S4. Du kennst die Antwort im Voraus, gemessen mit einer
völlig anderen Methode. Wenn beide übereinstimmen, funktioniert die neue Kette
end-to-end.

Toleranz: eine Auflösungszelle, also ~14 µm.

**Kontrolle 3 — der Moment der Wahrheit:**
Das Faserende bei ~1 m muss jetzt als **scharfer Peak an seiner echten
Position** stehen — und nicht mehr als Dreckband bei 21 cm.

**Und die Dauerprüfung:** Peakbreite ≲ 1,5 × Fenstergrenze (das Programm
druckt beide Zahlen nebeneinander).

---

## Danach

Wenn S6 sitzt, kannst du die Konfiguration hochdrehen:

```bash
python ofdr_sweep_calc.py --points 1000000 --speed 120 --dl 4
```

120 nm/s über 120 nm Span → Reichweite 3,48 m bei **7,0 µm** Auflösung, also
die volle Schärfe von heute mit 8-facher Reichweite. Der Haken: das Aux hat
dann nur noch 3,5 Punkte pro Fringe. Das ist über Nyquist, aber ohne Reserve —
`aux_check.py` wird es dir sagen. Deshalb erst danach, nicht davor.

---

## Wenn etwas nicht passt

| Symptom | Erster Verdacht | Test |
|---|---|---|
| Kontrast bricht stellenweise ein | Polarisation in den Aux-Armen | `aux_check.py`, Panel oben rechts |
| Monotonie-Alarm | Laser zittert zu stark | `aux_check.py`, Prüfung 2 — **nicht** weitermachen |
| τ_aux passt nicht zur verlegten Faser | Start/Endwellenlänge falsch eingegeben, oder Sweep nicht komplett im Puffer | Pufferrechnung nachziehen |
| 46,4-mm-Peak an falscher Stelle | τ_aux falsch, oder Kanäle vertauscht | Ch1–Ch4 Zuordnung prüfen |
| Peaks doppelt so weit weg wie erwartet | Reflexion/Transmission verwechselt | siehe die „÷2"-Erklärung |
| Flackern bricht am Dateiende ab | Sweep länger als der Puffer | Span kleiner oder schneller sweepen |
| Alles sieht falsch aus, keine Idee | erst messen, dann erklären | `ofdr_diagnose.py` hat vier fertige Tests |
