# Befund und Korrekturen

## Messung: aliased

| Größe | Wert |
|---|---|
| Schrittweite δλ | 0,9998 pm |
| Bandbreite | 49,908 nm = 6,268 THz |
| Auflösung δz = c/(2 n_g B) | 16,3 µm |
| **Nyquist-Reichweite** λ²/(4 n_g δλ) | **40,66 cm Reflektorabstand** (81,3 cm OPD) |
| Peak des Reflektogramms | 34,8 cm — kein Peak, ein breiter Berg |
| Energie in den obersten 20 % vor Nyquist | **68 %** |

Ch1 und Ch2 korrelieren mit **−0,93** (gegenphasige Kopplerports), die Summe
Ch1+Ch2 ist nahezu konstant (σ/µ = 0,10), die Differenz trägt den vollen Fringe
(σ/µ = 0,52). Das Signal ist also **echte Interferenz mit voller Modulation** —
kein Rauschen. Es ist nur jenseits von Nyquist.

Beste Einzelziel-Hypothese aus drei unabhängigen Verfahren
(erste Faltung, Chirp-Steigung über Teilbänder, dechirpter NUDFT-Scan):
**OPD ≈ 92 cm, entspricht 46 cm Reflektorabstand.** Das ist ~13 % jenseits der
Grenze. Belastbar wird die Zahl erst mit einem Sweep bei δλ ≤ 0,5 pm.

## Bugs in fft_reflectometer.py

1. **Letztes Sample ist Müll.** `Wavelength[-1] = 1519.95` (Rücksprung auf den
   Sweep-Anfang). Dadurch wird `np.mean(np.diff(wavelength))` zu **2,40e-6 nm**
   statt 1e-3 nm — die Frequenzachse ist um Faktor **416** falsch skaliert.
   → nicht-monotone Punkte verwerfen, `np.median` statt `np.mean`.

2. **Falsche Umrechnung im λ-Zweig (Default).** `length_nm = peak_frequency /
   refractive_index` gilt nur im Wellenzahl-Zweig. In der λ-Domäne ist
   f = n_g·ΔL/λ², also **ΔL = f·λ²/n_g**. Es fehlt λ² ≈ 2,39e6 nm².

3. **Faktor 2 fehlt.** Bei Reflexionsmessung ist der Reflektorabstand OPD/2.

4. **Keine DC-/Hüllkurven-Entfernung.** `argmax` über `freq > 0` greift Bin 1–2
   (Laserleistungshülle) statt des Signals. Ergebnis der alten Version:
   5,7e-9 m bzw. 1,3e-2 m — beides Artefakte.

5. **Kein Fenster.** Ohne Hann verschmieren Nebenkeulen schwache Reflexe.

6. **Ch2 ungenutzt.** `(Ch1−Ch2)/(Ch1+Ch2)` entfernt die Leistungshülle exakt
   und verdoppelt den Fringe-Kontrast. Geschenkt.

7. **Keine Nyquist-Prüfung.** Das Skript hat nie gewarnt, dass die Messung
   außerhalb des gültigen Bereichs liegt.

8. `plt.xlim(0, 3)` suggeriert 3 m Messbereich. Real sind es 40,7 cm.

9. Die gemeldeten Wellenlängen sind auf 1 pm **gerundet** (Schrittfolge
   0/1/2 pm). Für die Phase ist ein ideal gleichmäßiges Raster
   (`linspace(λ_0, λ_N, N)`) die bessere Annahme — die gerundeten Werte
   addieren bis zu ±0,5 pm Jitter, was bei OPD ≈ 1 m schon ~1,9 rad
   Phasenfehler bedeutet.

## Nächster Schritt

δλ = 0,5 pm → Reichweite 81 cm, 100 000 Punkte.
δλ = 0,2 pm → Reichweite 2,03 m, 250 000 Punkte (Sicherheitsfaktor 2 für 1 m).

Auflösung bleibt in allen Fällen 16,3 µm — die hängt nur an der Bandbreite.
