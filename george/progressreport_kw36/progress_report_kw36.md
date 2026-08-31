

Hi — quick update on what happened since the last report. 

---



## 2. No-fiber connector peak moved — but it's not a new mystery

first wanted to see if setup was exactlz like i left it so i compared to the alst one of the 20.8 , because on 20.8 i changed the connector in the last scan, i now have 2 internallz consistent families of the first peak (which i still think is from the first connector - but then whz dont i have peaks for the other connectors?) (587 mm-ish and 496 mm-ish),


## 3. Aux calibration reproducibility — mostly fine, one scan clearly bad

Ran the same no-fiber setup three times in a row to see how stable the
τ_aux calibration is scan-to-scan:

| scan | τ_aux (ns) | dL (m) | monotonicity | slow bow (pm p-p) |
|---|---|---|---|---|
| 1st | 17.6028 | 3.5948 | fails, 0.02% backward steps | 44.0 |
| 2nd | 17.6058 | 3.5954 | passes | 43.6 |
| 3rd | 17.0621 | 3.4844 | fails, 0.93% backward steps | 1421.5 |

this was the big question of the day - why is the aux calibration so bad, it should be 4m because the aux interferometer, I have a path difference of 4m which I calculated this waz... t aux is calculated this way ...

---

## 4. New fiber measurements + a length discrepancy worth flagging

Measured the ~1 m fiber and the "3 m" fiber (now labeled "3mB" to
distinguish it from the physically different 3 m fiber used last time,
which is now "3mA") on two different aux configurations (aux built to
~5 m and ~4 m dL). Same fiber, different aux config, same secondary peak
position to ~1 mm — good, that's exactly what should happen since the aux
only sets the axis calibration. Also measured 3mB + a 1 m fiber ("1mA")
in series: predicted 3527.5 mm from adding the two individual
measurements, got 3523.8 mm, 4 mm off — fine.

**But:** every fiber measured today reads *short* of its round-number
label by 9–13% (the "~1 m" fiber → 912 mm, "3 m" → 2616 mm). That's the
opposite direction from the fibers in the last report, which read *long*
by 6–9%. Dug into this a bit:

- Forcing the axis to the *nominal* aux length (4/5/6 m) instead of the
  fringe-counting-calibrated one makes the fiber lengths land within ~2%
  of their round-number labels — much better.
- But it makes the cross-scan connector-position consistency noticeably
  *worse* (spread goes from ~1 mm to ~10 mm).

So it looks like the fringe-counting calibration in the aux script has a
small, roughly constant ~10–15% multiplicative bias, and the nominal aux
length happens to approximately cancel it because the aux was built close
to its target — not that the aux fiber is physically ~50 cm short.  **Root cause still open.**

---

## 5. Problems / mix-ups sorted out today

- Repeat scan of the 3mB+1mA combo confirmed the aux4m build is
  reproducible (dL agreed to 0.4 mm between two scans) — not a one-off.

---

## Still open

- The ~10–15% calibration bias in the aux fringe-counting (see point 4) —
  haven't root-caused it yet. I dont understand this since when I measured the fibers the last time, the calibrated aux lenght said something around 4m - include measuremt aux lenght table here
- Still haven't gotten rid of the fixed connector reflection itself
  (587 mm / 496 mm family, whichever connector's in the slot right now).
- τ_aux calibration still needs to get tighter before I'd trust an
  absolute fiber length number in a report.

Plots attached: connector-family comparison (today vs. last session),
fiber lengths relative to each scan's own connector peak, and the
nominal-vs-calibrated aux correction comparison.

— Carolina
