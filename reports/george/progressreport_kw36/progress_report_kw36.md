Hi — quick update on what happened today

First, I wanted to check that the setup was exactly as I had left it, so I compared today’s scan with the last one from 20 August. That was the session where I changed the connector. Last time changing the connector gave me two families for the first peak (depending on which connector I use): one at around 587 mm and one at around 496 mm, today I was still at the 496mm position, which told me nothing had changed

Then I wondered:

The aux-arm length has to be known very precisely for the analysis to work at all. I therefore wanted to see what actually happens when it is known imprecisely. That is why I ran several repeat measurements and compared the calculated aux length.

τ_aux is the group delay between the two aux-arm paths and from it the aux lenght can be calculated.

τ_aux = Δφ_aux,unwrapped / (2π · Δf_sweep)

| Scan | τ_aux (ns) | dL (m) |
| 1st | 17.6028 | 3.5948 |
| 2nd | 17.6058 | 3.5954 |
| 3rd | 17.0621 | 3.4844 |

This was were the big question of the day showed up: no matter which fibre I used or in which order the fibres were connected, the reported aux length always came out around 0.5 m shorter than I would expect. It should read 4 m, since that is the path difference the aux interferometer was built to have, but its always around 3.5m . The measurements on the 20.8.26 gave me a value closer to 4m, but I couldnt bring it back??? It cant be the connector change, because the connector only affects the measurement MZI not the Reference MZI

I dont think, that the problem comes from the software.   . When I manually forced the aux length to the nominal 4 m—which I already know is what was built —the results all landed within approximately 2% of their round-number labels, which is much better.

| Fibre | Before aux correction | With aux correction |
| 1m #1  | 1062.4 mm | 1002.0 mm |20th august
| 1m #2  | 1079.5 mm | 1017.9 mm |20th august
| 1m #3  | 1067.1 mm | 1001.9 mm |20th august
| 1m #4  | 1086.0 mm | 1019.0 mm |20th august
| 1m aux5m  | 912.0 mm | 1015.7 mm | 31st August
| 3m | 3147.2 mm | 2915.4 mm |20th august
| 3m aux5m  | 2616.5 mm | 2913.4 mm | 31st August
| 3m aux4m  | 2615.5 mm | 2923.9 mm | 31st August
| 3m + 1m aux4m | 3523.8 mm | 3939.3 mm | 31st August

However, I still wouldn’t trust that result.

I think there is an error in calibration and I dont think that all the results are made up or something, since...

- I measured 3m and a 1 m fibre  individually and then connected them in series. From the two individual measurements, I predicted a total length of 3527.5 mm. The combined measurement gave 3523.8 mm, so the difference was only around 4 mm, which seems fine.

- I measured the approximately 1 m fibre and the 3 m fibre
using two different aux configurations, with the aux built to approximately 5 m and 4 m dL. In both cases, the same secondary peak appeared at the same position within approximately 1 mm, which is good.


I have now tested all the fibre and connector configurations I could think of, both for the aux length and the measurement length. At this point, I’m not really sure how to continue because I still don’t know how to fix the problem, that the aux fiber is calibrated rightly...

Plots attached:

- **Fibre comparison, aux corrected:** Here, I forced the aux length to 4 m instead of calculating it using the τ_aux formula above. All the peaks are then at approximately the positions I expect. I also forced the first large peak to be at zero because, as explained in my last progress report, my current theory is that its position has to be subtracted.
- **Fibre comparison relative to the connector:** Here, I again forced the first peak to be at zero, but this time I did not force the aux length to 4 m.

In general, the plots without a date in their labels are from 20 August, while the plots labelled “31.8” are from today.

— Carolina