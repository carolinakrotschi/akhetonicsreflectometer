# Instrument driver: Thorlabs TOF1550

`TunableFilters.py` — driver for the C-band digital tunable bandpass filter
used as the wavelength marker for the λ(t) calibration.

It lives here rather than in the repository's `Interface/` folder so that
`Interface/` stays exactly as it was. **If the filter should become a regular
GUI instrument**, this file belongs in `Interface/`, its configuration in
`Interface/config/`, and the device in the GUI's device registry.

Reasoning and measurements: `../report.md`, §2.

## Specification (User Guide TTN356097-D02)

| Property | Value |
|---|---|
| Tuning range | 1527 – 1567 nm |
| Passband (−3 dB) | 0.21 nm — measured here 0.15–0.23 nm |
| Setting error (absolute) | < 32 pm — measured here ~15 pm at 1550 nm |
| Repeatability | ±8 pm |
| Interface | USB-C, FTDI VCP, 115200 8N1 |

## Two deviations from the manual, found on the unit

1. Replies are **LF**-terminated, not CR as documented.
2. Replies are prose with a prefix: `OK. The center wavelength is 1550.000 nm.`
   The driver therefore extracts the first number from the reply.

## Port

Pinned in `config/TunableFilter_config.json` to **COM6**. It has to be pinned:
the unit enumerates through the generic FTDI VCP driver (VID 0403, "USB Serial
Port"), so it is indistinguishable from any other FTDI adapter, and it has no
identity command. Note that COM17 on this machine is a *Thorlabs TEC
controller* — the Thorlabs USB VID is not a reliable hint.

## Use

```python
from lina.interface.TunableFilters import TOF1550
f = TOF1550('tof1550-1')
f.open()
f.set_wavelength(1550.0)      # nm, returns the read-back value
print(f.get_wavelength())
f.close()
```

Set/read-back was exact to 0 pm at 1545, 1550.5 and 1560 nm. There is no
motion-status query on this instrument, so `set_wavelength` waits a fixed
`settle_s` (1 s by default, from the config) and then verifies by read-back.
