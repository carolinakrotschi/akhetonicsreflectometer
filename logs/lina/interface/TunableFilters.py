"""
Tunable optical filters.

TOF1550 — Thorlabs C-band digital tunable Gaussian bandpass filter.

Protocol (User Guide TTN356097-D02, §4.2): FTDI virtual COM port,
115200 8N1, no flow control. Commands are SCPI-like but abbreviated, one per
write, terminated with <NL>; replies are terminated with <CR> and start with
"OK" on success or "Err: ..." on failure. Multiple commands in one string are
NOT accepted, and exponential notation in parameters is NOT accepted.

Command set is small: WAVElength? / WAVElength X (nm), WAVElength:STep? /
WAVElength:STep X (pm), WAVElength:INcrease / :DEcrease, CHANnel? / CHANnel X
(ITU). There is no motion/settled status query, so set_wavelength() waits a
fixed settle time (``settle_s``, from config) and verifies by readback.
"""

import json
import re
import time
from pathlib import Path
from typing import Optional

import serial
import serial.tools.list_ports

# This driver lives under lina/ rather than Interface/, so it imports the base
# class from the repo's Interface package explicitly (the relative-import
# fallback used by the drivers inside Interface/ cannot work from here).
try:
    from Interface.LabDevice import LabDevice
except ImportError:                      # running with Interface/ on sys.path
    from LabDevice import LabDevice


_TF_CFG_CANDIDATES = [
    Path(__file__).resolve().parent / 'config' / 'TunableFilter_config.json',
    Path(__file__).resolve().parent / 'TunableFilter_config.json',
]


def _resolve_tf_config() -> Path:
    for p in _TF_CFG_CANDIDATES:
        if p.exists():
            return p
    raise FileNotFoundError(
        "TunableFilter_config.json not found. Tried:\n" +
        "\n".join(str(p) for p in _TF_CFG_CANDIDATES)
    )


class TOF1550(LabDevice):
    """Thorlabs TOF1550 digital tunable C-band bandpass filter.

    Specified 1527–1567 nm tuning range, 0.21 nm passband at -3 dB, <32 pm
    absolute setting error, ±8 pm repeatability — narrow and accurate enough
    to be used as a wavelength marker in a laser sweep.
    """

    WL_MIN_NM = 1527.0
    WL_MAX_NM = 1567.0

    def __init__(self, name: str = "tof1550-1", config_path: Optional[str] = None,
                 port: Optional[str] = None):
        self._baudrate = 115200
        self._timeout = 1.0
        self._settle_s = 1.0
        cfg_port = None
        try:
            cfg_path = Path(config_path).expanduser().resolve() \
                if config_path else _resolve_tf_config()
            with open(cfg_path) as fp:
                cfg = json.load(fp)
            dev_cfg = cfg.get('TOF1550', {}).get(name, {})
            cfg_port = (dev_cfg.get('serial_id') or '').strip() or None
            self._baudrate = int(dev_cfg.get('baudrate', self._baudrate))
            self._timeout = float(dev_cfg.get('timeout_s', self._timeout))
            self._settle_s = float(dev_cfg.get('settle_s', self._settle_s))
        except Exception:
            pass   # no config — port must be passed explicitly

        self.port = port or cfg_port
        if self.port and self.port.upper() in ('AUTO', 'NA', 'NONE'):
            self.port = None
        self.ser: Optional[serial.Serial] = None
        super().__init__(name, serial_id=self.port or 'AUTO')

    # ---------- discovery ----------

    @staticmethod
    def find_ports() -> list:
        """Candidate COM ports: the unit enumerates through the generic FTDI VCP
        driver (VID 0x0403, descriptor "USB Serial Port"), so it is NOT
        distinguishable from any other FTDI adapter by descriptor alone, and it
        has no identity command. Pin the port in TunableFilter_config.json;
        this list is only a fallback shortlist."""
        return [p.device for p in serial.tools.list_ports.comports()
                if p.vid in (0x0403, 0x1313)]

    # ---------- connection management ----------

    def open(self):
        if self.port is None:
            candidates = self.find_ports()
            if not candidates:
                raise ConnectionError(
                    "No Thorlabs FTDI COM port found for the TOF1550. Set "
                    "'serial_id' in TunableFilter_config.json to its COM port.")
            self.port = candidates[0]
            self.logger.info(f"Auto-selected {self.port}")
        try:
            self.ser = serial.Serial(
                port=self.port, baudrate=self._baudrate,
                bytesize=serial.EIGHTBITS, parity=serial.PARITY_NONE,
                stopbits=serial.STOPBITS_ONE, timeout=self._timeout,
                write_timeout=self._timeout,
                rtscts=False, dsrdtr=False, xonxoff=False)
        except serial.SerialException as ex:
            self.logger.error(f"Could not open {self.port}: {ex}")
            raise ConnectionError(f"Could not open {self.port}") from ex
        try:
            self.ser.reset_input_buffer()
            self.ser.reset_output_buffer()
        except Exception:
            pass
        wl = self.get_wavelength()
        self.logger.info(f"Connected on {self.port} @ {self._baudrate} bps, "
                         f"center {wl:.3f} nm")
        return True

    def close(self):
        if self.ser is not None:
            try:
                self.ser.close()
            except Exception as ex:
                self.logger.warning(f"Error closing port: {ex}")
            self.ser = None
            self.logger.info("Disconnected")

    # ---------- I/O primitives ----------

    def _send(self, command: str) -> str:
        """Send one <NL>-terminated command, return the reply with the leading
        "OK" stripped. Raises RuntimeError on an "Err: ..." reply."""
        if self.ser is None:
            raise RuntimeError(f"{self.id}: not connected")
        self.ser.reset_input_buffer()
        self.ser.write((command.strip() + "\n").encode("ascii"))
        self.ser.flush()
        # The manual specifies <CR>-terminated replies; the unit actually sends
        # <LF> ("OK. The center wavelength is 1550.000 nm.\n").
        raw = self.ser.read_until(b"\n").decode("ascii", "replace").strip()
        if not raw:
            raise RuntimeError(f"{self.id}: no reply to {command!r}")
        if raw.lower().startswith("err"):
            raise RuntimeError(f"{self.id}: {raw} (command {command!r})")
        return re.sub(r'^\s*OK[.:]?\s*', '', raw).strip()

    @staticmethod
    def _first_float(reply: str) -> float:
        m = re.search(r'[-+]?\d+(?:\.\d+)?', reply)
        if not m:
            raise ValueError(f"No number in reply {reply!r}")
        return float(m.group(0))

    # ---------- wavelength ----------

    def get_wavelength(self) -> float:
        """Current filter center wavelength in nm."""
        return self._first_float(self._send("WAVElength?"))

    def set_wavelength(self, wavelength_nm: float, wait: bool = True) -> float:
        """Set the center wavelength (nm) and return the value read back.

        There is no settled/busy query on this instrument, so ``wait`` just
        sleeps the configured ``settle_s`` before the readback.
        """
        wl = float(wavelength_nm)
        if not (self.WL_MIN_NM <= wl <= self.WL_MAX_NM):
            raise ValueError(
                f"{wl} nm is outside the TOF1550 range "
                f"({self.WL_MIN_NM}–{self.WL_MAX_NM} nm)")
        self._send(f"WAVElength {wl:.3f}")
        if wait:
            time.sleep(self._settle_s)
        return self.get_wavelength()

    def get_step_pm(self) -> float:
        """Increment/decrement step in pm."""
        return self._first_float(self._send("WAVElength:STep?"))

    def set_step_pm(self, step_pm: float):
        self._send(f"WAVElength:STep {float(step_pm):.0f}")

    def increase(self):
        """Increase the center wavelength by the configured step."""
        self._send("WAVElength:INcrease")

    def decrease(self):
        """Decrease the center wavelength by the configured step."""
        self._send("WAVElength:DEcrease")

    # ---------- ITU channel ----------

    def get_channel(self) -> int:
        """Closest ITU channel — NOTE: per the manual this query also SNAPS the
        center wavelength to that channel, so it is not a passive read."""
        return int(self._first_float(self._send("CHANnel?")))

    def set_channel(self, channel: int):
        self._send(f"CHANnel {int(channel)}")

    # ---------- state snapshot ----------

    def get_settings(self) -> dict:
        s = super().get_settings()
        s.update({
            "port": self.port,
            "wavelength_nm": self._safe(self.get_wavelength),
            "step_pm": self._safe(self.get_step_pm),
        })
        return s
