"""
Double Pulse Test waveform controller for BK Precision 4054B AFG.
Uses Siglent-style C1: command syntax (4054B is Siglent OEM).

Waveform structure:
    ___________         _______
   |     T1    |  T2   |  T3   |
___/           \_______/       \___...

SAFETY:
  - Output is forced OFF immediately on connect, and again after *RST.
  - BTWV IDLE,BOTTOM holds output at 0 V before firing and permanently after.
  - After fire(), the waveform ends on the tail of zeros; the instrument
    holds 0 V indefinitely — even if the USB cable is disconnected.
  - output_on() enables the output connector but NO signal appears until
    fire() is explicitly called.

Usage:
    afg = DPT_AFG()
    afg.setup_dpt(T1_us=100, T2_us=30, T3_us=30)
    afg.output_on()   # output connector live but held at 0 V
    input("Press Enter to fire...")
    afg.fire()        # one-shot DPT sequence, then returns to 0 V
    afg.close()       # turns output off and releases USB
"""

import pyvisa
import numpy as np
import struct
import time


# ---------------------------------------------------------------------------
# Waveform construction
# ---------------------------------------------------------------------------

def build_dpt_waveform(T1_us: float, T2_us: float, T3_us: float,
                       sample_rate: float = 1e6,
                       lead_us: float = 0.0,
                       tail_us: float = 10.0) -> tuple[np.ndarray, float, int]:
    """
    Build a DPT waveform array of signed 16-bit DAC values.

    Args:
        T1_us:       First pulse width  (µs)
        T2_us:       Inter-pulse gap    (µs)
        T3_us:       Second pulse width (µs)
        sample_rate: AWG sample rate    (Sa/s), default 1 MSa/s → 1 µs resolution
        lead_us:     Leading LOW time   (µs) before T1 — useful for scope pre-trigger validation
        tail_us:     Trailing LOW time  (µs) appended after T3

    Returns:
        (binary bytes, total period in seconds, number of points)
    """
    def us_to_samples(t_us: float) -> int:
        n = int(round(t_us * 1e-6 * sample_rate))
        if n < 1:
            raise ValueError(
                f"Timing {t_us} µs produces < 1 sample at {sample_rate/1e6:.1f} MSa/s. "
                "Increase sample_rate or timing values."
            )
        return n

    nlead = us_to_samples(lead_us) if lead_us > 0 else 0
    n1    = us_to_samples(T1_us)
    n2    = us_to_samples(T2_us)
    n3    = us_to_samples(T3_us)
    ntail = us_to_samples(tail_us)
    total = nlead + n1 + n2 + n3 + ntail

    if total < 8:
        raise ValueError(f"Total waveform is {total} points; instrument requires >= 8.")

    # Signed 16-bit little-endian: -32768 = LLEV = 0 V (LOW), +32767 = HLEV = 5 V (HIGH)
    samples = np.full(total, -32768, dtype=np.int16)
    # nlead samples remain -32768              # lead LOW (0V) before T1
    samples[nlead        : nlead+n1]         = 32767   # T1 HIGH (5V)
    # n2 samples remain -32768                # T2 LOW  (0V)
    samples[nlead+n1+n2  : nlead+n1+n2+n3]  = 32767   # T3 HIGH (5V)
    # ntail samples remain -32768             # tail LOW (0V) — instrument idles here after burst

    # Pack as raw little-endian bytes for write_raw transfer
    binary = struct.pack(f'<{total}h', *samples.tolist())
    period = total / sample_rate
    return binary, period, total


# ---------------------------------------------------------------------------
# Instrument interface
# ---------------------------------------------------------------------------

class DPT_AFG:
    """
    Controls BK Precision 4054B for single-shot Double Pulse Test sequences.
    Uses Siglent C1:/C2: command syntax.

    Output swings 0–5 V (High-Z), suitable for gate drivers.
    Output is NEVER non-zero until fire() is explicitly called.
    """

    RESOURCE = 'USB0::62700::60984::515I17148::0::INSTR'

    def __init__(self, resource_name: str | None = None, channel: int = 1):
        self.rm = pyvisa.ResourceManager('@py')
        self.channel = channel
        self._ch = f"C{channel}"

        if resource_name is None:
            resource_name = self.RESOURCE

        print(f"Opening: {resource_name}")
        self.inst = self.rm.open_resource(resource_name)
        self.inst.timeout           = 15_000
        self.inst.write_termination = '\n'
        self.inst.read_termination  = '\n'

        # ---- SAFETY: force output off before anything else ----
        self._write(f"{self._ch}:OUTP OFF")

        idn = self.inst.query("*IDN?").strip()
        print(f"Connected: {idn}")
        print(f"Output forcibly OFF.")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def setup_dpt(self,
                  T1_us: float,
                  T2_us: float,
                  T3_us: float,
                  sample_rate: float = 1e6,
                  lead_us: float = 10.0):
        """
        Upload and configure the DPT waveform. Output remains OFF.

        Args:
            T1_us:       First pulse width  (µs)
            T2_us:       Inter-pulse gap    (µs)
            T3_us:       Second pulse width (µs)
            sample_rate: AWG sample rate in Sa/s.
                         1e6  → 1 µs resolution  (default)
                         10e6 → 100 ns resolution (sub-µs pulses)
            lead_us:     Leading LOW time before T1 (µs), default 0.
                         Set e.g. 20 to add a 20µs dead zone before T1 for scope validation.
        """
        binary, period, n_pts = build_dpt_waveform(T1_us, T2_us, T3_us, sample_rate, lead_us)
        freq = 1.0 / period

        print(f"\nDPT waveform:")
        print(f"  T1 = {T1_us} µs,  T2 = {T2_us} µs,  T3 = {T3_us} µs")
        print(f"  {n_pts} points @ {sample_rate/1e6:.1f} MSa/s  →  period = {period*1e6:.1f} µs")
        print(f"  Playback frequency = {freq:.4f} Hz")

        # ---- SAFETY: output off before reset ----
        self._write(f"{self._ch}:OUTP OFF")
        # Force known polarity so LOW level corresponds to 0 V idle.
        self._write(f"{self._ch}:OUTP PLRT,NOR")

        self._write("*RST")
        time.sleep(2.0)

        # ---- SAFETY: output off again — RST restores default state ----
        self._write(f"{self._ch}:OUTP OFF")
        # Re-apply after reset in case the instrument recalls an inverted state.
        self._write(f"{self._ch}:OUTP PLRT,NOR")
        time.sleep(0.2)

        # ---- Upload waveform via binary write_raw ----
        # TYPE,8 = arbitrary, WAVEDATA followed by raw little-endian int16 bytes.
        # Must use write_raw (not write) — the binary data contains bytes that
        # would be corrupted by ASCII encoding or the write terminator logic.
        header = (
            f"{self._ch}:WVDT WVNM,DPT,"
            f"FREQ,{freq:.6g},"
            f"TYPE,8,"
            f"AMPL,5.0,"
            f"OFST,2.5,"
            f"PHASE,0,"
            f"MARK,OFF,"
            f"WAVEDATA,"
        ).encode('ascii')
        self.inst.write_raw(header + binary + b'\n')
        time.sleep(1.5)

        # ---- Switch to ARB mode, then select the uploaded waveform ----
        # BSWV WVTP,ARB must come first — it puts the channel in ARB mode.
        # ARWV NAME,DPT then selects our specific uploaded waveform.
        # Order matters: ARB mode first, then waveform selection.
        self._write(f"{self._ch}:BSWV WVTP,ARB")
        time.sleep(0.5)
        self._write(f"{self._ch}:ARWV NAME,DPT")
        time.sleep(0.5)

        # ---- Waveform voltage levels and frequency ----
        # Use HLEV/LLEV (absolute high/low voltage) instead of AMP+OFST.
        # This means the idle bottom is unambiguously 0 V — no offset to leak through.
        # DAC 0 → LLEV = 0 V, DAC 16383 → HLEV = 5 V.
        self._write(f"{self._ch}:BSWV HLEV,5.0")
        self._write(f"{self._ch}:BSWV LLEV,0.0")
        self._write(f"{self._ch}:BSWV FREQ,{freq:.6g}")
        self._write(f"{self._ch}:OUTP PLRT,NOR")

        # Diagnostic: confirm what frequency and sample rate the instrument reports
        bswv = self.inst.query(f"{self._ch}:BSWV?").strip()
        outp = self.inst.query(f"{self._ch}:OUTP?").strip()
        print(f"  Instrument reports: {bswv}")
        print(f"  Output config: {outp}")
        print(f"  Expected: FREQ={freq:.4f} Hz  →  {1/freq*1e6:.1f} µs per waveform cycle")
        print(f"  Expected timing: T1={T1_us}µs  T2={T2_us}µs  T3={T3_us}µs  gap between triggers = {T1_us+T2_us:.0f}µs  (each sample = {1e6/sample_rate:.2f} µs)")

        # ---- Load: High-Z ----
        self._write(f"{self._ch}:OUTP LOAD,HZ")

        # ---- Burst: 1 cycle, manual trigger, IDLE at bottom (0 V) ----
        # IDLE,BOTTOM: output is held at 0 V before AND after the burst fires.
        # This is hardware-latched — persists even if USB is disconnected.
        self._write(f"{self._ch}:BTWV STATE,ON")
        self._write(f"{self._ch}:BTWV GATE_NCYC,NCYC")
        self._write(f"{self._ch}:BTWV TIME,1")
        self._write(f"{self._ch}:BTWV TRSR,MAN")
        self._write(f"{self._ch}:BTWV DLAY,0")
        self._write(f"{self._ch}:BTWV IDLE,BOTTOM")

        # ---- SAFETY: confirm output still OFF ----
        self._write(f"{self._ch}:OUTP OFF")

        print("Setup complete.")
        print("  → Call output_on() to arm (output held at 0 V until fire()).")

    def output_on(self):
        """
        Enable the output connector. Signal stays at 0 V until fire() is called.
        The burst idle state (BOTTOM) ensures no pulse is emitted.
        """
        self._write(f"{self._ch}:OUTP ON")
        print(f"Channel {self.channel} output ON — held at 0 V. Call fire() when ready.")

    def output_off(self):
        """Disable the output connector (0 V, high impedance at the BNC)."""
        # Zero the offset before disabling — prevents a 2.5V glitch as the
        # instrument transitions from burst-idle to output-off state.
        self._write(f"{self._ch}:BSWV OFST,0.0")
        self._write(f"{self._ch}:OUTP OFF")
        print(f"Channel {self.channel} output OFF.")

    def fire(self):
        """
        Trigger one DPT pulse sequence: T1 HIGH → T2 LOW → T3 HIGH → 0 V forever.
        After the burst the instrument holds 0 V indefinitely (hardware state,
        independent of USB connection).
        """
        self._write(f"{self._ch}:BTWV MTRIG")
        print("DPT sequence fired. Output will hold 0 V after the burst completes.")

    def close(self):
        """Disable output and release the VISA resource."""
        try:
            self.output_off()
        finally:
            self.inst.close()
            self.rm.close()
            print("Connection closed.")

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------

    def query(self, cmd: str) -> str:
        """Send a query and print the response."""
        resp = self.inst.query(cmd).strip()
        print(f"  {cmd!r}  →  {resp!r}")
        return resp

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _write(self, cmd: str):
        self.inst.write(cmd)
        time.sleep(0.05)


# ---------------------------------------------------------------------------
# Main — edit T1/T2/T3 here and run directly
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # ---- Set your DPT timing here ----
    T1 = 55.0   # µs — first pulse  (sets inductor current level)
    T2 =  5.0   # µs — off time     (freewheeling / dead time)
    T3 =  5.0   # µs — second pulse (switching event under load current)

    # Sample rate determines rise/fall time (1 sample = 1 transition):
    #   1e6  Sa/s →  1 µs per sample (rise time ~1 µs)   ← previous setting
    #   25e6 Sa/s → 40 ns per sample (rise time ~9 ns, analog BW limited)
    # Memory constraint: 4054B holds ~16k points max per waveform.
    # At 25 MSa/s: 16384 points = 655 µs max waveform — sufficient for most DPT.
    SAMPLE_RATE = 25e6

    afg = DPT_AFG()
    afg.setup_dpt(T1, T2, T3, sample_rate=SAMPLE_RATE)
    afg.output_on()   # held at 0 V — safe to connect gate drivers now

    input("\nPress Enter to fire the DPT sequence (or Ctrl+C to abort)...")
    afg.fire()

    input("\nPress Enter to disable output and exit...")
    afg.close()
