"""
Test binary little-endian WVDT upload with TYPE,8 on BK Precision 4054B.
"""
import pyvisa
import struct
import time

rm = pyvisa.ResourceManager('@py')
inst = rm.open_resource('USB0::62700::60984::515I17148::0::INSTR')
inst.timeout = 15_000
inst.write_termination = '\n'
inst.read_termination  = '\n'

def q(cmd):
    try:
        resp = inst.query(cmd).strip()
    except Exception as e:
        resp = f"ERROR: {e}"
    print(f"  -> {resp}")
    return resp

def w(cmd, desc=""):
    print(f"  W: {desc or cmd[:120]}")
    inst.write(cmd)
    time.sleep(0.5)

def reset():
    w("C1:OUTP OFF")
    w("*RST"); time.sleep(2.0)
    w("C1:OUTP OFF")

# 8-point 50% square wave: 4x HIGH (+32767), 4x LOW (-32768)
# Little-endian signed 16-bit
samples = [32767, 32767, 32767, 32767, -32768, -32768, -32768, -32768]
binary_data = struct.pack(f'<{len(samples)}h', *samples)

# Build the command: header as ASCII, data as raw bytes, terminated with \n
header = "C1:WVDT WVNM,SQTEST,FREQ,1000,TYPE,8,AMPL,5.0,OFST,2.5,PHASE,0,MARK,OFF,WAVEDATA,".encode('ascii')
cmd_bytes = header + binary_data + b'\n'

print("=== 1. Reset ===")
reset()

print(f"\n=== 2. Upload via write_raw (binary, little-endian, TYPE,8) ===")
print(f"  Sending {len(samples)} samples as {len(binary_data)} bytes of binary data")
inst.write_raw(cmd_bytes)
time.sleep(1.5)

print("\n=== 3. Did it stick? ===")
q("C1:ARWV?")

print("\n=== 4. Select waveform and enable output ===")
w("C1:BSWV WVTP,ARB")
time.sleep(0.5)
w("C1:ARWV NAME,SQTEST")
time.sleep(0.5)
q("C1:ARWV?")
q("C1:BSWV?")

input("\n  Press Enter to enable output — watch screen for square wave ...")
w("C1:OUTP ON")
time.sleep(3)
input("  Square wave or staircase? Press Enter to exit...")
w("C1:OUTP OFF")

inst.close()
rm.close()

