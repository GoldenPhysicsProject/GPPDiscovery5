"""CERTIFIED (Arb) combined lower bound at the END of each character's effective prime-free window (test-space inclusion
then covers every smaller h)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s17_bounds import E_comb, C0, C1
from s17_build import CHARS, arb
for cn, f, e, ch in CHARS:
    n = {"chi3": 2, "chi4": 3, "chi5 quadratic": 2, "chi5 quartic": 4}[cn.split(" (")[0]]
    h = arb(n).log()/2
    print(f"{cn:45s} h_end = log({n})/2 = {h.str(10)} : LB = {(E_comb(h, e)+arb(f).log()-(C0 if e == 0 else C1)).str(10)}")
