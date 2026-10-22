"""Mori-route strip size at a window [-L,L] (L = HALF-width), using the team's Arb transcription of Mori's scalar
constants (GPPDiscovery5/s66/mori_consts.py, Mori VII/VIII zenodo 23265728). For given shift s (>= lambda_max of the
prime operator C_p, measured by cp_norm.py) and m = 1/2 it scans the cutoff Omega and finds the smallest Legendre
degree N at which Mori's hypotheses hold (C>0, u0E0<=C, u1E1<=C, n<m, c>2) and alpha+beta b^2 < 2^-target.
Pole weights: u0 = 0.999 C/E0, u1 = 0.999 C/E1 (saturating, as in Mori VIII where u0E0 ~ 0.13684 vs C ~ 0.13699).
Also prints the Astra necessary bound Omega > 2 pi e^{s+m}. Status: ESTIMATE (constants CERTIFIED arithmetic of
Mori's formulas; s and the cure factor are NUMERICAL inputs). usage: mori_strip.py L s "Om1,Om2,.." target_bits
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../s66'))
from mori_consts import consts
from flint import arb
L = sys.argv[1]; s = arb(sys.argv[2]); Oms = [int(x) for x in sys.argv[3].split(',')]; tb = int(sys.argv[4])
Lv = arb(L.split('/')[0])/int(L.split('/')[1]) if '/' in L else arb(L); m = arb(1)/2
print('Astra bound Omega >', (2*arb.pi()*(s+m).exp()).str(6))
rows = []
for Om in Oms:
    c = arb(Om)*Lv; C = (arb(Om)/(2*arb.pi())).log() - 1/arb(Om) - s - m
    E0 = (c+1)/(arb.pi()*c**2); E1 = (3*c**2+6*c+2)/(arb.pi()*c**3)
    if not C > 0: rows.append({"Omega": Om, "C": C.str(5), "ok": False}); print(rows[-1]); continue
    u0 = float((C/E0).mid())*0.999; u1 = float((C/E1).mid())*0.999
    lo, hi = int(float(c.mid())), int(4*float(c.mid()))+400
    def good(N):
        out, ok = consts(Lv, Om, s, m, u0, u1, N)
        return ok and bool(out['alpha_plus_beta_b2'] < arb(2)**(-tb)), out
    if not good(hi)[0]: rows.append({"Omega": Om, "C": C.str(5), "ok": False, "note": "no N<=%d" % hi}); print(rows[-1]); continue
    while hi-lo > 10:
        mid = (lo+hi)//2
        if good(mid)[0]: hi = mid
        else: lo = mid
    out = good(hi)[1]
    rows.append({"Omega": Om, "C": C.str(5), "u0": round(u0), "u1": round(u1), "N_min": hi, "N_over_OmL": hi/float(c.mid()),
                 "team_rule_N=1.671*Om*L": round(1.671*float(c.mid())), "mlog2_delta": out['mlog2_dc'].str(6), "ok": True})
    print(rows[-1], flush=True)
json.dump({"L": L, "s": sys.argv[2], "target_bits": tb, "rows": rows}, open('runs/tail/mori_strip_%s.json' % L.replace('/', '_'), 'w'), indent=1)
