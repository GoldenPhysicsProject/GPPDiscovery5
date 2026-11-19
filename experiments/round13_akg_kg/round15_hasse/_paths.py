"""Path setup for Round 15 scripts (AUDIT15 s7): run from any working directory.
Puts this folder, round14b_euler and round13_akg_kg on sys.path and chdirs here, so relative inputs/outputs
(../clues, ../round14b_euler/noghosts.json, *.json) resolve.  A scratch copy elsewhere: set GPP_R13=/path/to/round13_akg_kg."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
R13 = os.environ.get('GPP_R13') or os.path.join(HERE, '..')
for p in (os.path.join(R13, 'round14b_euler'), R13, HERE):
    if p not in sys.path: sys.path.insert(0, p)
os.chdir(HERE)
