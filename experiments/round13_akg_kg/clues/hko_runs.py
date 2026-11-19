"""Post-process hko_interlace.json: runs of >=3 consecutive same-piece zeros (a 'double interlacing defect')
vs the fakes' off-line zeros (test data from ../epstein_zeros_lowband.json and 14B/14F for DH)."""
import json
d = json.load(open('hko_interlace.json'))
off = {'EPS': [z[1] for z in json.load(open('../epstein_zeros_lowband.json'))['offline_right_half']],
       'DH': [85.6993484854, 114.16]}
res = {}
for k in ('EPS', 'DH'):
    ev = sorted([(x, 'A') for x in d[k]['zeros_A']] + [(x, 'B') for x in d[k]['zeros_B']])
    runs, i = [], 0
    while i < len(ev):
        j = i
        while j+1 < len(ev) and ev[j+1][1] == ev[i][1]: j += 1
        if j-i+1 >= 3: runs.append((round(ev[i][0], 3), round(ev[j][0], 3), ev[i][1], j-i+1))
        i = j+1
    T = max(ev)[0]
    singles = sum(1 for a in range(len(ev)-1) if ev[a][1] == ev[a+1][1])
    hit = []
    for g in off[k]:
        if g > T: continue
        dist = min(min(abs(g-r[0]), abs(g-r[1]), 0 if r[0] <= g <= r[1] else 9e9) for r in runs) if runs else None
        hit.append((g, round(dist, 3)))
    covered = [r for r in runs if any(r[0]-1.5 <= g <= r[1]+1.5 for g in off[k])]
    res[k] = dict(T=T, n_A=len(d[k]['zeros_A']), n_B=len(d[k]['zeros_B']), adjacent_same_piece_pairs=singles,
                  runs_ge3=runs, offline_heights=[g for g in off[k] if g <= T],
                  offline_to_nearest_run=hit, runs_within_1p5_of_an_offline_zero=len(covered))
    print(k, json.dumps(res[k]))
json.dump(res, open('hko_runs.json', 'w'), indent=1)
