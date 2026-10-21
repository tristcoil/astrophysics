"""Check 10 of the notes: per-plate least squares for alpha (p.312 '0".99'), and the effect-on-result claims."""
from pathlib import Path
import sys
import numpy as np, pandas as pd
D = (sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parents[1])).rstrip('/') + '/'
rd = lambda f, **k: pd.read_csv(D + f, comment='#', skipinitialspace=True, **k)
t9 = rd('dyson1920_table9_astrographic.csv', dtype={'eclipse_plate': str, 'comparison_plate': str})
t9 = t9[t9.eclipse_plate != 'Mean'].copy(); t9['eclipse_plate'] = t9.eclipse_plate.astype(int)
dy = rd('dyson1920_p311_astrographic_dy.csv'); ape = rd('dyson1920_p311_alpha_plus_e.csv')
DY = dy.pivot(index='eclipse_plate', columns='star', values='Dy_rev')
eq = {5: (-0.160, -1.107, -0.789), 4: (0.334, 0.472, 1.336), 3: (0.348, 0.360, 1.574), 6: (0.587, 1.099, 0.726), 10: (0.860, 1.321, 0.589)}
k = 12.3 * 50 / 15.8
def solve(p, d=None):
    rows, rhs = [], []
    for s in (5, 4, 3, 6, 10):
        y = DY.loc[p, s]
        if np.isnan(y): continue
        cd, ce, ca = eq[s]
        if d is None:
            rows.append([ca, ce, 1.0, cd]); rhs.append(y)
        else:
            rows.append([ca, ce, 1.0]); rhs.append(y - cd * d)
    sol = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)[0]
    return sol
dsets = {}
for _, r in t9.iterrows():
    dsets.setdefault(r.eclipse_plate, []).append(r.d_rev)
for label, pick in [('first-listed', lambda L: L[0]), ('second-listed', lambda L: L[-1]), ('mean d', lambda L: float(np.mean(L)))]:
    al = [solve(p, pick(dsets[p]))[0] for p in DY.index]
    print(f"d fixed ({label}): mean alpha {np.mean(al):.5f} r = {np.mean(al)*k:.3f}\"")
al = [solve(p)[0] for p in DY.index]
print(f"d free: mean alpha {np.mean(al):.5f} r = {np.mean(al)*k:.3f}\"")
# effect on result: recomputed values with 18_3 for plate 7 and mean d for plate 12
def val(p, d):
    y = DY.loc[p]
    if p == 15:
        return (y[3] + y[4] - y[5] - y[6] - 0.255 * d) / 2.973
    return (y[3] + y[4] - y[5] - 0.5 * (y[6] + y[10]) - 0.118 * d) / 3.041
for lab, d7, d12 in [('7:18_3, 12:mean', -0.054, -0.0555), ('7:mean, 12:mean', -0.024, -0.0555)]:
    vals = []
    for p in DY.index:
        d = dsets[p][0]
        if p == 7: d = d7
        if p == 12: d = d12
        vals.append(val(p, d))
    m = np.mean(vals)
    a1 = m - 0.243 * 0.082; a2 = m - 0.243 * 0.020
    print(f"recomputed ({lab}): mean {m:.5f}; alpha(e=0.082) {a1:.5f} r = {a1*k:.3f}\"; alpha(e=0.020) {a2:.5f} r = {a2*k:.3f}\"")
