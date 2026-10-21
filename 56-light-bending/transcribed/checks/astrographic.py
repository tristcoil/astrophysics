"""Independent internal checks (T3-audit). Usage: checks.py <dir with the CSVs>"""
from pathlib import Path
import sys, math
import numpy as np, pandas as pd
D = (sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parents[1])).rstrip('/') + '/'
rd = lambda f, **k: pd.read_csv(D + f, comment='#', skipinitialspace=True, **k)
t9 = rd('dyson1920_table9_astrographic.csv', dtype={'eclipse_plate': str, 'comparison_plate': str})
dy = rd('dyson1920_p311_astrographic_dy.csv')
ape = rd('dyson1920_p311_alpha_plus_e.csv')
rows = t9[t9.eclipse_plate != 'Mean'].copy()
mean = t9[t9.eclipse_plate == 'Mean'].iloc[0]
rows['eclipse_plate'] = rows.eclipse_plate.astype(int)
print("== Table IX: alpha_limb vs alpha_rev x 12.3 x 50/15.8")
k = 12.3 * 50 / 15.8
bad = 0
for _, r in rows.iterrows():
    calc = r.alpha_rev * k
    # interval of factor consistent with rounding of alpha (+-0.0005) and limb (+-0.005)
    ok = abs(calc - r.alpha_limb_arcsec) <= 0.005 + 1e-9
    ok2 = abs(calc - r.alpha_limb_arcsec) <= 0.005 + 0.0005 * k + 1e-9
    if not ok: bad += 1
    print(f"  plate {r.eclipse_plate:>2} {r.comparison_plate:>5} alpha {r.alpha_rev:+.3f} -> {calc:+.4f}  printed {r.alpha_limb_arcsec:+.2f}  diff {r.alpha_limb_arcsec-calc:+.4f} {'ok' if ok else ('within alpha rounding' if ok2 else 'FAIL')}")
print(f"  factor k = {k:.4f}; rows not exact after rounding: {bad}")
# range of k consistent with all rows (non-zero alpha)
lo, hi = 0, 1e9
for _, r in rows.iterrows():
    if r.alpha_rev == 0: continue
    lo = max(lo, (r.alpha_limb_arcsec - 0.005) / r.alpha_rev)
    hi = min(hi, (r.alpha_limb_arcsec + 0.005) / r.alpha_rev)
print(f"  factor range reproducing every printed limb value from printed alpha: [{lo:.4f}, {hi:.4f}]")
print("== Table IX means")
for c in ['d_rev', 'e_rev', 'alpha_rev', 'alpha_limb_arcsec']:
    s = rows[c].sum(); n = len(rows)
    pm = mean[c]
    print(f"  {c}: sum {s:+.3f} / {n} = {s/n:+.5f}  printed {pm}")
# alternative: per-eclipse-plate means (16 plates)
g = rows.groupby('eclipse_plate')[['e_rev', 'alpha_rev', 'alpha_limb_arcsec']].mean()
print("  per-plate (16) means:", {c: round(g[c].mean(), 5) for c in g.columns})
print("== p.311 combination coefficients")
eq = {5: (-0.160, -1.107, -0.789), 4: (0.334, 0.472, 1.336), 3: (0.348, 0.360, 1.574), 6: (0.587, 1.099, 0.726), 10: (0.860, 1.321, 0.589)}
A = [eq[5][i] + 0.5 * (eq[6][i] + eq[10][i]) for i in range(3)]
B = [eq[4][i] + eq[3][i] for i in range(3)]
C = [B[i] - A[i] for i in range(3)]
print("  (1)+((4)+(5))/2: d %.4f e %.4f a %.4f (printed +0.564 +0.103 -0.131)" % tuple(A))
print("  (2)+(3):         d %.4f e %.4f a %.4f (printed +0.682 +0.832 +2.910)" % tuple(B))
print("  difference:      d %.4f e %.4f a %.4f (printed 0.118 0.729 3.041)" % tuple(C))
print("  e/alpha printed-coefficient ratio 0.729/3.041 = %.5f" % (0.729 / 3.041))
# plate 15 variant: star 10 missing -> use (4) alone in place of mean of (4),(5)
A15 = [eq[5][i] + eq[6][i] for i in range(3)]
C15 = [B[i] - A15[i] for i in range(3)]
print("  plate 15 [(2)+(3)-(1)-(4)]: d %.4f e %.4f a %.4f -> e/alpha %.5f (printed 0.283)" % (C15[0], C15[1], C15[2], C15[1] / C15[2]))
print("== recompute alpha + 0.240e per plate")
DY = dy.pivot(index='eclipse_plate', columns='star', values='Dy_rev')
dmap = {}
for _, r in rows.iterrows():
    dmap.setdefault(r.eclipse_plate, []).append((r.comparison_plate, r.d_rev))
res = []
for _, r in ape.iterrows():
    p = int(r.eclipse_plate)
    y = DY.loc[p]
    opts = dmap[p]
    cand = [(c, d) for c, d in opts]
    if len(opts) > 1:
        cand.append(('mean', np.mean([d for _, d in opts])))
    for c, d in cand:
        if p == 15:
            rhs = y[3] + y[4] - y[5] - y[6]
            val = (rhs - 0.255 * d) / 2.973
            val_b = (rhs - C15[0] * d) / C15[2]
        else:
            rhs = y[3] + y[4] - y[5] - 0.5 * (y[6] + y[10])
            val = (rhs - 0.118 * d) / 3.041
            val_b = (rhs - C[0] * d) / C[2]
        diff = r.value_rev - val
        print(f"  plate {p:>2} d({c:>5})={d:+.3f} rhs {rhs:+.4f} -> {val:+.5f} (unrounded coeffs {val_b:+.5f}) printed {r.value_rev:+.3f} coeff {r.coeff_e:.3f} diff {diff:+.5f} {'ok' if abs(diff) <= 0.0005 + 1e-9 else 'OFF'}")
print("== p.312 means")
s = ape.value_rev.sum()
print(f"  sum of 16 values {s:.4f}, mean {s/len(ape):.5f} (printed 0.0435); mean coeff {ape.coeff_e.mean():.5f} (printed 243 [0.243])")
a1 = 0.0435 - 0.243 * 0.082
a2 = 0.0435 - 0.243 * 0.020
print(f"  alpha = 0.0435 - 0.243*0.082 = {a1:.5f} r -> limb {a1*k:.3f}\" ; rounded 0.024 -> {0.024*k:.3f}\" (printed 0.024 r = 0.93\")")
print(f"  alpha = 0.0435 - 0.243*0.020 = {a2:.5f} r -> limb {a2*k:.3f}\" ; rounded 0.039 -> {0.039*k:.3f}\" (printed 0.039 r = 1.52\")")
print("== n_stars and comparisons")
print(rows[['eclipse_plate', 'comparison_plate', 'n_stars']].to_string(index=False))
