"""Independent checks of Harvey 1979 Tables II/III and the combined result."""
from pathlib import Path
import sys, itertools
import numpy as np, pandas as pd
D = (sys.argv[1] if len(sys.argv) > 1 else str(Path(__file__).resolve().parents[1])).rstrip('/') + '/'
t2 = pd.read_csv(D + 'harvey1979_table2.csv', comment='#', skipinitialspace=True, dtype={'plate': str})
t3 = pd.read_csv(D + 'harvey1979_table3.csv', comment='#', skipinitialspace=True)
t1 = pd.read_csv(D + 'harvey1979_table1.csv', comment='#', skipinitialspace=True)
print(t1.to_string(index=False))
def wmean(x, s):
    w = 1 / s**2
    m = (w * x).sum() / w.sum(); e = 1 / np.sqrt(w.sum())
    chi2 = (w * (x - m)**2).sum()
    return m, e, chi2
res = {}
for inst in ('4-inch', '13-inch'):
    g = t2[(t2.instrument == inst) & (t2.plate != 'Mean')]
    pm = t2[(t2.instrument == inst) & (t2.plate == 'Mean')].iloc[0]
    m, e, chi2 = wmean(g['lambda'].values, g['lambda_se'].values)
    res[inst] = (m, e)
    n = len(g)
    print(f"{inst}: n={n} weighted mean {m:.4f} +- {e:.4f} (printed {pm['lambda']:.2f} +- {pm['lambda_se']:.2f}); chi2 {chi2:.2f} for {n-1} dof; x1.75 -> {1.75*m:.4f} +- {1.75*e:.4f}")
    # attainable range under +-0.005 rounding of each lambda and se (extreme search: corners for se, lambda shift)
    rng = np.random.default_rng(1)
    ms, es = [], []
    x0, s0 = g['lambda'].values, g['lambda_se'].values
    for _ in range(200000):
        x = x0 + rng.uniform(-0.005, 0.005, n); s = s0 + rng.uniform(-0.005, 0.005, n)
        mm, ee, _ = wmean(x, s); ms.append(mm); es.append(ee)
    # exact extremes for error: all se at -0.005 / +0.005
    emin = wmean(x0, s0 - 0.005)[1]; emax = wmean(x0, s0 + 0.005)[1]
    print(f"   attainable x1.75 (random search): mean {1.75*min(ms):.4f}..{1.75*max(ms):.4f}; error exact {1.75*emin:.4f}..{1.75*emax:.4f}")
print("Table III rows:")
print(t3.to_string(index=False))
z = t3.set_index('determination')
a = z.loc['4-inch plates measured on the Zeiss']; b = z.loc['Astrographic plates measured on the Zeiss']
x = np.array([a.displacement_arcsec, b.displacement_arcsec]); s = np.array([a.error_arcsec, b.error_arcsec])
m, e, chi2 = wmean(x, s)
print(f"combined Zeiss (printed T3): {m:.4f} +- {e:.4f}; chi2 {chi2:.3f} (1 dof) -> external {e*np.sqrt(chi2):.4f}; printed 1.87 +- 0.13")
print(f"  (1.87-1.75)/0.13 = {(1.87-1.75)/0.13:.3f}; ({m:.4f}-1.75)/{e:.4f} = {(m-1.75)/e:.3f}; (1.87-1.75)/0.105 = {0.12/0.105:.3f}")
print(f"  weighted mean of errors sum(s*w)/sum(w) = {(s/s**2).sum()/(1/s**2).sum():.4f}")
# combining lambda means first (unrounded, and printed)
for lab, (m4, e4), (m13, e13) in [('unrounded', res['4-inch'], res['13-inch']), ('printed', (1.09, 0.06), (0.88, 0.19))]:
    mm, ee, c2 = wmean(np.array([m4, m13]), np.array([e4, e13]))
    print(f"  lambda means combined ({lab}): {mm:.4f} +- {ee:.4f} -> {1.75*mm:.4f} +- {1.75*ee:.4f}; external {1.75*ee*np.sqrt(c2):.4f}")
allx = t2[t2.plate != 'Mean']
mm, ee, c2 = wmean(allx['lambda'].values, allx['lambda_se'].values)
print(f"  all 14 plates pooled: {mm:.4f} +- {ee:.4f} -> {1.75*mm:.4f} +- {1.75*ee:.4f}; chi2 {c2:.2f}/13 -> external x1.75 {1.75*ee*np.sqrt(c2/13):.4f}")
print(f"  1.75 x printed means: {1.75*1.09:.4f} {1.75*0.06:.4f} {1.75*0.88:.4f} {1.75*0.19:.4f}")
print(f"  Dyson PE->SE: 0.12/0.6745 = {0.12/0.6745:.4f}")
