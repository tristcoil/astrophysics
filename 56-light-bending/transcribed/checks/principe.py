"""Verify the Principe transcriptions (Dyson, Eddington & Davidson 1920) against the paper's own arithmetic.

Reads the CSVs in 56-light-bending/transcribed/ and prints every check with n, max deviation and PASS/FLAG.
FLAG = the printed numbers do not satisfy the relation within rounding (a property of the 1920 paper, see notes);
it never means the transcription was changed.  Run:  python -I verify.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

T = str(Path(__file__).resolve().parents[1]) + '/'
rd = lambda f: pd.read_csv(T + f, comment='#', skipinitialspace=True)
report = []


def check(name, computed, printed, tol, label=None, quiet=False):
    c = np.atleast_1d(np.asarray(computed, float))
    p = np.atleast_1d(np.asarray(printed, float))
    d = c - p
    ok = np.abs(d) <= tol + 1e-9
    mx = float(np.nanmax(np.abs(d))) if d.size else 0.0
    report.append((name, d.size, mx, tol, bool(ok.all())))
    if not quiet:
        print(f'[{"PASS" if ok.all() else "FLAG"}] {name}: n={d.size} max|dev|={mx:.4g} (tol {tol})')
        for i in np.where(~ok)[0]:
            lab = label[i] if label is not None else i
            print(f'        {lab}: computed {c[i]:.4f} printed {p[i]:.4f} diff {d[i]:+.4f}')
    return d


# =============================================================== exposures (pp.315-316)
E = rd('dyson1920_principe_exposures.csv')
print('=== exposures ===')
ec = E[(E.field == 'eclipse') & (E.place == 'Principe')]
t = ec.time_h * 3600 + ec.time_m * 60 + ec.time_s
start, end = t - ec.exposure_s / 2, t + ec.exposure_s / 2
gaps = start.values[1:] - end.values[:-1]
check('eclipse plates: uniform gap between exposures (mid-times to 1 s)', gaps, np.full(gaps.size, gaps.mean()), 0.5)
print(f'   gap = {gaps.mean():.2f} s; first exposure starts {start.values[0] - (2*3600+13*60+5):.1f} s after 2h13m5s;'
      f' last ends {(2*3600+18*60+7) - end.values[-1]:.1f} s before 2h18m7s; 16 exposures + 15 gaps = '
      f'{ec.exposure_s.sum() + gaps.sum() + 1.5 + 5:.1f} s = totality 302 s')


def jd(y, m, d):
    if m <= 2:
        y -= 1; m += 12
    A = y // 100; B = 2 - A + A // 4
    return int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5


def zd(lst_h, ra_h, dec_d, lat_d):
    H = np.radians((lst_h - ra_h) * 15); d = np.radians(dec_d); p = np.radians(lat_d)
    return np.degrees(np.arccos(np.sin(p) * np.sin(d) + np.cos(p) * np.cos(d) * np.cos(H)))


fields = {'check': (14 + 12/60 + 47/3600, 20.5), 'eclipse': (4 + 19/60 + 30/3600, 21 + 43/60)}
sites = {'Oxford': (51 + 46/60, -(5*60 + 3)/3600*15), 'Principe': (1 + 40/60, (29*60 + 32)/3600*15)}
mon = {'January': 1, 'February': 2, 'May': 5}
zc = []
for _, r in E.iterrows():
    ra, dec = fields[r.field]; lat, lon = sites[r.place]
    hh = r.time_h + r.time_m/60 + r.time_s/3600
    if r.time_scale == 'Loc. Sid. T.':
        lst = hh
    else:  # astronomical G.M.T. (day starts at noon)
        y, mn, dd = r.date.split()
        j = jd(int(y), mon[mn], int(dd)) + (12 + hh)/24
        lst = ((18.697374558 + 24.06570982441908*(j - 2451545.0)) + lon/15) % 24
    zc.append(zd(lst, ra, dec, lat))
check('Approx. Z.D. recomputed (Oxford: L.S.T.; Principe: astronomical G.M.T.)', zc, E.approx_zenith_distance_deg, 0.75,
      label=list(E.ref))

# =============================================================== Tables XI-XIII
t11 = rd('dyson1920_table11_check_measures.csv'); t12 = rd('dyson1920_table12_check_residuals.csv')
t13 = rd('dyson1920_table13_check_constants.csv')
pairs = ['q1_a1', 'w1_b1', 's2_c1', 'r1_d1', 'v1_e1']
X11, Y11, S11 = t11.x.values, t11.y.values, t11.star.values
print('\n=== Table XII internal ===')
for c in ['d1x', 'd1y']:
    m = t12[[f'{p}_{c}' for p in pairs]].mean(axis=1, skipna=True)
    check(f'Table XII mean_{c} = mean of the printed residuals of the row', m, t12[f'mean_{c}'], 0.0051, label=list(t12.star))
vals = np.concatenate([t12[[f'{p}_{c}' for p in pairs]].values.ravel() for c in ['d1x', 'd1y']])
vals = vals[np.isfinite(vals)]
print(f'   n = {vals.size}; mean |v| = {np.mean(np.abs(vals)):.4f}" (printed 0".21); Peters p.e. '
      f'0.8453*sum|v|/sqrt(n(n-30)) = {0.8453*np.abs(vals).sum()/np.sqrt(vals.size*(vals.size-30)):.4f}" (printed 0".22)')
check('Table XII p.e. of one measure (Peters, n=101, 30 unknowns)', 0.8453*np.abs(vals).sum()/np.sqrt(vals.size*(vals.size-30)),
      0.22, 0.005)
check('Table XII mean residual without regard to sign', np.mean(np.abs(vals)), 0.21, 0.005)
for p in pairs:
    for c in ['dx', 'dy']:
        assert np.array_equal(np.isfinite(t11[f'{p}_{c}'].values.astype(float)), np.isfinite(t12[f'{p}_d1{c[1]}'].values)), (p, c)
print('   blank-cell pattern of Table XII = that of Table XI (checked for all 10 columns)')

print('\n=== Table XI + Table XIII constants -> Table XII ===')
allc, allp, lab = [], [], []
eff = []
for p in pairs:
    row = t13[t13.comparison == p.replace('_', '-')].iloc[0]
    for c, (k1, k2) in zip(['dx', 'dy'], [('a', 'b'), ('d', 'e')]):
        D = t11[f'{p}_{c}'].values.astype(float); m = np.isfinite(D)
        r = D[m] - row[k1]*X11[m] - row[k2]*Y11[m]
        r = (r - r.mean()) * 0.003
        allc += list(r); allp += list(t12[f'{p}_d1{c[1]}'].values[m]); lab += [f'{p} {c} star {s}' for s in S11[m]]
        # effective constants used for Table XII
        A = np.c_[X11[m], Y11[m], np.ones(m.sum())]
        sol = np.linalg.lstsq(A, D[m] - t12[f'{p}_d1{c[1]}'].values[m]/0.003, rcond=None)[0]
        lsq = np.linalg.lstsq(A, D[m], rcond=None)[0]
        eff.append((p, c, k1, k2, sol[0], sol[1], row[k1], row[k2], lsq[0], lsq[1]))
d = check('Table XII = (Table XI - a x - b y - mean) x 0".003 with Table XIII constants', allc, allp, 0.0201, label=lab)
print(f'   rms {np.sqrt(np.mean(d**2)):.4f}"; {np.sum(np.abs(d) > 0.0051)} of {d.size} differ by more than 0".005')
print('   effective constants implied by XI->XII | Table XIII | exact LS fit of Table XI:')
for p, c, k1, k2, e1, e2, t1_, t2_, l1, l2 in eff:
    print(f'   {p} {c}: {k1}={e1:+7.2f} {k2}={e2:+7.2f} | {t1_:+6.1f} {t2_:+6.1f} | {l1:+6.1f} {l2:+6.1f}')
check('effective constants (from XI->XII) = Table XIII', [e[4] for e in eff] + [e[5] for e in eff],
      [e[6] for e in eff] + [e[7] for e in eff], 0.25, label=[f'{e[0]} {e[2]}' for e in eff] + [f'{e[0]} {e[3]}' for e in eff])
check('exact LS refit of Table XI = Table XIII constants', [e[8] for e in eff] + [e[9] for e in eff],
      [e[6] for e in eff] + [e[7] for e in eff], 0.05, label=[f'{e[0]} {e[2]}' for e in eff] + [f'{e[0]} {e[3]}' for e in eff])
# LS-refit Table XII residuals
lsr, lab2 = [], []
for p in pairs:
    for c in ['dx', 'dy']:
        D = t11[f'{p}_{c}'].values.astype(float); m = np.isfinite(D)
        A = np.c_[X11[m], Y11[m], np.ones(m.sum())]
        lsr += list((D[m] - A @ np.linalg.lstsq(A, D[m], rcond=None)[0]) * 0.003)
check('Table XII = residuals of an exact LS refit of Table XI', lsr, allp, 0.0051, quiet=True)
print(f'[FLAG] Table XII = residuals of an exact LS refit of Table XI: n={len(lsr)} max|dev|={np.max(np.abs(np.array(lsr)-allp)):.3f}" '
      f'rms {np.sqrt(np.mean((np.array(lsr)-np.array(allp))**2)):.3f}"')

print('\n=== Table XIII internal ===')
r = t13[t13.comparison != 'Mean']; mean = t13[t13.comparison == 'Mean'].iloc[0]
check("b'+d' column = b' + d'", r.b_prime + r.d_prime, r.b_prime_plus_d_prime, 0.05)
check("Mean a', e', b'+d' = column means", [r.a_prime.mean(), r.e_prime.mean(), r.b_prime_plus_d_prime.mean()],
      [mean.a_prime, mean.e_prime, mean.b_prime_plus_d_prime], 0.05)
check("refraction/aberration correction equal for b and d (b'-b = d'-d)", r.b_prime - r.b, r.d_prime - r.d, 0.05)
print("   corrections a'-a:", (r.a_prime - r.a).round(1).tolist(), " b'-b:", (r.b_prime - r.b).round(1).tolist(),
      " e'-e:", (r.e_prime - r.e).round(1).tolist())
v = r.e_prime.values - r.e_prime.mean()
check("p.e. of mean e' (Peters: 0.8453 sum|v| / (n sqrt(n-1)))", 0.8453*np.abs(v).sum()/(5*np.sqrt(4)), 2.1, 0.05)
# LS-based corrected means
corr_a = (r.a_prime - r.a).values; corr_e = (r.e_prime - r.e).values; corr_b = (r.b_prime - r.b).values
lsa = np.array([e[8] for e in eff if e[1] == 'dx']); lsb = np.array([e[9] for e in eff if e[1] == 'dx'])
lsd = np.array([e[8] for e in eff if e[1] == 'dy']); lse = np.array([e[9] for e in eff if e[1] == 'dy'])
print(f"   with LS-refit constants + the printed corrections: mean a' = {np.mean(lsa + corr_a):.1f} (31.2), "
      f"mean e' = {np.mean(lse + corr_e):.1f} (37.3), mean b'+d' = {np.mean(lsb + lsd + 2*corr_b):.1f} (-1.3)")

# =============================================================== eclipse-plate reductions
P = rd('dyson1920_principe_plates.csv')
t1 = rd('dyson1920_table1_stars.csv').set_index('no')   # other agent's transcription of Table I (read only)
const = {'X-G1': (3.8, 6.9, 38.6), 'X-H1': (3.8, 8.3, 38.6), 'W-D1': (4.9, 6.5, 39.7), 'W-I2': (4.0, 9.1, 38.8),
         'U-K2': (2.8, 8.9, 37.7)}
theta = {'X-G1': 163, 'X-H1': 10, 'W-D1': 91, 'W-I2': -30, 'U-K2': 240}
print('\n=== theoretical plate constants (p.321-326) ===')
check('X-G1 refraction: eclipse (-46.5,+8.2,-27.0) - G1 (-19.1,+0.7,-28.3) = (-27.4,+7.5,+1.3)',
      [-46.5 + 19.1, 8.2 - 0.7, -27.0 + 28.3], [-27.4, 7.5, 1.3], 0.05)
check('X-G1 total: (-27.4,+7.5,+1.3) + check scale (+31.2,-0.6,+37.3) = (+3.8,+6.9,+38.6)',
      [-27.4 + 31.2, 7.5 - 0.6, 1.3 + 37.3], [3.8, 6.9, 38.6], 0.05)
check("check scale b,d = (b'+d')/2 = -1.3/2", -1.3/2, -0.6, 0.051)
for comp, (a, b, e) in const.items():
    print(f'   {comp}: implied comparison-plate refraction+aberration a={-46.5+31.2-a:+.1f} b=d={8.2-0.6-b:+.1f} e={-27.0+37.3-e:+.1f}')

xs = P[(P.coord == 'x') & P.coordinate.notna()].groupby('star').coordinate.agg(['min', 'max'])
ys = P[(P.coord == 'y') & P.coordinate.notna()].groupby('star').coordinate.agg(['min', 'max'])
assert (xs['min'] == xs['max']).all() and (ys['min'] == ys['max']).all()
X = xs['min'].to_dict(); Y = ys['min'].to_dict()
check('plate coordinates x = 14 + 10 x(Table I), y = 14 + 10 y(Table I)', [X[s] for s in X] + [Y[s] for s in Y],
      [14 + 10*t1.loc[s, 'x'] for s in X] + [14 + 10*t1.loc[s, 'y'] for s in Y], 0.005)
AL = P[P.alpha_y.notna()].groupby('star').alpha_y.agg(['min', 'max']); assert (AL['min'] == AL['max']).all(); AL = AL['min'].to_dict()
EY = P[P.Ey.notna()].groupby('star').Ey.agg(['min', 'max']); assert (EY['min'] == EY['max']).all(); EY = EY['min'].to_dict()
TQ = P[P.three_quarter_Ex.notna()].groupby('star').three_quarter_Ex.agg(['min', 'max']); assert (TQ['min'] == TQ['max']).all()
TQ = TQ['min'].to_dict()
EX = P[P.Ex.notna()].set_index('star').Ex.to_dict()
print('   alpha_y, Ey, 3/4Ex, Ex identical wherever a star is printed in several tables (asserted)')
check('Ey = 184 alpha_y', [184*AL[s] for s in AL], [EY[s] for s in AL], 0.5)
check('3/4 Ex = 0.75 Ex (U-K2 Ex)', [0.75*EX[s] for s in EX if s in TQ], [TQ[s] for s in EX if s in TQ], 0.5)
ex_full = {s: EX.get(s, TQ.get(s, np.nan)/0.75) for s in set(EX) | set(TQ)}
check('Ey x 0".003 vs Table I Principe y displacement (2 decimals)', [EY[s]*0.003 for s in EY], [t1.loc[s, 'principe_dy'] for s in EY], 0.0081)
check('Ex x 0".003 vs Table I Principe x displacement (2 decimals)', [ex_full[s]*0.003 for s in ex_full],
      [t1.loc[s, 'principe_dx'] for s in ex_full], 0.0081)
from scipy.optimize import least_squares
st = sorted(ex_full)
def resf(pp):
    out = []
    for s in st:
        ddx, ddy = X[s] - pp[0], Y[s] - pp[1]; r2 = ddx*ddx + ddy*ddy
        out += [10*ddy/r2 - AL[s], 10*ddx/r2 - ex_full[s]/184]
    return out
sun = least_squares(resf, [19.5, 11.5]).x
fit = np.array(resf(sun))
check('alpha_y and Ex/184 = 10 (r - r_sun)/|r - r_sun|^2 (one Sun position)', fit, np.zeros(fit.size), 0.0075)
print(f'   fitted Sun centre: x = {sun[0]:.3f}, y = {sun[1]:.3f} (plate units); 1.75*15.78/50/0.003 = {1.75*15.78/50/0.003:.1f} (kappa 184)')
check('relative y displacement 5-3 = 1".2 and 11-3 = 0".8 (p.321)', [(EY[3]-EY[5])*0.003, (EY[3]-EY[11])*0.003], [1.2, 0.8], 0.05)

dev_rows = []
for comp in ['X-G1', 'X-H1', 'W-D1', 'W-I2', 'U-K2']:
    a, b, e = const[comp]; d_ = b; th = theta[comp]
    print(f'\n=== {comp}: a={a} b=d={b} e={e} theta={th} ===')
    sx = P[(P.comparison == comp) & (P.coord == 'x')].set_index('star')
    sy = P[(P.comparison == comp) & (P.coord == 'y')].set_index('star')
    stars = list(sx.index); assert stars == list(sy.index)
    xx = np.array([X[s] for s in stars]); yy = np.array([Y[s] for s in stars])
    if sx.term_ax_or_dx.notna().any():
        check(f'{comp} x: printed terms 3.8x, 6.9y', list(sx.term_ax_or_dx) + list(sx.term_by_or_ey), list(a*xx) + list(b*yy), 0.51,
              label=[f'ax {s}' for s in stars] + [f'by {s}' for s in stars])
        check(f'{comp} y: printed terms 6.9x, 38.6y', list(sy.term_ax_or_dx) + list(sy.term_by_or_ey), list(d_*xx) + list(e*yy), 0.6,
              label=[f'dx {s}' for s in stars] + [f'ey {s}' for s in stars])
        check(f'{comp} x: Delta1x = Delta x - printed terms (exact)', sx.delta - sx.term_ax_or_dx - sx.term_by_or_ey, sx.delta1, 0)
        check(f'{comp} y: Delta1y = Delta y - printed terms (exact)', sy.delta - sy.term_ax_or_dx - sy.term_by_or_ey, sy.delta1, 0)
    check(f'{comp} x: Delta1x = Delta x - (a x + b y)', sx.delta - a*xx - b*yy, sx.delta1, 1.0, label=stars)
    check(f'{comp} y: Delta1y = Delta y - (d x + e y)', sy.delta - d_*xx - e*yy, sy.delta1, 1.0, label=stars)
    check(f'{comp} x: orientation column = theta y', sx.orientation_term, th*yy, 0.51 if comp not in ('X-G1', 'W-D1') else 1.7, label=stars)
    check(f'{comp} y: orientation column = -theta x', sy.orientation_term, -th*xx, 0.51, label=stars)
    if comp != 'U-K2':
        tq = np.array([TQ[s] for s in stars])
        check(f'{comp} x: Delta2x = Delta1x - 3/4Ex (exact)', sx.delta1 - tq, sx.delta2, 0)
        check(f'{comp} x: Delta3x = Delta2x - theta y (exact)', sx.delta2 - sx.orientation_term, sx.delta3, 0)
        check(f'{comp} x: Delta4x = Delta3x - 1/4Ex', sx.delta3 - np.array([ex_full[s] for s in stars])/4, sx.delta4, 0.75, label=stars)
        check(f'{comp} y: Delta3y = Delta1y + theta x (exact)', sy.delta1 - sy.orientation_term, sy.delta3, 0)
        check(f'{comp} y: Delta4y = Delta3y - Ey (exact)', sy.delta3 - sy.Ey, sy.delta4, 0)
        th_ls = np.polyfit(yy, sx.delta2.values.astype(float), 1)[0]
        check(f"{comp}: printed theta = LS solution of Delta2x = c + theta y", th_ls, th, 0.5)
    else:
        check(f'{comp} x: Delta4x = Delta1x - 240y - Ex (exact)', sx.delta1 - sx.orientation_term - sx.Ex, sx.delta4, 0)
        check(f'{comp} y: Delta4y = Delta1y + 240x - Ey (exact)', sy.delta1 - sy.orientation_term - sy.Ey, sy.delta4, 0, label=stars)
    check(f'{comp} x: Resid = Delta4x - round(mean Delta4x)', sx.delta4 - np.round(sx.delta4.mean()), sx.resid, 0, label=stars)
    check(f'{comp} y: Resid = Delta4y - round(mean Delta4y)', sy.delta4 - np.round(sy.delta4.mean()), sy.resid, 0, label=stars)
    if comp == 'X-G1':
        check('X-G1: mean Delta4x printed -5942', sx.delta4.mean(), -5942, 0.5)
    al = np.array([AL.get(s, EY[s]/184) for s in stars]); rr = sy.resid.values.astype(float)
    n = len(stars); c22 = np.sum(al**2) - al.sum()**2/n; rhs = np.sum(al*rr) - al.sum()*rr.sum()/n
    dev_rows.append((comp, al.sum(), np.sum(al**2), rr.sum(), np.sum(al*rr), c22, rhs, rhs/c22))

print('\n=== normal equations, delta kappa, limb deflections (pp.323-325) ===')
printed_ne = {'X-G1': (2.83, 4.83, -1, 64, 3.23, 64, 20, 1.94), 'X-H1': (2.83, 4.83, 1, -105, 3.23, -105, -33, 1.44),
              'W-D1': (3.42, 5.21, -1, -62, 2.87, -61, -21, 1.55), 'W-I2': (3.42, 5.21, 2, -24, 2.87, -25, -9, 1.67)}
for comp, s1, s2, r1, r2, c22, rhs, dk in dev_rows:
    if comp not in printed_ne:
        continue
    pn = printed_ne[comp]
    check(f'{comp}: sum alpha, sum alpha^2, eliminated coefficient', [s1, s2, c22], [pn[0], pn[1], pn[4]], 0.005)
    check(f'{comp}: rhs1 = sum resid', r1, pn[2], 0.5)
    check(f'{comp}: rhs2 = sum alpha*resid (from printed residuals)', r2, pn[3], 0.5)
    check(f'{comp}: eliminated rhs from the printed normal equations', pn[3] - pn[0]*pn[2]/5, pn[5], 0.5)
    check(f'{comp}: eliminated rhs from the printed residuals', rhs, pn[5], 0.5)
    check(f'{comp}: delta kappa = eliminated rhs / coefficient (printed numbers)', pn[5]/pn[4], pn[6], 0.5)
    check(f'{comp}: limb deflection 1.75 (184 + dk)/184 with printed eliminated equation', 1.75*(184 + pn[5]/pn[4])/184, pn[7], 0.005)
    check(f'{comp}: limb deflection from the printed residuals (exact LS)', 1.75*(184 + dk)/184, pn[7], 0.005)
check('mean of the four limb deflections', np.mean([1.94, 1.44, 1.55, 1.67]), 1.65, 0.005)
check('final result 1.65 - 0.04', 1.65 - 0.04, 1.61, 0.005)

print('\n=== U-K2 joint solution (p.326) ===')
sx = P[(P.comparison == 'U-K2') & (P.coord == 'x')].set_index('star'); sy = P[(P.comparison == 'U-K2') & (P.coord == 'y')].set_index('star')
stars = list(sx.index); xx = np.array([X[s] for s in stars]); yy = np.array([Y[s] for s in stars])
A = np.zeros((10, 4)); A[:5, 0] = 1; A[:5, 1] = yy; A[:5, 3] = sx.Ex.values/184
A[5:, 2] = 1; A[5:, 1] = -xx; A[5:, 3] = sy.Ey.values/184
bb = np.r_[sx.resid.values, sy.resid.values].astype(float)
sol = np.linalg.lstsq(A, bb, rcond=None)[0]; vv = bb - A @ sol; Ninv = np.linalg.inv(A.T @ A)
pe_peters = 0.8453*np.abs(vv).sum()/np.sqrt(10*6)*np.sqrt(Ninv[3, 3])*1.75/184
check('U-K2: delta theta (LS of x and y residuals for dc, dtheta, df, dkappa)', sol[1], 2, 0.5)
check('U-K2: delta kappa', sol[3], 121, 0.5)
check('U-K2: limb deflection', 1.75*(184 + sol[3])/184, 2.90, 0.005)
check('U-K2: probable error (Peters)', pe_peters, 0.87, 0.005)
print(f'   LS: dtheta {sol[1]:+.2f}, dkappa {sol[3]:+.1f}, deflection {1.75*(184+sol[3])/184:.2f}", p.e. {pe_peters:.2f}" '
      f'(Gauss {0.6745*np.sqrt(vv@vv/6)*np.sqrt(Ninv[3,3])*1.75/184:.2f}"), weight of dkappa {1/Ninv[3,3]:.2f}')

print('\n=== p.328 residual table ===')
R = rd('dyson1920_p328_principe_residuals.csv').set_index('star')
cmap = {'G': 'X-G1', 'H': 'X-H1', 'D': 'W-D1', 'I': 'W-I2'}
for c in ['x', 'y']:
    cc, pp, ll = [], [], []
    for k, comp in cmap.items():
        s = P[(P.comparison == comp) & (P.coord == c)].set_index('star').resid
        for star in R.index:
            if star in s.index:
                cc.append(0.003*s[star]); pp.append(R.loc[star, f'{c}_{k}']); ll.append(f'{c}_{k} {star}')
            else:
                assert not np.isfinite(R.loc[star, f'{c}_{k}'])
    check(f'p.328 {c} residuals = 0.003 x Resid of pp.323-325', cc, pp, 0.0051, label=ll)
    cols = [f'{c}_{k}' for k in cmap]; full = R[cols].notna().all(axis=1)
    check(f'p.328 {c} Mean = mean of the four', R.loc[full, cols].mean(axis=1), R.loc[full, f'{c}_mean'], 0.0051)
xu = np.concatenate([P[(P.comparison == comp) & (P.coord == 'x')].resid.values for comp in cmap.values()]).astype(float)
yu = np.concatenate([P[(P.comparison == comp) & (P.coord == 'y')].resid.values for comp in cmap.values()]).astype(float)
check('p.328 average y-residual (0".003 x mean|Resid|, 20 values)', 0.003*np.mean(np.abs(yu)), 0.22, 0.005)
check('p.328 p.e. of y (Peters, 20 values, 4 unknowns)', 0.003*0.8453*np.abs(yu).sum()/np.sqrt(20*16), 0.21, 0.005)
check('p.328 p.e. of x (Peters, 20 values, 4 unknowns)', 0.003*0.8453*np.abs(xu).sum()/np.sqrt(20*16), 0.25, 0.005)
check('p.328 p.e. of kappa = 0.21/sqrt(3)', 0.21/np.sqrt(3), 0.12, 0.005)
check('p.328 p.e. at the limb = 0.12 x 50/15.78', 0.21/np.sqrt(3)*50/15.78, 0.38, 0.005)
check('p.328 total p.e. sqrt(0.38^2/2 + 0.14^2)', np.hypot(0.38/np.sqrt(2), 0.14), 0.30, 0.005)

print('\n=== error budget p.327 ===')
k2l = 1.75/184
for plate, comp, pr in [('X', 'X-G1', (0.68, 0.45, 4.0, 0.18, 0.35, 0.092)), ('W', 'W-D1', (0.68, 0.22, 3.5, 0.08, 0.38, 0.029))]:
    sy = P[(P.comparison == comp) & (P.coord == 'y')].set_index('star'); sx = P[(P.comparison == comp) & (P.coord == 'x')].set_index('star')
    stars = list(sy.index); xx = np.array([X[s] for s in stars]); yy = np.array([Y[s] for s in stars])
    al = np.array([AL[s] for s in stars]); alc = al - al.mean(); Sw = np.sum(alc**2)
    per_th = np.sum(alc*xx)/Sw*k2l
    pe_th = (0.20/0.003)/np.sqrt(np.sum((yy - yy.mean())**2))
    q = np.array([ex_full[s] for s in stars])/4
    slope = np.polyfit(yy, q, 1)[0]
    check(f'p.327 plate {plate}: 10 units of e -> deflection', np.sum(alc*10*yy)/Sw*k2l, pr[0], 0.0051)
    check(f'p.327 plate {plate}: 10 units of orientation -> deflection', 10*per_th, pr[1], 0.0051)
    check(f'p.327 plate {plate}: p.e. 0".20 in x -> p.e. of orientation', pe_th, pr[2], 0.05)
    check(f'p.327 plate {plate}: ... -> p.e. of deflection', pe_th*per_th, pr[3], 0.0051)
    check(f'p.327 plate {plate}: p.e. 0".20 in y -> p.e. of deflection', (0.20/0.003)/np.sqrt(Sw)*k2l, pr[4], 0.0051)
    check(f'p.327 plate {plate}: neglected 1/4 Ex in the orientation -> deflection', slope*per_th, pr[5], 0.0011)
check('p.327 scale p.e.: 2.1/10 x 0".68', 0.21*0.68, 0.14, 0.005)
check("p.327 largest e' (45.0) instead of mean (37.3): (45.0-37.3)/10 x 0.68", 0.77*0.68, 0.52, 0.005)
check('p.327 mean correction (0.029 + 0.092)/2', (0.029 + 0.092)/2, 0.06, 0.005)

print('\n=== orientation: printed theta vs least squares, effect on the deflections ===')
defl_ls = []
for comp in ['X-G1', 'X-H1', 'W-D1', 'W-I2']:
    sx = P[(P.comparison == comp) & (P.coord == 'x')].set_index('star'); sy = P[(P.comparison == comp) & (P.coord == 'y')].set_index('star')
    stars = list(sx.index); xx = np.array([X[s] for s in stars]); yy = np.array([Y[s] for s in stars])
    al = np.array([AL[s] for s in stars]); alc = al - al.mean()
    th_ls = np.polyfit(yy, sx.delta2.values.astype(float), 1)[0]
    d4 = sy.delta1.values + th_ls*xx - sy.Ey.values
    dk = np.sum(alc*(d4 - d4.mean()))/np.sum(alc**2)
    defl_ls.append(1.75*(184 + dk)/184)
    print(f'   {comp}: theta printed {theta[comp]:+d}, LS {th_ls:+.1f}; deflection with LS theta {defl_ls[-1]:.2f}"')
print(f'   mean with LS orientations: {np.mean(defl_ls):.3f}" (paper 1.65)')

print('\n=== SUMMARY ===')
for name, n, mx, tol, ok in report:
    print(f'{"PASS" if ok else "FLAG"} | n={n:3d} | max dev {mx:.4g} | tol {tol} | {name}')
print(f'{sum(r[4] for r in report)} PASS, {sum(not r[4] for r in report)} FLAG, {sum(r[1] for r in report)} values compared')
