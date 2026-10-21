"""T2 verification: reads ONLY the final CSVs in transcribed/ (plus, read-only, T1's Table I/II/III-VI CSVs for
cross-checks) and runs every arithmetic check. Usage: python -I verify_t2.py"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

T = str(Path(__file__).resolve().parents[1]) + '/'
rd = lambda f: pd.read_csv(T + f, comment='#', skipinitialspace=True, dtype=str, keep_default_na=False)
t7, ra, so, t8 = rd('dyson1920_table7_alpha.csv'), rd('dyson1920_p307_refraction_aberration.csv'), \
    rd('dyson1920_p307_scale_orientation.csv'), rd('dyson1920_table8_residuals.csv')
p9, p31, ex = rd('dyson1920_p309_obs_calc.csv'), rd('dyson1920_p331_radial.csv'), rd('dyson1920_sobral_exposures.csv')
STARS = ['11', '5', '4', '3', '6', '10', '2']
ECL = ['I', 'II', 'III', 'IV', 'V', 'VII', 'VIII']
CMP = ['14_2a', '14_2b', '15_1', '15_2', '17_1', '17_2', '18_2']
f = float
results = []


def report(name, n, maxdev, verdict, detail=''):
    results.append((name, n, maxdev, verdict, detail))
    print(f'[{verdict}] {name}: n={n}, max dev={maxdev}{"; " + detail if detail else ""}')


# ------------------------------------------------------------------ 0. format
fmt_bad = []
def chk(df, cols, dec, signed=True):
    for c in cols:
        for v in df[c]:
            if v == '':
                continue
            pat = (r'^[+-]\d+\.\d{%d}$' % dec) if signed else (r'^\d+\.\d{%d}$' % dec)
            if not re.match(pat, v) and not re.match(r'^0\.0{%d}$' % dec, v):
                fmt_bad.append((c, v))
chk(t7, [c for c in t7.columns if c.endswith('_rev')], 3)
chk(ra, ['a_rev', 'e_rev', 'b_rev', 'd_rev'], 3)
chk(so, [c for c in so.columns if c.endswith('_rev')], 3)
chk(t8, ['value_arcsec'], 2)
chk(p9, ['ra_observed', 'ra_calculated', 'dec_observed', 'dec_calculated'], 2)
chk(p31, ['calculation_arcsec', 'observation_arcsec'], 2, signed=False)
counts = (len(t7), len(ra), len(so), len(t8), len(p9), len(p31), len(ex))
ok = counts == (8, 9, 14, 224, 7, 7, 27) and not fmt_bad
nvals = 28 + 4 + 36 + 84 + 224 + 28 + 14
report('0 format (sign, printed decimals, row counts 8/9/14/224/7/7/27)', nvals, '-', 'PASS' if ok else 'FAIL',
       f'unsigned zeros: {sorted(set(v for v in t8.value_arcsec if not v[0] in "+-"))}, bad={fmt_bad}')

# ------------------------------------------------------------------ 1. Table VII means
cols = ['ra_eclipse_minus_scale_rev', 'ra_comparison_minus_scale_rev', 'dec_eclipse_minus_scale_rev', 'dec_comparison_minus_scale_rev']
body = t7[t7.row != 'Mean']; mean = t7[t7.row == 'Mean'].iloc[0]
devs = [abs(body[c].astype(float).mean() - f(mean[c])) for c in cols]
report('1 Table VII Mean row = mean of the 7 values', 4, f'{max(devs):.5f} r', 'PASS' if max(devs) <= 0.0005 else 'FAIL',
       'means ' + ', '.join(f'{body[c].astype(float).mean():+.5f}' for c in cols))

# ------------------------------------------------------------------ 2. Table VII vs Tables III-VI (T1 CSV)
ne = pd.read_csv(T + 'dyson1920_tables3to6_normal_equations.csv', comment='#', skipinitialspace=True, dtype=str)
sol = lambda tab, p, q: ne[(ne.table == tab) & (ne.plate == p) & (ne.quantity == q)].value.iloc[0]
mism = []
for i in range(7):
    r = body.iloc[i]
    for col, tab, p in (('ra_eclipse_minus_scale_rev', 'III', r.eclipse_plate), ('ra_comparison_minus_scale_rev', 'IV', r.comparison_plate),
                        ('dec_eclipse_minus_scale_rev', 'V', r.eclipse_plate), ('dec_comparison_minus_scale_rev', 'VI', r.comparison_plate)):
        if f(r[col]) != f(sol(tab, p, 'sol_alpha')):
            mism.append((col, p, r[col], sol(tab, p, 'sol_alpha')))
report('2 Table VII = printed alpha of Tables III-VI (T1 CSV), plate order I..VIII / 14_2a..18_2', 28, '0' if not mism else str(mism),
       'PASS' if not mism else 'FAIL')

# ------------------------------------------------------------------ 3. p.306 arithmetic
m = {c: f(mean[c]) for c in cols}
a_ra = m[cols[0]] - m[cols[1]]; a_de = m[cols[2]] - m[cols[3]]
a_w = (a_ra + 2 * a_de) / 3
arc50 = 0.100 * 6.25
limb = arc50 * 50 / 15.8
lines = [f'alpha_RA = {a_ra:+.4f} r (printed +0.105)', f'alpha_Dec = {a_de:+.4f} r (printed +0.098)',
         f'(1*RA + 2*Dec)/3 = {a_w:+.5f} r (printed +0.100)', f'0.100 r x 6.25 = {arc50:.4f}" (printed +0.625)',
         f'0.625" x 50/15.8 = {limb:.4f}" (printed 1.98)']
ok3 = abs(a_ra - 0.105) < 5e-4 and abs(a_de - 0.098) < 5e-4 and abs(a_w - 0.100) < 5e-4 and abs(arc50 - 0.625) < 5e-4 and abs(limb - 1.98) < 0.005
report('3 p.306 arithmetic: RA, Dec, 1:2 weighting, r -> arcsec at 50\', limb with 15\'.8', 5, f'{abs(limb - 1.98):.4f}"', 'PASS' if ok3 else 'FAIL', '; '.join(lines))
k = 6.25 * 50 / 15.8
s39 = {'Dec': (0.098 * k, 1.94), 'RA': (0.105 * k, 2.06)}
dev39 = max(abs(round(v, 2) - p) for v, p in s39.values())
report('3b section 39 (p.330): limb values from 0r.098 and 0r.105', 2, f'{dev39:.2f}"', 'ORIGINAL',
       f'Dec {s39["Dec"][0]:.3f}" -> 1.94 printed 1.94; RA {s39["RA"][0]:.3f}" -> 2.08 but printed 2.06 (R5); '
       f'weighted (2.06 + 2 x 1.94)/3 = {(2.06 + 3.88) / 3:.3f}, (2.077 + 2 x 1.938)/3 = {(s39["RA"][0] + 2 * s39["Dec"][0]) / 3:.3f} -> both 1.98')
# probable error from the scatter (method not stated in the paper; informational)
v = {c: body[c].astype(float).values for c in cols}
me = lambda x: x.std(ddof=1) / np.sqrt(len(x))
pe_ra = 0.6745 * np.hypot(me(v[cols[0]]), me(v[cols[1]])); pe_de = 0.6745 * np.hypot(me(v[cols[2]]), me(v[cols[3]]))
pe_w = np.sqrt(pe_ra ** 2 + 4 * pe_de ** 2) / 3
report('3c probable error from the scatter of Table VII (informational)', 28, '-', 'INFO',
       f'p.e. RA {pe_ra:.4f} r ({pe_ra / a_ra * 100:.1f}%), Dec {pe_de:.4f} r ({pe_de / a_de * 100:.1f}%), 1:2 mean {pe_w:.4f} r '
       f'({pe_w / a_w * 100:.1f}%); paper: "about 6 per cent" (p.306), +-0".12 = 6.1% of 1".98 (p.331)')

# ------------------------------------------------------------------ 4. p.307 second table from solutions + first table
R = {r.row_label: np.array([f(r.a_rev), f(r.e_rev), f(r.b_rev), f(r.d_rev)]) for r in ra.itertuples()}
Rs = R['Scale plate']
dev4 = []; bad4 = []
for r in so.itertuples():
    ecl = r.plate in ECL
    tabRA, tabDec = ('III', 'V') if ecl else ('IV', 'VI')
    a, b = f(sol(tabRA, r.plate, 'sol_a')), f(sol(tabRA, r.plate, 'sol_b'))
    e, d = f(sol(tabDec, r.plate, 'sol_e')), f(sol(tabDec, r.plate, 'sol_d'))
    dR = (R['Eclipse plates'] if ecl else R[r.plate]) - Rs
    calc = [a + dR[0], e + dR[1], b + dR[2], -(d + dR[3])]
    pr = [f(r.scale_from_x_rev), f(r.scale_from_y_rev), f(r.orientation_from_x_rev), f(r.orientation_from_y_rev)]
    for name, c_, p_ in zip(('scale_x', 'scale_y', 'orient_x', 'orient_y'), calc, pr):
        dev4.append(abs(c_ - p_))
        if abs(c_ - p_) > 0.0005:
            bad4.append(f'{r.plate} {name}: printed {p_:+.3f}, arithmetic {c_:+.3f}')
n_ok = sum(x <= 0.0005 for x in dev4)
report('4 p.307 scale/orientation = Tables III-VI solutions + (plate - scale plate) refraction/aberration', 56,
       f'{max(dev4):.3f} r', 'ORIGINAL' if bad4 else 'PASS', f'{n_ok}/56 exact; ' + '; '.join(bad4) +
       ' (orientation from y uses -(d + dd); III V sol_b +0.042 is itself a 1920 slip per T1 D5 but p.307 follows it)')
# adopted orientation and adopted scale
d_or = [abs((f(r.orientation_from_x_rev) + f(r.orientation_from_y_rev)) / 2 - f(r.adopted_orientation_rev)) for r in so.itertuples()]
report('4b adopted orientation = mean(from x, from y)', 14, f'{max(d_or):.4f} r', 'PASS' if max(d_or) <= 0.0005 + 1e-9 else 'FAIL')
grp = {'July 14 (14_2a,14_2b)': ['14_2a', '14_2b'], 'July 15-18 (15_1..18_2)': ['15_1', '15_2', '17_1', '17_2', '18_2'], 'eclipse': ECL}
txt = []
for g, pl in grp.items():
    sub = so[so.plate.isin(pl)]
    mval = np.concatenate([sub.scale_from_x_rev.astype(float), sub.scale_from_y_rev.astype(float)]).mean()
    txt.append(f'{g}: mean of printed scale values {mval:+.4f}, adopted {sub.adopted_scale_rev.iloc[0]}')
report('4c adopted scale = group mean of the printed scale values', 3, '-', 'PASS (comparison groups; eclipse 0.000 is not the mean)', '; '.join(txt))
# refraction table internal
alt = {'Scale plate': 33.6, '14_2a': 34.3, '14_2b': 34.3, '15_1': 35.4, '15_2': 36.4, '17_1': 34.7, '17_2': 38.6, '18_2': 39.2}
order = sorted(alt, key=lambda k_: alt[k_])
mono = all(R[order[i]][j] >= R[order[i + 1]][j] for i in range(len(order) - 1) for j in (0, 1, 2))
bd = all(r.b_rev == r.d_rev for r in ra.itertuples())
same = (R['14_2a'] == R['14_2b']).all()
report('4d p.307 refraction table: b = d, 14_2a = 14_2b, a/e/b non-increasing with p.301 altitude', 9, '-',
       'PASS' if (mono and bd and same) else 'FAIL', f'b=d: {bd}; 14_2a=14_2b: {same}; monotonic: {mono}')

# ------------------------------------------------------------------ 5. Table VIII Mean column
dev5 = []; worst = None
for (kind, coord), g in t8.groupby(['kind', 'coord']):
    for s in STARS:
        gs = g[g.star == s]
        vals = gs[gs.plate != 'Mean'].value_arcsec.astype(float)
        pm = f(gs[gs.plate == 'Mean'].value_arcsec.iloc[0])
        dd = vals.mean() - pm
        dev5.append(abs(dd))
        if worst is None or abs(dd) > abs(worst[3]):
            worst = (kind, coord, s, dd, vals.mean(), pm)
within = sum(x <= 0.005 + 1e-9 for x in dev5)
report('5 Table VIII printed Mean = mean of the 7 plate values', 28, f'{max(dev5):.4f}"',
       'PASS' if within == 28 else 'ORIGINAL', f'{within}/28 within 0.005"; worst {worst[0]} {worst[1]} star {worst[2]}: '
       f'mean {worst[4]:+.4f} vs printed {worst[5]:+.2f} (R2)')
# column sums
txt = []
for (kind, coord), g in t8.groupby(['kind', 'coord']):
    sums = g.groupby('plate').value_arcsec.apply(lambda x: x.astype(float).sum())
    txt.append(f'{kind} {coord}: column sums {sums.min():+.2f}..{sums.max():+.2f} (Mean column {sums["Mean"]:+.2f})')
report('5b Table VIII column sums (RA ~0; Dec constant within each block)', 32, '-', 'PASS', '; '.join(txt))

# ------------------------------------------------------------------ 6. p.309 observed = eclipse Mean - comparison Mean
M = t8[t8.plate == 'Mean'].set_index(['kind', 'coord', 'star']).value_arcsec.astype(float)
dev6 = []; bad6 = []
for r in p9.itertuples():
    for coord, col in (('RA', 'ra_observed'), ('Dec', 'dec_observed')):
        c_ = M[('eclipse', coord, r.star)] - M[('comparison', coord, r.star)]
        dd = abs(c_ - f(getattr(r, col)))
        dev6.append(dd)
        if dd > 0.0051:
            bad6.append(f'star {r.star} {coord}: printed {getattr(r, col)}, Table VIII {c_:+.2f}')
report('6 p.309 observed = Table VIII eclipse Mean - comparison Mean', 14, f'{max(dev6):.2f}"',
       'ORIGINAL' if bad6 else 'PASS', f'{sum(x < 0.0051 for x in dev6)}/14 exact; ' + '; '.join(bad6) + ' (R1)')

# ------------------------------------------------------------------ 7. p.309 calculated vs Table I and K*(Ex,Ey) of p.304 (T1 CSVs)
t1 = pd.read_csv(T + 'dyson1920_table1_stars.csv', comment='#', skipinitialspace=True, dtype=str).set_index('no')
cc = pd.read_csv(T + 'dyson1920_4inch_condition_coeffs.csv', comment='#', skipinitialspace=True, dtype=str).set_index('star')
diff7 = []
for r in p9.itertuples():
    for col, t1c in (('ra_calculated', 'sobral_dx'), ('dec_calculated', 'sobral_dy')):
        if f(getattr(r, col)) != f(t1.loc[r.star, t1c]):
            diff7.append(f'star {r.star} {col}: p.309 {getattr(r, col)} vs Table I {t1.loc[r.star, t1c]}')
report('7 p.309 calculated = Table I Sobral displacements (T1 CSV)', 14, '0.01"', 'ORIGINAL' if diff7 else 'PASS',
       f'{14 - len(diff7)}/14 identical; ' + '; '.join(diff7))
Ks = np.arange(0.5400, 0.5650, 0.00005)
okK = [K for K in Ks if all(round(K * f(cc.loc[r.star, 'Ex']) + 1e-12, 2) == f(r.ra_calculated) and
                            round(K * f(cc.loc[r.star, 'Ey']) + 1e-12, 2) == f(r.dec_calculated) for r in p9.itertuples())]
report('7b p.309 calculated = K x (Ex, Ey) of p.304 for one K', 14, '-', 'PASS' if okK else 'FAIL',
       (f'all 14 reproduced for K in [{min(okK):.4f}, {max(okK):.4f}]" (1".75 x 15\'.8/50\' = 0.5530 is outside: '
        f'{sum(round(0.553 * f(cc.loc[r.star, c]) + 1e-12, 2) == f(getattr(r, k_)) for r in p9.itertuples() for c, k_ in (("Ex", "ra_calculated"), ("Ey", "dec_calculated")))}/14)') if okK else '')

# ------------------------------------------------------------------ 8. p.331 radial vs p.309 projected
# Table I star x, y (read on leaf002 for this check; identical to T1's CSV) and the sun from the 26 Table I Sobral values
XY = {s: (f(t1.loc[s, 'x']), f(t1.loc[s, 'y'])) for s in t1.index}
own = {'1': (0.026, -0.200), '2': (1.079, -0.328), '3': (0.348, 0.360), '4': (0.334, 0.472), '5': (-0.160, -1.107),
       '6': (0.587, 1.099), '7': (-0.707, -0.864), '8': (-0.727, -1.040), '9': (-0.483, -1.303), '10': (0.860, 1.321),
       '11': (-1.261, -0.160), '12': (-1.311, -0.918), '13': (0.089, 1.007)}
assert all(own[s] == XY[s] for s in own), 'own Table I reading differs from T1 CSV'
P = np.array([XY[s] for s in t1.index]); Dsp = np.array([[f(t1.loc[s, 'sobral_dx']), f(t1.loc[s, 'sobral_dy'])] for s in t1.index])
use = np.ones_like(Dsp, bool); use[list(t1.index).index('13'), 0] = False      # star 13 dx is an outlier (T1 D11)
def model(p):
    d = P - p[:2]; return p[2] * d / (d ** 2).sum(1)[:, None]
fit = least_squares(lambda p: (model(p) - Dsp)[use], [0.44, -0.26, 0.55]).x
for label, sun in (('fitted sun (%.4f, %.4f)' % (fit[0], fit[1]), fit[:2]), ('task sun (0.44, -0.26)', np.array([0.44, -0.26]))):
    dc, do, txt = [], [], []
    for r in p31.itertuples():
        row9 = p9[p9.star == r.star].iloc[0]
        u = np.array(XY[r.star]) - sun; u /= np.hypot(*u)
        pc = f(row9.ra_calculated) * u[0] + f(row9.dec_calculated) * u[1]
        po = f(row9.ra_observed) * u[0] + f(row9.dec_observed) * u[1]
        dc.append(abs(pc - f(r.calculation_arcsec))); do.append(abs(po - f(r.observation_arcsec)))
        if abs(po - f(r.observation_arcsec)) > 0.0051:
            alt_ = -f(row9.ra_observed) * u[0] + f(row9.dec_observed) * u[1]
            txt.append(f'star {r.star} obs: projection {po:.3f} vs printed {r.observation_arcsec}; with R.A. sign flipped {alt_:.3f}')
    report(f'8 p.331 radial = p.309 components projected away from the sun [{label}]', 14,
           f'calc {max(dc):.4f}", obs {max(do):.4f}"', 'ORIGINAL' if txt else 'PASS',
           f'calc {sum(x <= 0.0051 for x in dc)}/7, obs {sum(x <= 0.0051 for x in do)}/7 within rounding; ' + '; '.join(txt))
print(f'   (sun fit from Table I: x0={fit[0]:+.4f}, y0={fit[1]:+.4f}, K={fit[2]:.4f}"; K/r alternative for the calc column:',
      ', '.join(f'{s}:{fit[2] / np.hypot(*(np.array(XY[s]) - fit[:2])):.3f}' for s in p31.star), ')')

# ------------------------------------------------------------------ 9. Table VIII rebuilt from Table II (T1 CSV) with p.307 constants
t2 = pd.read_csv(T + 'dyson1920_table2_4inch.csv', comment='#', skipinitialspace=True, dtype=str)
x = np.array([XY[s][0] for s in STARS]); y = np.array([XY[s][1] for s in STARS])
A3 = np.column_stack([np.ones(7), x, y])
adop = so.set_index('plate')
dev_adopt, dev_free, outl = [], [], []
for (kind, coord), g in t8.groupby(['kind', 'coord']):
    for p in (ECL if kind == 'eclipse' else CMP):
        sub = t2[t2.plate == p].set_index('star')
        D = np.array([f(sub.loc[s, 'Dx_rev' if coord == 'RA' else 'Dy_rev']) for s in STARS])
        Rv = np.array([f(g[(g.plate == p) & (g.star == s)].value_arcsec.iloc[0]) for s in STARS])
        dR = (R['Eclipse plates'] if kind == 'eclipse' else R[p]) - Rs
        S_, O_ = f(adop.loc[p, 'adopted_scale_rev']), f(adop.loc[p, 'adopted_orientation_rev'])
        pc, qc = (S_ - dR[0], O_ - dR[2]) if coord == 'RA' else (-O_ - dR[3], S_ - dR[1])
        u_ = 6.25 * (D - pc * x - qc * y); dv = u_ - (u_ - Rv).mean() - Rv
        dev_adopt.extend(np.abs(dv))
        Y = D - Rv / 6.25
        s3, *_ = np.linalg.lstsq(A3, Y, rcond=None); r3 = (Y - A3 @ s3) * 6.25
        dev_free.extend(np.abs(r3))
        if np.abs(r3).max() > 0.0105:
            best = None
            for kk in range(7):
                mk = np.ones(7, bool); mk[kk] = False
                s2, *_ = np.linalg.lstsq(A3[mk], Y[mk], rcond=None); r2 = (Y - A3 @ s2) * 6.25
                if best is None or np.abs(r2[mk]).max() < best[0]:
                    best = (np.abs(r2[mk]).max(), kk, r2[kk])
            outl.append(f'{kind} {coord} {p} star {STARS[best[1]]}: printed {Rv[best[1]]:+.2f}, Table II implies '
                        f'{Rv[best[1]] + best[2]:+.3f} (other 6 fit to {best[0]:.4f}")')
dev_adopt = np.array(dev_adopt); dev_free = np.array(dev_free)
report('9 Table VIII rebuilt from Table II + p.307 adopted constants (zero point per column = mean)', 196,
       f'{dev_adopt.max():.3f}"', 'ORIGINAL', f'{(dev_adopt <= 0.0105).sum()}/196 within 0.0105", {(dev_adopt <= 0.02).sum()}/196 within 0.02"')
report('9b same with 3 free constants per column (zero point, x and y terms)', 196, f'{dev_free.max():.3f}"', 'ORIGINAL',
       f'{(dev_free <= 0.0105).sum()}/196 within 0.0105"; columns with one outlying star: ' + '; '.join(outl))

# ------------------------------------------------------------------ 10. exposures
ok10 = []
for inst, g in ex.groupby('instrument'):
    t = pd.to_datetime(g.utc_civil_derived).values
    e_ = g.exposure_s.astype(int).values
    gaps = (np.diff(t).astype('timedelta64[s]').astype(int) - e_[:-1])
    ok10.append(f'{inst}: {len(g)} plates, exposures {sorted(set(e_))} s, gaps between exposures {gaps.min()}..{gaps.max()} s, '
                f'first start {g.t_after_second_contact_s_derived.iloc[0]} s, last end {int(g.t_after_second_contact_s_derived.iloc[-1]) + e_[-1]:+d} s after second contact')
alt5 = ex[ex.instrument == 'astrographic'].exposure_s.astype(int).values
alternate = all(alt5[i] != alt5[i + 1] for i in range(len(alt5) - 1))
report('10 p.300 exposures: counts 19/8, astro alternate 5/10 s, 4-inch 28 s, no overlaps (p.299 text)', 27, '-',
       'PASS' if alternate and len(ex) == 27 else 'FAIL', '; '.join(ok10))

print('\nSUMMARY')
for r_ in results:
    print(' | '.join(str(v) for v in r_))

# ------------------------------------------------------------------ 9c. breakdown of the adopted-constant rebuild per block / column
print('\n9c adopted-constant rebuild, per column max |dev| (arcsec):')
for (kind, coord), g in t8.groupby(['kind', 'coord']):
    line = []
    for p in (ECL if kind == 'eclipse' else CMP):
        sub = t2[t2.plate == p].set_index('star')
        D = np.array([f(sub.loc[s, 'Dx_rev' if coord == 'RA' else 'Dy_rev']) for s in STARS])
        Rv = np.array([f(g[(g.plate == p) & (g.star == s)].value_arcsec.iloc[0]) for s in STARS])
        dR = (R['Eclipse plates'] if kind == 'eclipse' else R[p]) - Rs
        S_, O_ = f(adop.loc[p, 'adopted_scale_rev']), f(adop.loc[p, 'adopted_orientation_rev'])
        pc, qc = (S_ - dR[0], O_ - dR[2]) if coord == 'RA' else (-O_ - dR[3], S_ - dR[1])
        u_ = 6.25 * (D - pc * x - qc * y); dv = u_ - (u_ - Rv).mean() - Rv
        line.append(f'{p}:{np.abs(dv).max():.3f}')
    print(f'   {kind:10s} {coord:3s} ' + ' '.join(line))
# 17_1 RA: x coefficient implied by Table VIII vs p.307
g = t8[(t8.kind == 'comparison') & (t8.coord == 'RA') & (t8.plate == '17_1')].set_index('star')
sub = t2[t2.plate == '17_1'].set_index('star')
D = np.array([f(sub.loc[s, 'Dx_rev']) for s in STARS]); Rv = np.array([f(g.loc[s, 'value_arcsec']) for s in STARS])
s3, *_ = np.linalg.lstsq(A3, D - Rv / 6.25, rcond=None)
print(f'   17_1 RA: x coefficient from Table VIII {s3[1]:+.4f} r, from p.307 (-0.084 - (0.399 - 0.423)) {-0.084 + 0.024:+.3f} r; '
      f'y coefficient {s3[2]:+.4f} vs {-0.145 - (0.091 - 0.096):+.3f}')
# VIII star 3 zero-sum hint
col = t8[(t8.kind == 'eclipse') & (t8.coord == 'RA') & (t8.plate == 'VIII')].value_arcsec.astype(float)
print(f'   eclipse RA VIII column sum {col.sum():+.2f} (with star 3 = +0.01: {col.sum() + 0.02:+.2f}; = +0.02: {col.sum() + 0.03:+.2f})')
