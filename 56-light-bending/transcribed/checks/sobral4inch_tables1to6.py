"""Verification of the T1 transcription (Dyson, Eddington & Davidson 1920, Sobral 4-inch chain).

Reads ONLY the six CSV files in transcribed/ and tests them against the paper's own arithmetic:
  0. format: explicit signs, printed precision
  A. Table I: displacements radial from one Sun position with magnitude k/r, per site (fit x0, y0, k)
  B. p.304 condition coefficients: x, y = Table I; (Ex, Ey) = (s - s0)/|s - s0|^2; k*Ex vs Table I
  C. Table II: offset convention (value - offset small and positive); column sums vs printed rhs1
  D. LHS of Tables III-VI rebuilt from the p.304 coefficients (normal equations + Gauss reduction)
  E. RHS of Tables III-VI rebuilt from Table II minus offsets: every printed line and solution,
     E1 with the values as printed, E2 with the three Table II values the 1920 computation used
  F. internal arithmetic of every printed table (each line from the printed lines above it)
  G. single-digit / transposition / sign search for every column that fails E1
  H. hand-computation model: products rounded to 3 decimals before summing (rhs2-rhs4)
  I. cross-checks with later pages (typed here only for checking, read from the full-resolution scans):
     Table VII (p.306) alpha values and means; p.307 scale/orientation; Table VIII (p.308) R.A. residuals
usage: python -I verify_t1.py      (prints a report, writes nothing)
Companion checks (need astropy / tesseract): gaia_xy_check.py, sun_check.py, ocr_check.py
"""
from pathlib import Path
import re
from decimal import Decimal, ROUND_HALF_UP

import numpy as np
import pandas as pd

T = str(Path(__file__).resolve().parents[1]) + '/'
rd = dict(comment='#', skipinitialspace=True)
STR = dict(dtype=str, keep_default_na=False, **rd)

t1 = pd.read_csv(T + 'dyson1920_table1_stars.csv', **rd)
t1s = pd.read_csv(T + 'dyson1920_table1_stars.csv', **STR)
cc = pd.read_csv(T + 'dyson1920_4inch_condition_coeffs.csv', **rd).set_index('star')
ccs = pd.read_csv(T + 'dyson1920_4inch_condition_coeffs.csv', **STR).set_index('star')
t2s = pd.read_csv(T + 'dyson1920_table2_4inch.csv', **STR)
offs = pd.read_csv(T + 'dyson1920_table2_offsets.csv', **STR).set_index('plate')
nes = pd.read_csv(T + 'dyson1920_tables3to6_normal_equations.csv', **STR)
lhs = pd.read_csv(T + 'dyson1920_tables3to6_lhs.csv', **STR)
nev = nes.set_index(['table', 'plate', 'quantity'])['value'].astype(float)

STARS = [11, 5, 4, 3, 6, 10, 2]
ECL = ['I', 'II', 'III', 'IV', 'V', 'VII', 'VIII']
CMP = ['14_2a', '14_2b', '15_1', '15_2', '17_1', '17_2', '18_2']
QRA = ['rhs%d' % i for i in range(1, 11)] + ['sol_alpha', 'sol_a', 'sol_b']
QDE = ['rhs%d' % i for i in range(1, 11)] + ['sol_alpha', 'sol_e', 'sol_d']
TABS = (('III', ECL, 'Dx'), ('IV', CMP, 'Dx'), ('V', ECL, 'Dy'), ('VI', CMP, 'Dy'))
TOL = 0.002
# the three Table II values the 1920 reduction used where they differ from the printed ones (section E2)
USED = {('I', 'Dx', 2): -0.798, ('14_2a', 'Dx', 2): -0.346}
USED_OFF = {('VIII', 'Dy'): -1.332}
SUMMARY = []


def rep(name, n, maxdev, verdict, extra=''):
    SUMMARY.append((verdict, name, n, maxdev, extra))
    print(f'  -> [{verdict}] {name}: n={n}, max|dev|={maxdev:.4f} {extra}')


def table2(plate, comp, used=False):
    sub = t2s[t2s.plate == plate].set_index('star').loc[[str(s) for s in STARS]]
    v = sub[comp + '_rev'].astype(float).to_numpy().copy()
    o = float(offs.loc[plate, comp + '_offset'])
    if used:
        for i, s in enumerate(STARS):
            v[i] = USED.get((plate, comp, s), v[i])
        o = USED_OFF.get((plate, comp), o)
    return v, o


# ================================================================== 0. format
print('=== 0. Format: explicit sign and printed precision in every numeric cell ===')
p3, p2 = re.compile(r'^[+-]\d\.\d{3}$'), re.compile(r'^[+-]\d\.\d{2}$')
cells = ([(v, p3) for c in ('x', 'y') for v in t1s[c]] +
         [(v, p2) for c in ('sobral_dx', 'sobral_dy', 'principe_dx', 'principe_dy') for v in t1s[c]] +
         [(v, re.compile(r'^\d\.\d$')) for v in t1s['photog_mag']] +
         [(v, p3) for c in ('x', 'y', 'Ex', 'Ey') for v in ccs[c]] +
         [(v, p3) for c in ('Dx_rev', 'Dy_rev') for v in t2s[c]] +
         [(v, p3) for c in ('Dx_offset', 'Dy_offset') for v in offs[c]] +
         [(v, p3) for v in nes['value']] +
         [(v, p3) for c in ('coef1', 'coef2', 'coef3', 'coef4') for v in lhs[c] if v != ''])
bad = [v for v, p in cells if not p.match(v)]
rep('format of all numeric cells', len(cells), 0.0, 'PASS' if not bad else 'FAIL', str(bad) if bad else '')
counts = {'table1': len(t1s), 'cond': len(ccs), 'table2': len(t2s), 'offsets': len(offs), 'ne': len(nes), 'lhs': len(lhs)}
ok = counts == {'table1': 13, 'cond': 7, 'table2': 98, 'offsets': 14, 'ne': 364, 'lhs': 40}
rep('row counts (13, 7, 98, 14, 364, 40)', 6, 0.0, 'PASS' if ok else 'FAIL', str(counts))

# ================================================================== A. Table I
print('\n=== A. Table I: radial displacement pattern k (s - s0)/|s - s0|^2 for each site ===')
xy = t1[['x', 'y']].to_numpy(float)
nos = t1['no'].tolist()
inside = bool(np.all(np.abs(xy) <= 1.4))
rep('Table I coordinates inside the 4-inch/astrographic field (|x|,|y| <= 1.4 units = 70 arcmin)', 26,
    float(np.abs(xy).max()), 'PASS' if inside else 'FAIL', '(max |coordinate| shown)')


def model(p, s):
    dx, dy = s[:, 0] - p[0], s[:, 1] - p[1]
    r2 = dx * dx + dy * dy
    return np.concatenate([p[2] * dx / r2, p[2] * dy / r2])


def fit(s, d, w=None, p0=(0.44, -0.26, 0.55), npar=3):
    w = np.ones(len(d)) if w is None else w
    p = np.array(p0, float)
    for _ in range(200):
        full = lambda q: model(np.r_[q, 1.0] if npar == 2 else q, s)
        f0 = full(p)
        J = np.empty((len(d), npar))
        for k in range(npar):
            dp = np.zeros(npar); dp[k] = 1e-7
            J[:, k] = (full(p + dp) - f0) / 1e-7
        st = np.linalg.lstsq(J * w[:, None], (d - f0) * w, rcond=None)[0]
        p += st
        if np.max(np.abs(st)) < 1e-13:
            break
    return p, d - full(p)


SUN = {}
for site, outlier in (('sobral', (13, 'dx')), ('principe', (8, 'dx'))):
    d = np.concatenate([t1[site + '_dx'].to_numpy(float), t1[site + '_dy'].to_numpy(float)])
    p_all, r_all = fit(xy, d)
    w = np.ones(26); io = nos.index(outlier[0]) + (0 if outlier[1] == 'dx' else 13); w[io] = 0
    p, r = fit(xy, d, w)
    SUN[site] = p
    keep = w > 0
    print(f'  {site:8s} all 26: x0={p_all[0]:+.4f} y0={p_all[1]:+.4f} k={p_all[2]:.4f}"  rms {np.sqrt(np.mean(r_all**2)):.4f}" '
          f'max {np.abs(r_all).max():.4f}" (star {nos[int(np.argmax(np.abs(r_all))) % 13]})')
    print(f'  {site:8s} w/o star {outlier[0]} {outlier[1]}: x0={p[0]:+.4f} y0={p[1]:+.4f} k={p[2]:.4f}" '
          f'(= 1.75" x {p[2] / 1.75 * 50:.2f}\'/50\')  rms(25) {np.sqrt(np.mean(r[keep]**2)):.4f}" max(25) '
          f'{np.abs(r[keep]).max():.4f}"; star {outlier[0]} {outlier[1]}: printed {d[io]:+.2f}, pattern {d[io] - r[io]:+.4f}')
    dv = np.stack([t1[site + '_dx'], t1[site + '_dy']], 1).astype(float)
    rv = xy - p[:2]
    cosang = np.sum(dv * rv, 1) / np.linalg.norm(dv, axis=1) / np.linalg.norm(rv, axis=1)
    rep(f'Table I {site} displacements vs radial k/r pattern (25 values; rounding alone gives rms 0.0029")', 25,
        float(np.abs(r[keep]).max()), 'PASS' if np.abs(r[keep]).max() <= 0.0075 else 'FAIL',
        f'rms {np.sqrt(np.mean(r[keep]**2)):.4f}"; all 13 vectors point away from the Sun (min cos {cosang.min():.4f})')
    rep(f'Table I {site} star {outlier[0]} {outlier[1]} (printed {d[io]:+.2f})', 1, float(abs(r[io])),
        'ORIGINAL', f'differs from the paper\'s own pattern ({d[io] - r[io]:+.3f}); kept as printed')
sh = SUN['principe'][:2] - SUN['sobral'][:2]
print(f'  Sun moved Sobral -> Principe by dx={50 * sh[0]:+.2f}\' dy={50 * sh[1]:+.2f}\' (Horizons, sun_check.py: '
      f'+5.12\', +0.90\' between the mid-totality times)')

# ================================================================== B. condition coefficients
print('\n=== B. p.304 condition coefficients ===')
t1i = t1.set_index('no')
dev = [abs(cc.loc[s, c] - t1i.loc[s, c]) for s in STARS for c in ('x', 'y')]
rep('x, y printed on p.304 equal Table I', 14, max(dev), 'PASS' if max(dev) == 0 else 'FAIL')
s7 = cc.loc[STARS, ['x', 'y']].to_numpy(float)
E = np.concatenate([cc.loc[STARS, 'Ex'].to_numpy(float), cc.loc[STARS, 'Ey'].to_numpy(float)])
p2, r2 = fit(s7, E, p0=(0.44, -0.26), npar=2)
rep('Ex, Ey = (s - s0)/|s - s0|^2 for one Sun position', 14, float(np.abs(r2).max()),
    'PASS' if np.abs(r2).max() <= 0.0015 else 'FAIL',
    f'Sun x0={p2[0]:+.4f} y0={p2[1]:+.4f} vs Table I Sobral fit {SUN["sobral"][0]:+.4f} {SUN["sobral"][1]:+.4f}')
pred = SUN['sobral'][2] * E
obs = np.concatenate([t1i.loc[STARS, 'sobral_dx'].to_numpy(float), t1i.loc[STARS, 'sobral_dy'].to_numpy(float)])
rep('k_Sobral x (Ex, Ey) vs Table I Sobral displacements of the 7 stars', 14, float(np.abs(obs - pred).max()),
    'PASS' if np.abs(obs - pred).max() <= 0.0075 else 'FAIL')

# ================================================================== C. Table II offsets and sums
print('\n=== C. Table II: offsets below the line ===')
neg, sums = [], []
for plate in ECL + CMP:
    for comp in ('Dx', 'Dy'):
        v, o = table2(plate, comp)
        u = np.round(v - o, 3)
        if np.any(u < 0):
            neg.append(f'{plate} {comp} star {STARS[int(np.argmin(u))]} {u.min():+.3f}')
        tab = ('III' if plate in ECL else 'IV') if comp == 'Dx' else ('V' if plate in ECL else 'VI')
        sums.append((tab, plate, comp, round(nev[(tab, plate, 'rhs1')] - u.sum(), 4)))
rep('value - printed offset >= 0 ("small and positive")', 196, 0.0, 'ORIGINAL' if neg else 'PASS',
    'exception: ' + '; '.join(neg) + ' (printed offset -1.322 vs column minimum -1.332)' if neg else '')
exact = [s for s in sums if s[3] == 0]
other = [s for s in sums if s[3] != 0]
print('  columns whose printed rhs1 != sum(value - offset):', other)
rep('printed rhs1 = sum(Table II value - offset), exact', 28, max(abs(s[3]) for s in sums),
    'PASS' if len(other) == 4 else 'FAIL', f'{len(exact)}/28 exact; the 4 others explained in sections E2/G')

# ================================================================== D. LHS
print('\n=== D. Left-hand sides of Tables III-VI rebuilt from the p.304 coefficients ===')
x = cc.loc[STARS, 'x'].to_numpy(float); y = cc.loc[STARS, 'y'].to_numpy(float)
A = {'Dx': np.column_stack([np.ones(7), y, x, cc.loc[STARS, 'Ex'].to_numpy(float)]),     # c, b, a, alpha
     'Dy': np.column_stack([np.ones(7), x, y, cc.loc[STARS, 'Ey'].to_numpy(float)])}     # f, d, e, alpha


def chain(N, r):
    """Gauss reduction in the printed order; returns the 10 lhs rows, rhs1..rhs10 and (alpha, u3, u2)."""
    N1 = N[1:, 1:] - np.outer(N[0, 1:], N[0, 1:]) / N[0, 0]
    N2 = N1[1:, 1:] - np.outer(N1[0, 1:], N1[0, 1:]) / N1[0, 0]
    N3 = N2[1, 1] - N2[0, 1] ** 2 / N2[0, 0]
    r1 = r[1:] - N[0, 1:] * r[0] / N[0, 0]
    r2 = r1[1:] - N1[0, 1:] * r1[0] / N1[0, 0]
    r3 = r2[1] - N2[0, 1] * r2[0] / N2[0, 0]
    al = r3 / N3
    u3 = (r2[0] - N2[0, 1] * al) / N2[0, 0]
    u2 = (r1[0] - N1[0, 1] * u3 - N1[0, 2] * al) / N1[0, 0]
    rows = [N[0, 0:], N[1, 1:], N[2, 2:], N[3, 3:], N1[0, 0:], N1[1, 1:], N1[2, 2:], N2[0, 0:], N2[1, 1:], [N3]]
    return rows, np.array([*r, *r1, *r2, r3, al, u3, u2])


for tab, _, comp in TABS:
    sub = lhs[lhs.table == tab].set_index('row')
    rows, _ = chain(A[comp].T @ A[comp], np.zeros(4))
    devs, big = [], []
    for i in range(10):
        pr = np.array([float(c) for c in sub.loc[str(i + 1), ['coef1', 'coef2', 'coef3', 'coef4']] if c != ''])
        for k, (a_, b_) in enumerate(zip(pr, rows[i])):
            devs.append(abs(a_ - b_))
            if abs(a_ - b_) > 0.0015:
                big.append(f'row {i + 1} coef{k + 5 - len(pr)} {a_:+.3f} vs {b_:+.4f}')
    verdict = 'ORIGINAL' if big else 'PASS'
    rep(f'Table {tab} printed LHS (20 coefficients) vs rebuilt', len(devs), max(devs), verdict, '; '.join(big))

# ================================================================== E. RHS
print('\n=== E. Right-hand sides and solutions of Tables III-VI rebuilt from Table II ===')
FAILCOLS = []
for label, used in (('E1 as printed', False), ('E2 with the 3 values used in 1920', True)):
    print(f'  --- {label} ---')
    for tab, plates, comp in TABS:
        q = QRA if comp == 'Dx' else QDE
        alld = []
        for plate in plates:
            v, o = table2(plate, comp, used)
            _, calc = chain(A[comp].T @ A[comp], A[comp].T @ (v - o))
            pr = np.array([nev[(tab, plate, k)] for k in q])
            d = pr - calc
            alld.append(d)
            flag = np.abs(d) > TOL
            if flag.any():
                print(f'    {tab:3s} {plate:6s}: ' + ', '.join(f'{q[i]} {d[i]:+.4f}' for i in np.where(flag)[0]))
                if not used:
                    FAILCOLS.append((tab, plate, comp))
        alld = np.array(alld)
        n_in = int(np.sum(np.abs(alld) <= TOL))
        rep(f'{label}: Table {tab} rhs1-rhs10 + 3 solutions', alld.size, float(np.abs(alld).max()),
            'PASS' if n_in == alld.size else 'ORIGINAL', f'{n_in}/{alld.size} within {TOL}; '
            f'rhs1-rhs4 max {np.abs(alld[:, :4]).max():.4f}')
v, o = table2('I', 'Dx')
print(f'  plate I R.A.: sum(Dx - offset) = {np.sum(v - o):.3f} with the printed star 2 value -0.789; Table III prints '
      f'rhs1 = {nev[("III", "I", "rhs1")]:+.3f}; with -0.798: {np.sum(v - o) - 0.009:.3f}')

# ================================================================== F. stepwise
print('\n=== F. Internal arithmetic of the printed tables ===')
STEPFAIL = []
for tab, plates, comp in TABS:
    sub = lhs[lhs.table == tab].set_index('row')
    g = lambda r, k: float(sub.loc[str(r), 'coef%d' % k])
    u3n, u2n = ('sol_a', 'sol_b') if comp == 'Dx' else ('sol_e', 'sol_d')
    devs = []
    for plate in plates:
        P = {k: nev[(tab, plate, k)] for k in (QRA if comp == 'Dx' else QDE)}
        calc = {'rhs5': P['rhs2'] - g(1, 2) / g(1, 1) * P['rhs1'], 'rhs6': P['rhs3'] - g(1, 3) / g(1, 1) * P['rhs1'],
                'rhs7': P['rhs4'] - g(1, 4) / g(1, 1) * P['rhs1'], 'rhs8': P['rhs6'] - g(5, 3) / g(5, 2) * P['rhs5'],
                'rhs9': P['rhs7'] - g(5, 4) / g(5, 2) * P['rhs5'], 'rhs10': P['rhs9'] - g(8, 4) / g(8, 3) * P['rhs8'],
                'sol_alpha': P['rhs10'] / g(10, 4),
                u3n: (P['rhs8'] - g(8, 4) * P['sol_alpha']) / g(8, 3)}
        calc[u2n] = (P['rhs5'] - g(5, 3) * P[u3n] - g(5, 4) * P['sol_alpha']) / g(5, 2)
        for k, c in calc.items():
            devs.append(abs(P[k] - c))
            if abs(P[k] - c) > TOL:
                STEPFAIL.append(f'{tab} {plate} {k}: printed {P[k]:+.3f}, from printed lines {c:+.4f}')
    rep(f'Table {tab}: rhs5-rhs10 and solutions from the printed lines above', len(devs), max(devs),
        'PASS' if max(devs) <= TOL else 'ORIGINAL')
for s in STEPFAIL:
    print('    ', s)

# ================================================================== G. digit search
print('\n=== G. Single-change search for every column failing E1 (rhs1-rhs4 must then agree within 0.0011) ===')


def variants(s):
    out = set()
    digs = [i for i, ch in enumerate(s) if ch.isdigit()]
    for i in digs:
        out |= {s[:i] + d + s[i + 1:] for d in '0123456789' if d != s[i]}
    for a, b in zip(digs, digs[1:]):
        t = list(s); t[a], t[b] = t[b], t[a]; out.add(''.join(t))
    out.add(('-' if s[0] == '+' else '+') + s[1:])
    out.discard(s)
    return out


for tab, plate, comp in FAILCOLS:
    q = QRA if comp == 'Dx' else QDE
    sub = t2s[t2s.plate == plate].set_index('star').loc[[str(s) for s in STARS]]
    vals = list(sub[comp + '_rev'])
    o_str = offs.loc[plate, comp + '_offset']
    pr = np.array([nev[(tab, plate, k)] for k in q])
    hits = []

    def score(vv, oo, prx):
        _, calc = chain(A[comp].T @ A[comp], A[comp].T @ (vv - oo))
        return np.abs(prx[:4] - calc[:4]).max(), np.abs(prx - calc).max()

    base = np.array([float(v_) for v_ in vals])
    for i, s in enumerate(vals):
        for vv in variants(s):
            b_ = base.copy(); b_[i] = float(vv)
            m4, m13 = score(b_, float(o_str), pr)
            if m4 <= 0.0011:
                hits.append(f'Table II star {STARS[i]} {comp} {s} -> {vv} (rhs1-4 {m4:.4f}, all 13 {m13:.4f})')
    for vv in variants(o_str):
        m4, m13 = score(base, float(vv), pr)
        if m4 <= 0.0011:
            hits.append(f'offset {o_str} -> {vv} (rhs1-4 {m4:.4f}, all 13 {m13:.4f})')
    for k in ('rhs1', 'rhs2', 'rhs3', 'rhs4'):
        for vv in variants(nes.set_index(['table', 'plate', 'quantity']).loc[(tab, plate, k), 'value']):
            prx = pr.copy(); prx[q.index(k)] = float(vv)
            m4, _ = score(base, float(o_str), prx)
            if m4 <= 0.0011:
                hits.append(f'printed {k} -> {vv} (rhs1-4 {m4:.4f})')
    print(f'  {tab} {plate}: ' + ('; '.join(hits[:4]) if hits else 'no single change in Table II, the offset or rhs1-4'))

# ================================================================== H. rounded products
print('\n=== H. rhs2-rhs4 as hand-computed: each product rounded to 3 decimals, then summed (E2 values) ===')
n_ex = n_rp = 0
big = []
for tab, plates, comp in TABS:
    cols = ('y', 'x', 'Ex') if comp == 'Dx' else ('x', 'y', 'Ey')
    for plate in plates:
        v, o = table2(plate, comp, used=True)
        u = [Decimal(f'{a:.3f}') - Decimal(f'{o:.3f}') for a in v]
        for k, col in zip(('rhs2', 'rhs3', 'rhs4'), cols):
            c = [Decimal(ccs.loc[str(s), col]) for s in STARS]
            ex = sum(ci * ui for ci, ui in zip(c, u)).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP)
            rp = sum((ci * ui).quantize(Decimal('0.001'), rounding=ROUND_HALF_UP) for ci, ui in zip(c, u))
            prv = Decimal(nes.set_index(['table', 'plate', 'quantity']).loc[(tab, plate, k), 'value'])
            n_ex += ex == prv; n_rp += rp == prv
            if abs(prv - rp) > Decimal('0.001'):
                big.append(f'{tab} {plate} {k}: printed {prv}, rounded products {rp} ({float(prv - rp):+.3f})')
print(f'  exact match of the 84 printed rhs2-rhs4: rounded exact sum {n_ex}/84, sum of rounded products {n_rp}/84')
rep('rhs2-rhs4 vs sum of rounded products: |dev| <= 0.001', 84, max(abs(float(b.split('(')[1][:6])) for b in big),
    'ORIGINAL', f'{84 - len(big)}/84 within 0.001; others: ' + '; '.join(big))

# ================================================================== I. cross-checks with pp.306-308
print('\n=== I. Cross-checks with Table VII (p.306), the p.307 tables and Table VIII (p.308) ===')
T7 = {'III': '+0.098 +0.126 +0.107 +0.148 +0.140 +0.073 +0.145', 'IV': '+0.042 +0.024 -0.015 +0.018 +0.020 +0.005 +0.008',
      'V': '+0.126 +0.139 +0.114 +0.111 +0.137 +0.139 +0.136', 'VI': '+0.044 +0.007 +0.021 +0.010 +0.040 +0.060 +0.036'}
T7MEAN = {'III': 0.120, 'IV': 0.015, 'V': 0.129, 'VI': 0.031}
d7 = []
for tab, plates, _ in TABS:
    pa = np.array([nev[(tab, p, 'sol_alpha')] for p in plates])
    d7 += list(np.abs(pa - np.array(T7[tab].split(), float)))
    d7.append(abs(round(pa.mean(), 3) - T7MEAN[tab]))
rep('alpha of Tables III-VI = Table VII columns, and the 4 printed means', 32, max(d7), 'PASS' if max(d7) < 1e-9 else 'FAIL')
m = {t: np.mean([nev[(t, p, 'sol_alpha')] for p in pl]) for t, pl, _ in TABS}
ra, de = m['III'] - m['IV'], m['V'] - m['VI']
print(f'  eclipse - comparison: R.A. {ra:+.4f} r (paper +0.105), Decl. {de:+.4f} r (paper +0.098); weights 1:2 -> '
      f'{(ra + 2 * de) / 3:+.4f} r = {(ra + 2 * de) / 3 * 6.25:+.3f}" at 50\' -> limb {(ra + 2 * de) / 3 * 6.25 * 50 / 15.8:.3f}" '
      f'(paper +0.100 r = +0.625", 1".98)')

CORR = {'ecl': (0.240, 0.168, 0.062, 0.062), 'scale': (0.423, 0.207, 0.096, 0.096),
        '14_2a': (0.409, 0.207, 0.091, 0.091), '14_2b': (0.409, 0.207, 0.091, 0.091), '15_1': (0.390, 0.207, 0.087, 0.087),
        '15_2': (0.370, 0.202, 0.087, 0.087), '17_1': (0.399, 0.207, 0.091, 0.091), '17_2': (0.337, 0.202, 0.077, 0.077),
        '18_2': (0.327, 0.202, 0.072, 0.072)}
P307 = {'I': (-0.025, -0.010, -0.237, -0.265, 0.000, -0.251), 'II': (-0.009, -0.024, -0.045, -0.107, 0.000, -0.076),
        'III': (0.006, 0.053, 0.014, 0.025, 0.000, 0.020), 'IV': (-0.056, 0.006, -0.027, -0.040, 0.000, -0.034),
        'V': (-0.055, -0.006, 0.008, 0.031, 0.000, 0.020), 'VII': (0.050, 0.017, 0.032, 0.050, 0.000, 0.041),
        'VIII': (-0.014, 0.031, 0.008, 0.006, 0.000, 0.007), '14_2a': (0.010, 0.004, 0.081, 0.033, 0.013, 0.057),
        '14_2b': (0.008, 0.030, -0.035, -0.049, 0.013, -0.042), '15_1': (-0.063, -0.085, 0.155, 0.086, -0.084, 0.120),
        '15_2': (-0.065, -0.075, 0.003, -0.035, -0.084, -0.016), '17_1': (-0.118, -0.077, -0.116, -0.174, -0.084, -0.145),
        '17_2': (-0.072, -0.109, 0.062, 0.029, -0.084, 0.046), '18_2': (-0.093, -0.087, -0.014, -0.074, -0.084, -0.044)}
d307, odd = [], []
for plate in ECL + CMP:
    ac, ec, bc, dc = CORR['ecl' if plate in ECL else plate]
    as_, es, bs, ds = CORR['scale']
    tra, tde = ('III', 'V') if plate in ECL else ('IV', 'VI')
    calc = (nev[(tra, plate, 'sol_a')] + ac - as_, nev[(tde, plate, 'sol_e')] + ec - es,
            nev[(tra, plate, 'sol_b')] + bc - bs, -(nev[(tde, plate, 'sol_d')] + dc - ds))
    for name, pv, cv in zip(('scale x', 'scale y', 'orient x', 'orient y'), P307[plate][:4], calc):
        d307.append(abs(pv - cv))
        if abs(pv - cv) > 0.0006:
            odd.append(f'{plate} {name}: p.307 {pv:+.3f}, from solutions {cv:+.3f}')
rep('p.307 scale/orientation = printed a, e, b, -d + refraction/aberration differences', 56, max(d307),
    'PASS' if len(odd) <= 1 else 'FAIL', '; '.join(odd) + ' (p.307 inconsistency, not part of T1)' if odd else '')

T8RA = {'I': '-0.18 -0.45 +0.08 -0.23 -0.14 +0.17 +0.75', 'II': '-0.51 -0.81 +0.11 -0.11 +0.23 +0.06 +1.03',
        'III': '-0.46 -0.38 -0.08 -0.19 -0.09 +0.14 +1.06', 'IV': '-0.07 -0.58 -0.11 -0.05 -0.11 -0.18 +1.09',
        'V': '-0.04 -0.60 -0.11 -0.02 -0.13 -0.11 +1.01', 'VII': '-0.72 -0.36 -0.16 -0.02 +0.13 +0.14 +0.98',
        'VIII': '-0.43 -0.62 -0.18 -0.01 -0.08 -0.01 +1.30', '14_2a': '-0.19 -0.42 -0.01 +0.14 +0.02 +0.17 +0.31',
        '14_2b': '-0.24 +0.16 +0.03 +0.09 -0.18 -0.06 +0.18', '15_1': '-0.23 -0.36 -0.01 +0.28 +0.26 +0.20 -0.16',
        '15_2': '-0.28 -0.32 +0.05 +0.10 +0.06 +0.18 +0.22', '17_1': '+0.11 -0.24 -0.04 -0.03 +0.13 +0.13 -0.04',
        '17_2': '-0.19 -0.33 +0.23 +0.21 +0.03 -0.02 +0.08', '18_2': '-0.02 -0.26 +0.08 -0.01 +0.14 +0.15 -0.06'}


def t8(plate, u):
    ac, ec, bc, dc = CORR['ecl' if plate in ECL else plate]
    as_, es, bs, ds = CORR['scale']
    S, O = P307[plate][4], P307[plate][5]
    r = u - ((S - (ac - as_)) * x + (O - (bc - bs)) * y)
    return 6.25 * (r - r.mean())


for label, used in (('as printed', False), ('with the values used in 1920', True)):
    res8 = []
    for plate in ECL + CMP:
        v, o = table2(plate, 'Dx', used)
        res8.append(np.array(T8RA[plate].split(), float) - t8(plate, v - o))
    res8 = np.array(res8)
    worst = np.unravel_index(np.argmax(np.abs(res8)), res8.shape)
    print(f'  Table VIII R.A. rebuilt from Table II {label}: rms {np.sqrt(np.mean(res8**2)):.3f}", worst '
          f'{(ECL + CMP)[worst[0]]} star {STARS[worst[1]]} {res8[worst]:+.3f}"')
    if used:
        clean = [(p, res8[i]) for i, p in enumerate(ECL + CMP) if p not in ('14_2b', '17_1')]
        cm = np.array([r_ for _, r_ in clean])
        rep('Table VIII R.A. (p.308) rebuilt from Table II, 84 values on 12 plates (14_2b and 17_1 excluded: '
            'Table VIII anomalies)', cm.size, float(np.abs(cm).max()), 'PASS' if np.abs(cm).max() <= 0.035 else 'FAIL',
            f'rms {np.sqrt(np.mean(cm**2)):.4f}" (Table VIII rounding 0.01")')
        i1 = (ECL + CMP).index('I'); v, o = table2('I', 'Dx')
        print(f'  plate I star 2: Table VIII +0.75"; rebuilt {t8("I", v - o)[-1]:+.3f}" with -0.789, '
              f'{t8("I", np.r_[v[:-1], -0.798] - o)[-1]:+.3f}" with -0.798')
        v, o = table2('14_2a', 'Dx')
        print(f'  plate 14_2a star 2: Table VIII +0.31"; rebuilt {t8("14_2a", v - o)[-1]:+.3f}" with -0.396, '
              f'{t8("14_2a", np.r_[v[:-1], -0.346] - o)[-1]:+.3f}" with -0.346')

# ================================================================== summary
print('\n=== SUMMARY (PASS = within tolerance; ORIGINAL = deviation traced to the 1920 print/arithmetic, kept as printed) ===')
for verdict, name, n, mx, extra in SUMMARY:
    print(f'  {verdict:8s} | {name} | n={n} | max|dev|={mx:.4f}' + (f' | {extra}' if extra else ''))
