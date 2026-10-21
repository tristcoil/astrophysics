"""Verify every transcribed Wallal CSV with its internal arithmetic, reading the files exactly as
pandas.read_csv(path, comment='#', skipinitialspace=True). Prints a report used in wallal1928_notes.md."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

D = Path(__file__).resolve().parents[1]
R_SUN = 0.2661                     # deg, Sun's radius adopted by Campbell & Trumpler
U28 = 27084 / 206264.806           # rad per 1928 coordinate unit (200 mm)
U23 = 4508 / 3600                  # deg per 1923 coordinate unit (100 mm)
P6 = ['M15_M7', 'M16_M11', 'M17_M9', 'N18_N11', 'N19_N13', 'N20_N9']
P4 = ['M15_M7', 'M17_M9', 'N18_N11', 'N20_N9']


def rd(name):
    return pd.read_csv(D / name, comment='#', skipinitialspace=True)


def hdr(t):
    print('\n' + '=' * 8 + ' ' + t)


def stats(devs, tol):
    devs = np.asarray(devs, float)
    return f'n={len(devs)} max={devs.max():.4f} n>{tol}={int((devs > tol + 1e-9).sum())}'


t4 = rd('wallal1928_table4_displacements.csv')
t5 = rd('wallal1928_table5_checkfield.csv')
t6 = rd('wallal1928_table6_radial.csv')
t7 = rd('wallal1928_table7_groups.csv')
t7a = rd('wallal1928_table7_allplates.csv')
t8 = rd('wallal1928_table8_checkfield_groups.csv')
t9 = rd('wallal1928_table9_limb.csv')
t10 = rd('wallal1928_table10_combined.csv')
t11 = rd('wallal1928_table11_groups.csv')
t12 = rd('wallal1928_table12_best15.csv')
ar = rd('wallal1928_areas.csv')
c23 = rd('wallal1923_table2_displacements.csv')
t1 = rd('wallal1928_table1_plates.csv')
t2 = rd('wallal1928_table2_refraction_aberration.csv')
t3 = rd('wallal1928_table3_plate_constants.csv')
for n, t in [('t1', t1), ('t2', t2), ('t3', t3), ('t4', t4), ('t5', t5), ('t6', t6), ('t7', t7), ('t7a', t7a), ('t8', t8),
             ('t9', t9), ('t10', t10), ('t11', t11), ('t12', t12), ('areas', ar), ('1923 II', c23)]:
    print(f'{n}: {t.shape[0]} rows x {t.shape[1]} cols')

# ------------------------------------------------------------------ Table 6
hdr('Table 6 internal')
assert t6['no'].is_unique and set(t6['no']) == set(range(1, 148)) - {66, 90}
assert (np.diff(t6['distance_deg']) >= 0).all()
wsum_dev, mean_dev = [], []
for _, r in t6.iterrows():
    w = np.array([r[f'wt_{p}'] for p in P6], float)
    d = np.array([r[f'Dr_{p}'] for p in P6], float)
    ok = ~np.isnan(d)
    assert (np.isnan(w) == np.isnan(d)).all()
    wsum_dev.append(abs(w[ok].sum() - r['wt_total']))
    mean_dev.append(abs((w[ok] * d[ok]).sum() / w[ok].sum() - r['Dr_mean']))
print('wt_total = sum(w):', stats(wsum_dev, 1e-9))
print('Dr_mean = sum(w Dr)/sum(w):', stats(mean_dev, 0.005))
print('values: Dr cells', int(t6[[f"Dr_{p}" for p in P6]].notna().sum().sum()), '; page counts', t6['page'].value_counts().to_dict())

# ------------------------------------------------------------------ Table 4
hdr('Table 4 internal + vs Table 6')
t4i = t4.set_index('no'); t6i = t6.set_index('no')
assert set(t4i.index) == set(t6i.index)
wdev, mx, my, nm = [], [], [], 0
rad_devs, rad_fail, wt_mismatch, rej_conflict = [], [], 0, 0
mean_rad = []
for no, r in t4i.iterrows():
    x, y = r['x'], r['y']; rr = math.hypot(x, y)
    sw = sx = sy = 0.0
    for p in P6:
        dx, dy, w, rej = r[f'Dx_{p}'], r[f'Dy_{p}'], r[f'wt_{p}'], r[f'rej_{p}']
        dr6, w6 = t6i.loc[no, f'Dr_{p}'], t6i.loc[no, f'wt_{p}']
        if np.isnan(dx):
            assert np.isnan(dr6); continue
        if rej == 1:
            if not np.isnan(dr6): rej_conflict += 1
            continue
        sw += w; sx += w * dx; sy += w * dy
        if w != w6: wt_mismatch += 1
        calc = (dx * x + dy * y) / rr
        dtan = abs(-dx * y + dy * x) / rr
        tol = 0.005 * (abs(x) + abs(y)) / rr + 0.0051 + dtan * 0.0005 * (abs(x) + abs(y)) / rr ** 2
        dev = abs(calc - dr6); rad_devs.append(dev)
        if dev > tol:
            rad_fail.append((int(no), p, round(calc, 4), dr6, round(dev, 4)))
    wdev.append(abs(sw - r['wt_mean']))
    mx.append(abs(sx / sw - r['Dx_mean'])); my.append(abs(sy / sw - r['Dy_mean']))
    calc = (r['Dx_mean'] * x + r['Dy_mean'] * y) / rr
    dtan = abs(-r['Dx_mean'] * y + r['Dy_mean'] * x) / rr
    tol = 0.005 * (abs(x) + abs(y)) / rr + 0.0051 + dtan * 0.0005 * (abs(x) + abs(y)) / rr ** 2
    mean_rad.append((abs(calc - t6i.loc[no, 'Dr_mean']), int(no), tol))
    if r['wt_mean'] != t6i.loc[no, 'wt_total']: nm += 1
print('wt_mean = sum(w):', stats(wdev, 1e-9), '; wt_mean vs Table 6 wt_total mismatches:', nm)
print('Dx_mean:', stats(mx, 0.005), '; >0.0101:', [int(t4i.index[i]) for i, v in enumerate(mx) if v > 0.0101])
print('Dy_mean:', stats(my, 0.005), '; >0.0101:', [int(t4i.index[i]) for i, v in enumerate(my) if v > 0.0101])
print('per-plate weights T4 vs T6 mismatches:', wt_mismatch, '; rejected-in-T4 but present-in-T6:', rej_conflict)
print('radial projection vs Table 6:', stats(rad_devs, 0.0101), '; beyond tolerance:', rad_fail)
print('mean radial vs Table 6 Dr_mean: max', round(max(m for m, _, _ in mean_rad), 4), '; beyond rounding tolerance:', [(n, round(m, 4)) for m, n, t in mean_rad if m > t])
print('rejected (parenthesised) measures:', int(sum((t4[f'rej_{p}'] == 1).sum() for p in P6)))

# distances vs x, y
dd = []
for no, r in t4i.iterrows():
    d = math.degrees(math.atan(math.hypot(r['x'], r['y']) * U28))
    dd.append(d - t6i.loc[no, 'distance_deg'])
dd = np.array(dd)
print(f'Table 6 distance vs gnomonic(x,y): mean {dd.mean():+.4f} rms {np.sqrt((dd**2).mean()):.4f} max {np.abs(dd).max():.4f}')

# RA/Dec vs x, y
try:
    from astropy.coordinates import SkyCoord, FK4
    import astropy.units as u
    from scipy.optimize import least_squares
    ra = np.asarray(15 * (t4['ra1900_h'] + t4['ra1900_min'] / 60), float)
    dec = np.asarray(np.where(t4['dec1900_sign'] == '-', -1, 1) * (t4['dec1900_deg'] + t4['dec1900_arcmin'] / 60), float)
    c = SkyCoord(ra * u.deg, dec * u.deg, frame=FK4(equinox='B1900', obstime='B1900')).transform_to(
        FK4(equinox='B1922', obstime='B1922'))
    a, dl = c.ra.rad, c.dec.rad
    X, Y = t4['x'].values, t4['y'].values
    mask = t4['no'].values != 129

    def model(p):
        a0, d0, s, th = p
        ca = np.cos(a - a0); Dn = np.sin(dl) * np.sin(d0) + np.cos(dl) * np.cos(d0) * ca
        xi = np.cos(dl) * np.sin(a - a0) / Dn; eta = (np.sin(dl) * np.cos(d0) - np.cos(dl) * np.sin(d0) * ca) / Dn
        return ((xi * np.cos(th) - eta * np.sin(th)) / (U28 * s), (xi * np.sin(th) + eta * np.cos(th)) / (U28 * s))

    sol = least_squares(lambda p: np.concatenate([(model(p)[0] - X)[mask], (model(p)[1] - Y)[mask]]),
                        [math.radians(177.65), math.radians(1.02), 1.0, 0.0])
    mxp, myp = model(sol.x)
    rx, ry = mxp - X, myp - Y
    print(f'RA/Dec(1900->1922 FK4) vs x,y, 4-parameter gnomonic fit without star 129: rms x {np.sqrt((rx[mask]**2).mean()):.4f} '
          f'y {np.sqrt((ry[mask]**2).mean()):.4f} units, max |dx| {np.abs(rx[mask]).max():.4f} |dy| {np.abs(ry[mask]).max():.4f}; '
          f'fitted centre RA {math.degrees(sol.x[0])/15:.4f} h Dec {math.degrees(sol.x[1]):.4f} deg')
    i129 = int(np.where(~mask)[0][0])
    print(f'star 129 residual: dx {rx[i129]:+.4f} dy {ry[i129]:+.4f} units ({ry[i129]*7.523*60:+.0f} arcmin in y)')
    # BD+4 2583 (1855: 12h04m16.4s +4d51.5m, BD catalogue via VizieR I/122) as the measured star?
    cb = SkyCoord('12h04m16.4s', '+4d51.5m', frame=FK4(equinox='B1855', obstime='B1855')).transform_to(
        FK4(equinox='B1922', obstime='B1922'))
    a_s, d_s = a.copy(), dl.copy(); a[i129], dl[i129] = cb.ra.rad, cb.dec.rad
    mxp, myp = model(sol.x)
    print(f'  with BD+4 2583 instead: dx {mxp[i129]-X[i129]:+.4f} dy {myp[i129]-Y[i129]:+.4f} units')
    a[:], dl[:] = a_s, d_s
except Exception as e:  # pragma: no cover
    print('RA/Dec check skipped:', e)

# ------------------------------------------------------------------ Table 5
hdr('Table 5 internal')
assert len(t5) == 75 and list(t5['no']) == [f'{i}c' for i in range(1, 76)]
wdev, mx, my = [], [], []
for _, r in t5.iterrows():
    sw = sx = sy = 0.0
    for p in P4:
        if np.isnan(r[f'vx_{p}']) or r[f'rej_{p}'] == 1: continue
        w = r[f'wt_{p}']; sw += w; sx += w * r[f'vx_{p}']; sy += w * r[f'vy_{p}']
    wdev.append(abs(sw - r['wt_mean'])); mx.append(abs(sx / sw - r['vx_mean'])); my.append(abs(sy / sw - r['vy_mean']))
print('wt_mean = sum(w):', stats(wdev, 1e-9), '; vx_mean:', stats(mx, 0.005), '; vy_mean:', stats(my, 0.005))
print('rejected (parenthesised) residuals:', int(sum((t5[f'rej_{p}'] == 1).sum() for p in P4)))

# ------------------------------------------------------------------ Table 7 from Table 6
hdr('Table 7 vs Table 6')
sizes = [5, 6, 9, 25, 20, 19, 19, 20, 22]
groups, i = [], 0
for s in sizes:
    groups.append(t6.iloc[i:i + s]); i += s
series = {'M15-M7': ['M15_M7'], 'M16-M11': ['M16_M11'], 'M17-M9': ['M17_M9'], 'N18-N11': ['N18_N11'],
          'N19-N13': ['N19_N13'], 'N20-N9': ['N20_N9'], 'Mean_M': ['M15_M7', 'M16_M11', 'M17_M9'],
          'Mean_N': ['N18_N11', 'N19_N13', 'N20_N9']}
wd, rdv, ddv = [], [], []
for _, r in t7.iterrows():
    g = groups[int(r['group']) - 1]
    sw = swd = swr = 0.0
    for p in series[r['series']]:
        ok = g[f'Dr_{p}'].notna()
        sw += g.loc[ok, f'wt_{p}'].sum(); swr += (g.loc[ok, f'wt_{p}'] * g.loc[ok, f'Dr_{p}']).sum()
        swd += (g.loc[ok, f'wt_{p}'] * g.loc[ok, 'distance_deg']).sum()
    wd.append(abs(sw - r['weight'])); rdv.append(abs(swr / sw - r['Dr'])); ddv.append(abs(swd / sw - r['dist_deg']))
print('weights:', stats(wd, 1e-6), '; Dr:', stats(rdv, 0.005), '; dist (weighted):', stats(ddv, 0.005))
hdr('Table 7 all plates')
for _, r in t7a.iterrows():
    g = groups[int(r['group']) - 1]
    W = g['wt_total'].sum(); Dr = (g['wt_total'] * g['Dr_mean']).sum() / W; dist = (g['wt_total'] * g['distance_deg']).sum() / W
    ein = 1.75 * R_SUN / r['dist_deg'] - 0.009 * r['dist_deg']
    print(f"group {int(r['group'])}: n {len(g)}/{int(r['n_stars'])} W {W:.2f}/{r['weight']} Dr {Dr:+.4f}/{r['Dr']:+.2f} "
          f"dist {dist:.3f}/{r['dist_deg']} einstein {ein:+.4f}/{r['einstein_rel']:+.2f} "
          f"O-T {r['Dr']-r['einstein_rel']:+.2f}/{r['obs_minus_theory']:+.2f} pe(0.33/sqrtW) {0.33/math.sqrt(r['weight']):.3f}/{r['pe']}")

# ------------------------------------------------------------------ Table 8 from Table 5
hdr('Table 8 vs Table 5')
rows8 = []
for _, r in t5.iterrows():
    rr = math.hypot(r['x'], r['y'])
    rows8.append((math.degrees(math.atan(rr * U28)), r['wt_mean'], (r['vx_mean'] * r['x'] + r['vy_mean'] * r['y']) / rr))
rows8.sort()
i = 0
for _, r in t8.iterrows():
    acc, j = 0.0, i
    while j < len(rows8) and acc + rows8[j][1] <= r['weight'] + 1e-6:
        acc += rows8[j][1]; j += 1
    g = rows8[i:j]; W = sum(x[1] for x in g)
    print(f"group {int(r['group'])}: n {len(g)} W {W:.2f}/{r['weight']} Dr {sum(x[1]*x[2] for x in g)/W:+.4f}/{r['Dr']:+.2f} "
          f"dist {sum(x[0]*x[1] for x in g)/W:.3f}/{r['dist_deg']}")
    i = j
print('stars used', i, 'of', len(rows8))

# ------------------------------------------------------------------ Table 9
hdr('Table 9')
pl = t9[t9['weight'].notna()]
for lab, sel in [('Mean of obj. M', pl['plates'].str.startswith('M')), ('Mean of obj. N', pl['plates'].str.startswith('N')),
                 ('Mean of all plates', pl['plates'].notna())]:
    s = pl[sel]; W = s['weight'].sum()
    print(f"{lab}: weighted beta {(s['weight']*s['beta']).sum()/W:.4f} printed {t9.loc[t9['plates']==lab,'beta'].iloc[0]}"
          + (f"; alpha {(s['weight']*s['alpha']).sum()/W:+.4f} printed {t9.loc[t9['plates']==lab,'alpha'].iloc[0]}" if lab.startswith('Mean of all') else ''))
W = pl['weight'].sum()
print(f'p.e.(mean) = 0.33/sqrt(sum w = {W:.2f}) = {0.33/math.sqrt(W):.3f} (printed .15); '
      f"p.e. unit weight from beta scatter = {0.6745*math.sqrt((pl['weight']*(pl['beta']-1.82)**2).sum()/5):.3f} (printed .37 on p.151)")
for _, r in pl.iterrows():
    p = r['plates'].replace('-', '_')
    ok = t6[f'Dr_{p}'].notna()
    rr = t6.loc[ok, 'distance_deg'].values; d = t6.loc[ok, f'Dr_{p}'].values; w = t6.loc[ok, f'wt_{p}'].values
    A = np.column_stack([-rr, R_SUN / rr]); N = A.T @ (w[:, None] * A)
    a, b = np.linalg.solve(N, A.T @ (w * d)); wb = 1 / np.linalg.inv(N)[1, 1]
    print(f"refit {r['plates']}: n {ok.sum()} beta {b:.3f}/{r['beta']} weight {wb:.3f}/{r['weight']} alpha {a:+.4f}/{r['alpha']}")

# ------------------------------------------------------------------ Table 10
hdr('Table 10 vs Table 6 and 1923 Table II')
c23i = c23.set_index('no')
ow, oo, od, oe = [], [], [], []
for _, r in t10.iterrows():
    parts, dist = [], None
    if not np.isnan(r['n15']):
        g = c23i.loc[int(r['n15'])]; d15 = g['dist'] * U23
        parts.append((g['weight'], g['D1'] + 0.060 * d15)); dist = d15
    if not np.isnan(r['n5']):
        g = t6i.loc[int(r['n5'])]; parts.append((g['wt_total'] / 3, g['Dr_mean'])); dist = g['distance_deg']
    W = sum(p[0] for p in parts)
    ow.append(abs(W - r['weight'])); oo.append(abs(sum(p[0] * p[1] for p in parts) / W - r['Dr']))
    od.append(abs(dist - r['dist_deg']))
    oe.append(abs(1.75 * R_SUN / r['dist_deg'] - 0.009 * r['dist_deg'] - r['einstein_rel']))
print('weight:', stats(ow, 0.005), '; Dr:', stats(oo, 0.0101), '; dist:', stats(od, 0.005), '; einstein:', stats(oe, 0.005))
in5 = set(t6.loc[t6['distance_deg'] <= 3.45, 'no']); got5 = set(t10['n5'].dropna().astype(int))
print('5-ft stars <= 3.45 deg in Table 6:', len(in5), '; listed in Table 10:', len(got5), '; symmetric diff:', sorted(in5 ^ got5))
print('15-ft stars listed:', int(t10['n15'].notna().sum()), 'of', len(c23), '; missing:', sorted(set(c23['no']) - set(t10['n15'].dropna().astype(int))))

# ------------------------------------------------------------------ Table 11
hdr('Table 11 vs Tables 10 and 7')
sz = [10, 17, 14, 12, 14, 18, 18]
i = 0
for k, s in enumerate(sz):
    g = t10.iloc[i:i + s]; i += s; r = t11.iloc[k]
    W = g['weight'].sum()
    print(f"group {k+1}: n {s}/{int(r['n_stars'])} W {W:.2f}/{r['weight']} Dr {(g['weight']*g['Dr']).sum()/W:+.4f}/{r['Dr']:+.2f} "
          f"dist {(g['weight']*g['dist_deg']).sum()/W:.3f}/{r['dist_deg']} pe(0.20/sqrtW printed) {0.2/math.sqrt(r['weight']):.3f}/{r['pe']}")
for k in range(4):
    r = t11.iloc[7 + k]; a = t7a.iloc[5 + k]
    print(f"group {8+k}: from Table 7 group {6+k}: W {a['weight']/3:.2f}/{r['weight']} Dr {a['Dr']:+.2f}/{r['Dr']:+.2f} "
          f"n {int(a['n_stars'])}/{int(r['n_stars'])} dist {a['dist_deg']}/{r['dist_deg']}")
ce, cc, oe2, oc2 = [], [], [], []
for _, r in t11.iterrows():
    x = r['dist_deg']
    ce.append(abs(1.75 * R_SUN / x - 0.009 * x - r['einstein_rel']))
    cc.append(abs(0.55 * (1 - math.sqrt(math.sin(math.radians(x / 2)))) - 0.0735 * x - r['courvoisier_rel']))
    oe2.append(abs(r['Dr'] - r['einstein_rel'] - r['obs_minus_einstein'])); oc2.append(abs(r['Dr'] - r['courvoisier_rel'] - r['obs_minus_courvoisier']))
print('einstein_rel:', stats(ce, 0.005), '; courvoisier_rel:', stats(cc, 0.005), '; O-E:', stats(oe2, 1e-6), '; O-C:', stats(oc2, 1e-6))

# ------------------------------------------------------------------ Table 12
hdr('Table 12')
o1, o2, o3 = [], [], []
t4n = t4.set_index('no')
for _, r in t12.iterrows():
    o1.append(max(abs(r['Dx_obs'] - r['Dx_theory'] - r['Dx_omt']), abs(r['Dy_obs'] - r['Dy_theory'] - r['Dy_omt'])))
    g = c23i.loc[r['n15']]; X, Y = g['X'], g['Y']; rr = math.hypot(X, Y); dist = t6i.loc[r['n5'], 'distance_deg']
    th = 1.75 * R_SUN / dist - 0.009 * dist
    o2.append(max(abs(th * X / rr - r['Dx_theory']), abs(th * Y / rr - r['Dy_theory'])))
    w15 = g['weight']; cx, cy = g['Ax'] + 0.060 * X * U23, g['Ay'] + 0.060 * Y * U23
    m = t4n.loc[r['n5']]; w5 = m['wt_mean'] / 3
    o3.append(max(abs((w15 * cx + w5 * m['Dx_mean']) / (w15 + w5) - r['Dx_obs']), abs((w15 * cy + w5 * m['Dy_mean']) / (w15 + w5) - r['Dy_obs'])))
print('omt = obs - theory:', stats(o1, 1e-6), '; theory (1923 X,Y; Table 6 r):', stats(o2, 0.005), '; obs from both instruments:', stats(o3, 0.0101))

# ------------------------------------------------------------------ areas
hdr('8 areas vs Table 4')
area_of = {}
for _, r in t4.iterrows():
    inner = abs(r['x']) < 0.5 and abs(r['y']) < 0.5
    q = {(True, True): 0, (False, True): 1, (False, False): 2, (True, False): 3}[(r['x'] < 0, r['y'] > 0)]
    area_of.setdefault((['V', 'VI', 'VII', 'VIII'] if inner else ['I', 'II', 'III', 'IV'])[q], []).append(r)
for _, r in ar[ar['n_stars'].notna()].iterrows():
    g = pd.DataFrame(area_of[r['area']]); W = g['wt_mean'].sum()
    print(f"area {r['area']}: n {len(g)}/{int(r['n_stars'])} Dx {(g['wt_mean']*g['Dx_mean']).sum()/W:+.3f}/{r['Dx_obs']:+.2f} "
          f"Dy {(g['wt_mean']*g['Dy_mean']).sum()/W:+.3f}/{r['Dy_obs']:+.2f}  corr-obs ({r['Dx_corr']-r['Dx_obs']:+.2f},{r['Dy_corr']-r['Dy_obs']:+.2f})")

# ------------------------------------------------------------------ 1923 Table II
hdr('1923 Table II')
dd, d1, th = [], [], []
for _, r in c23.iterrows():
    X, Y = r['X'], r['Y']; rr = math.hypot(X, Y)
    dd.append(abs(rr - r['dist']))
    c1 = (r['Ax'] * X + r['Ay'] * Y) / rr; t = abs(-r['Ax'] * Y + r['Ay'] * X) / rr
    tol = 0.0051 + 0.005 * (abs(X) + abs(Y)) / rr + t * 0.005 * math.sqrt(2) / rr
    d1.append((abs(c1 - r['D1']), tol, int(r['no'])))
    th.append(abs(1.75 * R_SUN / (r['dist'] * U23) - 0.079 * r['dist'] - r['theory']))
print('dist = sqrt(X^2+Y^2):', stats(dd, 0.0101), '; D1 = radial(Ax,Ay):', stats([x[0] for x in d1], 0.0101),
      'beyond rounding tol:', [x[2] for x in d1 if x[0] > x[1]], '; theory:', stats(th, 0.005))
t3_23 = [(8, 9.09, 0.64, 0.64, 0.69, 0.70), (11, 19.42, 1.06, 0.35, 0.46, 0.37), (10, 20.15, 1.40, 0.30, 0.39, 0.24),
         (8, 22.41, 1.66, 0.16, 0.22, 0.17), (9, 21.10, 1.90, 0.17, 0.21, 0.13), (8, 24.67, 2.00, 0.15, 0.17, 0.11),
         (11, 21.32, 2.22, 0.08, 0.08, 0.08), (13, 21.37, 2.55, -0.09, -0.14, 0.02), (14, 22.78, 2.97, -0.04, -0.08, -0.03)]
srt = c23.sort_values(['dist', 'no'], kind='stable').reset_index(drop=True)
i = 0
for k, (n, W, dist, d1m, d2m, thm) in enumerate(t3_23):
    g = srt.iloc[i:i + n]; i += n; w = g['weight']
    print(f"Table III group {k+1}: n {n} W {w.sum():.2f}/{W} dist {(w*g['dist']).sum()/w.sum()*U23:.3f}/{dist} "
          f"D1 {(w*g['D1']).sum()/w.sum():+.3f}/{d1m:+.2f} D2 {(w*g['D2']).sum()/w.sum():+.3f}/{d2m:+.2f} "
          f"theory {(w*g['theory']).sum()/w.sum():+.3f}/{thm:+.2f}")

# ------------------------------------------------------------------ Table 1
hdr('Table 1 simultaneous M/N exposures')
pairs = [('M15', 'N18'), ('M16', 'N19'), ('M17', 'N20'), ('M7', 'N9'), ('M9', 'N11'), ('M11', 'N13')]
t1['plate_f'] = t1['plate']
for a, b in pairs:
    A = t1[t1['plate'] == a].reset_index(drop=True); B = t1[t1['plate'] == b].reset_index(drop=True)
    for k in range(len(A)):
        ra_, rb_ = A.iloc[k], B.iloc[k]
        diffs = [c for c in ['date_1922', 'lmt_h', 'lmt_min', 'exposure_s', 'star_field', 'ha_w_h', 'ha_w_min', 'temperature_F', 'barometer_in']
                 if str(ra_[c]) != str(rb_[c])]
        if diffs:
            print(f'{a} vs {b} line {k+1}:', {c: (ra_[c], rb_[c]) for c in diffs})
print('Table 1 rows', len(t1))
