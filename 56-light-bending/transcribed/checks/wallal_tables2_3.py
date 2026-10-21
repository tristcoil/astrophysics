"""Tables 2 and 3 of the 1928 paper, read from the transcribed CSVs: internal arithmetic.
 A) Check F. corr  = det - dR_chk - dA_chk - I   (I: d and e' get d_i; e and f' get e_i)
 B) Eclipse F. adopt = Check F. corr + dR_ecl + dA_ecl   (coefficients d..k, d'..k')
 C) Eclipse F. corr  = adopt - dR_ecl - dA_ecl            (b, c, b', c')
 D) M16-M11, N19-N13: corr = adopt - dR_ecl - dA_ecl      (b..k, b'..k')
 dR = R(Wallal plate) - R(Tahiti plate), dA likewise (Table 2).
 Variants: 'printed' uses Table 2 as printed (no g_a, h'_a rows -> 0);
           'model' sets the unprinted g_a = i_a and h'_a = k'_a;
           'model+ia' additionally replaces the printed i_a = -.01 of "M7 and N9, Check F." by +.01 (= its k'_a)."""
from pathlib import Path
import sys
import pandas as pd

D = str(Path(__file__).resolve().parents[1]) + '/'
t2 = pd.read_csv(D + 'wallal1928_table2_refraction_aberration.csv', comment='#', skipinitialspace=True)
t3 = pd.read_csv(D + 'wallal1928_table3_plate_constants.csv', comment='#', skipinitialspace=True)


def key(plates, field):
    return f'{plates}|{field}'


R, A, I = {}, {}, {}
for _, r in t2.iterrows():
    k = key(r['plates'], r['field'])
    tgt = {'refraction': R, 'aberration': A, 'inclination': I}[r['block']]
    tgt.setdefault(r['coef'], {})[k if r['block'] != 'inclination' else r['plates']] = float(r['value'])
COEF = ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', "b'", "c'", "d'", "e'", "f'", "g'", "h'", "i'", "k'"]
W_ECL_A = key('M15, M16, M17, N18, N19, N20', 'eclipse')
MAP = {  # pair: (R wallal chk, R tahiti chk, R wallal ecl, R tahiti ecl, A wallal chk, A tahiti chk, A tahiti ecl)
    'M15-M7': (key('M15 and N18', 'check'), key('M7 and N9', 'check'), key('M15 and N18', 'eclipse'), key('M7 and N9', 'eclipse'),
               key('M15, N18', 'check'), key('M7 and N9', 'check'), key('M7 and N9', 'eclipse')),
    'M17-M9': (key('M17 and N20', 'check'), key('M9 and N11', 'check'), key('M17 and N20', 'eclipse'), key('M9 and N11', 'eclipse'),
               key('M17, N20', 'check'), key('M9 and N11', 'check'), key('M9 and N11', 'eclipse')),
    'N18-N11': (key('M15 and N18', 'check'), key('M9 and N11', 'check'), key('M15 and N18', 'eclipse'), key('M9 and N11', 'eclipse'),
                key('M15, N18', 'check'), key('M9 and N11', 'check'), key('M9 and N11', 'eclipse')),
    'N20-N9': (key('M17 and N20', 'check'), key('M7 and N9', 'check'), key('M17 and N20', 'eclipse'), key('M7 and N9', 'eclipse'),
               key('M17, N20', 'check'), key('M7 and N9', 'check'), key('M7 and N9', 'eclipse')),
    'M16-M11': (None, None, key('M16, N19', 'eclipse'), key('M11, N13', 'eclipse'), None, None, key('M11, N13', 'eclipse')),
    'N19-N13': (None, None, key('M16, N19', 'eclipse'), key('M11, N13', 'eclipse'), None, None, key('M11, N13', 'eclipse')),
}


def run(variant):
    Av = {c: dict(v) for c, v in A.items()}
    if variant in ('model', 'model+ia'):
        Av['g'] = dict(Av['i']); Av["h'"] = dict(Av["k'"])
    if variant == 'model+ia':
        Av['i'][key('M7 and N9', 'check')] = 0.01; Av['g'][key('M7 and N9', 'check')] = 0.01

    def a(c, k):
        return Av.get(c, {}).get(k, 0.0)

    devs, fails = [], []
    for pair, (rwc, rtc, rwe, rte, awc, atc, ate) in MAP.items():
        sub = t3[t3['pair'] == pair].set_index('coef')
        for c in COEF:
            row = sub.loc[c]
            dE = R[c][rwe] - R[c][rte] + a(c, W_ECL_A) - a(c, ate)
            if rwc:
                inc = I["d_i=e'_i"][pair] if c in ('d', "e'") else (I["e_i=f'_i"][pair] if c in ('e', "f'") else 0.0)
                exp = row['check_det'] - (R[c][rwc] - R[c][rtc]) - (a(c, awc) - a(c, atc)) - inc
                devs.append(abs(exp - row['check_corr']))
                if devs[-1] > 0.0151: fails.append(f'A {pair} {c}: det {row["check_det"]:+.2f} -> corr {exp:+.3f}, printed {row["check_corr"]:+.2f}')
                if c in ('b', 'c', "b'", "c'"):
                    exp = row['eclipse_adopt'] - dE; devs.append(abs(exp - row['eclipse_corr']))
                    if devs[-1] > 0.0151: fails.append(f'C {pair} {c}: adopt -> corr {exp:+.3f}, printed {row["eclipse_corr"]:+.2f}')
                else:
                    exp = row['check_corr'] + dE; devs.append(abs(exp - row['eclipse_adopt']))
                    if devs[-1] > 0.0151: fails.append(f'B {pair} {c}: corr -> adopt {exp:+.3f}, printed {row["eclipse_adopt"]:+.2f}')
            else:
                exp = row['eclipse_adopt'] - dE; devs.append(abs(exp - row['eclipse_corr']))
                if devs[-1] > 0.0151: fails.append(f'D {pair} {c}: adopt -> corr {exp:+.3f}, printed {row["eclipse_corr"]:+.2f}')
    print(f'{variant}: {len(devs)} relations, max dev {max(devs):.3f}, n>0.005: {sum(d > 0.0051 for d in devs)}, '
          f'n>0.015: {len(fails)}')
    for f in fails:
        print('   ', f)


for v in ('printed', 'model', 'model+ia'):
    run(v)
# aberration symmetries in every printed column
for k in A['b']:
    bad = []
    if abs(A['b'][k] - A["c'"][k]) > 1e-9: bad.append("b_a != c'_a")
    if abs(A['c'][k] + A["b'"][k]) > 1e-9: bad.append("c_a != -b'_a")
    if abs(A['f'][k] + A['d'][k]) > 1e-9: bad.append('f_a != -d_a')
    if abs(A['i'][k] - A["k'"][k]) > 1e-9: bad.append(f"i_a {A['i'][k]:+.2f} != k'_a {A[chr(107)+chr(39)][k]:+.2f}")
    print('aberration column', k, ':', 'b=c\', c=-b\', f=-d, i=k\' all hold' if not bad else bad)
for p in I['p']:
    print(f"inclination {p}: d_i/p = {I[chr(100)+'_i=e'+chr(39)+'_i'][p]/I['p'][p]:.4f}, e_i/q = {I['e_i=f'+chr(39)+'_i'][p]/I['q'][p]:.4f}")
# scale check-eclipse table of p.138 from Table 3 corr columns
for pair in ['M15-M7', 'M17-M9', 'N18-N11', 'N20-N9']:
    s = t3[t3['pair'] == pair].set_index('coef')
    print(f"{pair}: b check-ecl {s.loc['b','check_corr']-s.loc['b','eclipse_corr']:+.2f}, c' check-ecl {s.loc[chr(99)+chr(39),'check_corr']-s.loc[chr(99)+chr(39),'eclipse_corr']:+.2f}")
