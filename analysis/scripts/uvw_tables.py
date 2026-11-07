#!/usr/bin/env python3
"""uvw_tables.py - comparisons for the capped-ECP test (U), substituted boundary cages (V) and mid-grain sites (W).

python3 uvw_tables.py tables/cf_all_runs.csv tables/uvw_tests.csv
Each row: test, ion, site, run, reference run, A20r2/split/S2/S4/S6/E1/tilt of both, and % changes."""
import sys, csv
src, out = sys.argv[1:3]
R = {r['run']: r for r in csv.DictReader(open(src)) if r['block'] == 'CASSCF'}
P = {k.split('_')[0]: k for k in R}
f = lambda k, c: float(R[k][c]) if R[k].get(c, '') not in ('', 'nan') else float('nan')
SITES = ['bulk4f5', 'GBA13', 'bulk4g2', 'GBB9']
pairs = []
for ion, base in (('Nd', 1), ('Dy', 9)):
    for i, s in enumerate(SITES):
        pairs.append(('capped ECP vs bare pc (def2-ECP level)', ion, s, P['U%02d' % (base + 2 * i + 1)], P['U%02d' % (base + 2 * i)]))
DKH = dict(Nd=dict(bulk4f5='T03', GBA13='R12', bulk4g2='R10', GBB9='T09'), Dy=dict(bulk4f5='T06', GBA13='R24', bulk4g2='R22', GBB9='T12'))
for ion, base in (('Nd', 1), ('Dy', 9)):
    for i, s in enumerate(SITES):
        pairs.append(('def2-ECP pc vs DKH production', ion, s, P['U%02d' % (base + 2 * i)], P[DKH[ion][s]]))
for v, ref, ion, s in (('V02', 'R12', 'Nd', 'GBA13 Zr for nearest RE'), ('V03', 'T09', 'Nd', 'GBB9 B for nearest Fe'), ('V04', 'T09', 'Nd', 'GBB9 Zr for nearest RE'),
                       ('V05', 'R24', 'Dy', 'GBA13 B for nearest Fe'), ('V06', 'R24', 'Dy', 'GBA13 Zr for nearest RE'), ('V07', 'T12', 'Dy', 'GBB9 B for nearest Fe'),
                       ('V08', 'T12', 'Dy', 'GBB9 Zr for nearest RE')):
    pairs.append(('substitution vs unsubstituted GB site', ion, s, P[v], P[ref]))
for w, ref, ion, s in (('W01', 'R10', 'Nd', 'GB class C (2) vs bulk 4g'), ('W02', 'T03', 'Nd', 'GB class D (6) vs bulk 4f'),
                       ('W03', 'R22', 'Dy', 'GB class C (2) vs bulk 4g'), ('W04', 'T06', 'Dy', 'GB class D (6) vs bulk 4f')):
    pairs.append(('mid-grain site vs parent bulk', ion, s, P[w], P[ref]))
cols = ['A20r2', 'split', 'S2', 'S4', 'S6', 'E1', 'g_angle', 'gz']
with open(out, 'w', newline='') as fh:
    wr = csv.writer(fh)
    wr.writerow(['test', 'ion', 'site', 'run', 'reference'] + [c + x for c in cols for x in ('_run', '_ref')] + ['d_' + c + '_pct' for c in cols[:6]])
    for t, ion, s, a, b in pairs:
        va = [f(a, c) for c in cols]; vb = [f(b, c) for c in cols]
        row = [t, ion, s, a, b] + ['%.3f' % x for pair in zip(va, vb) for x in pair]
        row += ['%.2f' % (100 * (x / y - 1)) if y == y and x == x and y != 0 else '' for x, y in zip(va[:6], vb[:6])]
        wr.writerow(row)
print('wrote', out, len(pairs), 'rows')
