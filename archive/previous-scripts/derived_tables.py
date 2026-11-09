#!/usr/bin/env python3
"""derived_tables.py - matched-site comparison (Table 2) and sensitivity tests (Table 3, Table S5) from cf_all_runs.csv.

python3 derived_tables.py tables/cf_all_runs.csv tables/
writes tables/matched_pairs_R25.csv and tables/sensitivity_tests.csv
Matched pairs (R = 25 A): GB class A (atom 13) vs bulk 4f (atom 5); GB class B (atom 9) vs bulk 4g (atom 2).
"""
import sys, os, csv

src, outdir = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(src)))
CAS = {r['run']: r for r in rows if r['block'] == 'CASSCF'}
NEV = {r['run']: r for r in rows if r['block'] == 'NEVPT2'}
IONS = ['Ce', 'Pr', 'Nd', 'Gd', 'Tb', 'Dy']
pick = lambda prefix: next((k for k in CAS if k.startswith(prefix)), None)
num = lambda x: float(x) if x not in ('', None, 'nan') else None
pct = lambda a, b: 100 * (b / a - 1) if a not in (None, 0) and b is not None else ''

# ------------------------------------------------------------------ matched pairs
KEYS = ['split', 'A20r2', 'S2', 'S4', 'S6', 'f_ax', 'E1', 'tun', 'Jspan', 'gz', 'g_angle', 'ax2_angle', 'orb_d1']
PCT = ['split', 'A20r2', 'S2', 'S4', 'S6', 'E1']
R25 = {ion: (pick('R%02d_%s_bulk_R25' % (4 * i + 2, ion)), pick('R%02d_%s_GB_R25' % (4 * i + 4, ion))) for i, ion in enumerate(IONS)}
out = []
for i, ion in enumerate(IONS):
    for pair, bulk, gb in (('A (4f)', pick('T%02d_%s_bulk4f' % (i + 1, ion)), R25[ion][1]),
                           ('B (4g)', R25[ion][0], pick('T%02d_%s_GBclassB' % (i + 7, ion)))):
        r = dict(ion=ion, pair=pair)
        ok = bulk in CAS and gb in CAS and CAS[gb]['normal'] == 'True' and not (ion == 'Ce' and pair.startswith('A'))
        for k in KEYS:
            r[k + '_b'] = CAS[bulk].get(k, '') if ok else ''; r[k + '_g'] = CAS[gb].get(k, '') if ok else ''
        for k in PCT:
            r[k + '_pct'] = pct(num(r[k + '_b']), num(r[k + '_g'])) if ok else ''
        out.append(r)
with open(os.path.join(outdir, 'matched_pairs_R25.csv'), 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(out[0])); wr.writeheader(); wr.writerows(out)

# ------------------------------------------------------------------ sensitivity tests
T = []
def add(test, ion, site, ref, var, refrow, varrow):
    r = dict(test=test, ion=ion, site=site, ref=ref, var=var)
    for k in ('split', 'A20r2', 'E1'):
        a, b = num(refrow.get(k)), num(varrow.get(k)); r[k + '_ref'], r[k + '_var'], r[k + '_pct'] = a, b, pct(a, b)
    T.append(r)
for i, ion in enumerate(IONS[1:] + IONS[:1]):                       # radius 18 -> 25 A (Ce last; Ce GB at 25 A failed)
    j = IONS.index(ion)
    for env, site, o in (('bulk', 'bulk atom 2 (4g)', 1), ('GB', 'GB atom 13 (A)', 3)):
        a, b = 'R%02d_%s_%s_R18_FeNeutral' % (4 * j + o, ion, env), 'R%02d_%s_%s_R25_FeNeutral' % (4 * j + o + 1, ion, env)
        if ion == 'Ce' and env == 'GB': continue
        add('radius 18->25 A', ion, site, a, b, CAS[a], CAS[b])
for ion, (b1, b2, g1, g2) in (('Nd', ('R10', 'R25', 'R12', 'R26')), ('Dy', ('R22', 'R27', 'R24', 'R28'))):
    for site, a, b in (('bulk atom 2 (4g)', pick(b1), pick(b2)), ('GB atom 13 (A)', pick(g1), pick(g2))):
        add('charge partition (B -1, Fe -5/14)', ion, site, a, b, CAS[a], CAS[b])
for t, ion, site in (('T13', 'Nd', 'bulk atom 2 (4g)'), ('T14', 'Nd', 'GB atom 13 (A)'), ('T15', 'Dy', 'bulk atom 2 (4g)'), ('T16', 'Dy', 'GB atom 13 (A)')):
    r = pick(t); add('SC-NEVPT2', ion, site, r, r, CAS[r], NEV[r])
for ion in ['Ce', 'Pr', 'Nd', 'Tb', 'Dy', 'Gd']:                     # single non-neutral supercell (published embedding)
    for env, site in (('bulk', 'bulk atom 2 (4g)'), ('GB', 'GB atom 13 (A)')):
        ref = R25[ion][0] if env == 'bulk' else (R25[ion][1] if ion != 'Ce' else pick('R03_Ce_GB_R18'))
        var = 'PUB_%s_%s' % (ion, env)
        if var in CAS: add('single cell (net +9/+21) vs 25 A sphere', ion, site, ref, var, CAS[ref], CAS[var])
for site, ref, var in (('bulk atom 2 (4g)', 'PUB_Nd_bulk', pick('C07')), ('GB atom 13 (A)', 'PUB_Nd_GB', pick('C08'))):
    add('basis TZVP vs QZVP (single cell)', 'Nd', site, ref, var, CAS[ref], CAS[var])
for site, casci, cas in (('bulk atom2', 'T20', 'R02'), ('GBclassB atom9', 'T21', 'T07'), ('GB atom13', 'T22', 'R03')):
    a, b = pick(casci), pick(cas); add('Ce CASSCF vs CASCI(free-ion orbitals)', 'Ce', site, a, b, CAS[a], CAS[b])
with open(os.path.join(outdir, 'sensitivity_tests.csv'), 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(T[0])); wr.writeheader(); wr.writerows(T)
print('wrote matched_pairs_R25.csv (%d rows) and sensitivity_tests.csv (%d rows)' % (len(out), len(T)))
