#!/usr/bin/env python3
"""make_supplementary.py - writes the LaTeX tables of the Supplementary Material from the analysis data.
python3 make_supplementary.py <analysis_dir> <bulk_CONTCAR> <GB_CONTCAR> <site_coordination.csv> <out.tex-fragment>"""
import os, sys, csv, numpy as np, pandas as pd
AN, BULK, GB, SITES, OUT = sys.argv[1:6]


def contcar(p):
    L = open(p).read().splitlines(); s = float(L[1])
    lat = np.array([[float(x) for x in L[i].split()] for i in (2, 3, 4)]) * s
    el = L[5].split(); n = [int(x) for x in L[6].split()]
    fr = np.array([[float(x) for x in L[8 + i].split()[:3]] for i in range(sum(n))])
    return lat, sum(([e] * k for e, k in zip(el, n)), []), fr


def bulk_label(i, f):
    # Nd sublattice labels in P4_2/mnm: 4f (x,x,0) family, 4g (x,-x,0) family; z = 1/2 layer via 4_2
    x, y, z = f
    if abs(z) < 0.1 or abs(z - 1) < 0.1: return '4f' if abs(((x - y) + 0.5) % 1 - 0.5) < 0.05 else '4g'
    return '4g' if abs(((x - y) + 0.5) % 1 - 0.5) < 0.05 else '4f'


GBCLASS = {**{a: 'A' for a in (5, 8, 13, 16)}, **{a: 'B' for a in (1, 4, 9, 12)}, **{a: 'C' for a in (2, 3, 10, 11)}, **{a: 'D' for a in (6, 7, 14, 15)}}


def coord_table(path, title, label, gb=False):
    lat, el, fr = contcar(path); cart = fr @ lat
    rows = []; cnt = {}
    for i, (e, f, c) in enumerate(zip(el, fr, cart), 1):
        cnt[e] = cnt.get(e, 0) + 1
        tag = ''
        if e == 'Nd': tag = ('class ' + GBCLASS[i]) if gb else bulk_label(i, f)
        rows.append('%d & %s & %.5f & %.5f & %.5f & %.4f & %.4f & %.4f & %s \\\\' % (i, e, *f, *c, tag))
    head = ('\\begin{longtable}{rlrrrrrrl}\n\\caption{%s Lattice vectors (\\AA): $a_1$ = (%.6f, %.6f, %.6f), $a_2$ = (%.6f, %.6f, %.6f), $a_3$ = (%.6f, %.6f, %.6f). '
            'Fractional ($x,y,z$) and Cartesian ($X,Y,Z$, \\AA) coordinates; atom indices are those used in the main text.}\\label{%s}\\\\\n'
            '\\toprule\n\\# & El. & $x$ & $y$ & $z$ & $X$ & $Y$ & $Z$ & Site \\\\\n\\midrule\n\\endfirsthead\n'
            '\\toprule\n\\# & El. & $x$ & $y$ & $z$ & $X$ & $Y$ & $Z$ & Site \\\\\n\\midrule\n\\endhead\n') % (title, *lat[0], *lat[1], *lat[2], label)
    return head + '\n'.join(rows) + '\n\\bottomrule\n\\end{longtable}\n'


out = []
out.append(coord_table(BULK, 'Relaxed bulk Nd$_8$Fe$_{56}$B$_4$ supercell (68 atoms, space group $P4_2/mnm$).', 'tab:S1'))
out.append(coord_table(GB, 'Relaxed $90^\\circ$ [001] twist grain-boundary supercell Nd$_{16}$Fe$_{112}$B$_8$ (136 atoms). Boundary planes at $z = 0$ and $z = 1/2$.', 'tab:S2', gb=True))

# S3 site coordination
sc = list(csv.DictReader(open(SITES)))
lat_b, el_b, fr_b = contcar(BULK)
r = ['\\begin{longtable}{llllcl}\n\\caption{Coordination of every rare-earth site: number of Fe neighbors within 3.3~\\AA{} and the eight shortest RE--Fe distances (periodic images included).}\\label{tab:S3}\\\\\n\\toprule\nEnvironment & Atom & Site & $z$ (\\AA) & $N_{\\mathrm{Fe}}$ & shortest RE--Fe distances (\\AA) \\\\\n\\midrule\n\\endfirsthead\n\\toprule\nEnvironment & Atom & Site & $z$ (\\AA) & $N_{\\mathrm{Fe}}$ & shortest RE--Fe distances (\\AA) \\\\\n\\midrule\n\\endhead']
for s in sc:
    a = int(s['atom_1based']); z = float(s['xyz'].split()[2])
    site = bulk_label(a, fr_b[a - 1]) if s['environment'] == 'bulk' else 'class ' + GBCLASS[a]
    dd = ', '.join('%.3f' % float(x) for x in s['nearest_16_Fe_A'].split()[:8])
    r.append('%s & %d & %s & %.3f & %s & %s \\\\' % (s['environment'], a, site, z, s['Fe_within_3p3A'], dd))
r.append('\\bottomrule\n\\end{longtable}\n'); out.append('\n'.join(r))

# S4 crystal-field parameters (production matched runs)
d = pd.read_csv(AN + '/tables/cf_all_runs.csv'); d = d[d.block == 'CASSCF'].set_index('run')
w = pd.read_csv(AN + '/tables/wybourne_Bkq_all_runs.csv'); w = w[w.block == 'CASSCF'].set_index('run')
IONS = ['Ce', 'Pr', 'Nd', 'Gd', 'Tb', 'Dy']
runs = []
for i, n in zip(IONS, range(1, 7)):
    for lab, rn in (('bulk $4f$ (5)', 'T%02d_%s_bulk4f_atom5_R25' % (n, i)), ('bulk $4g$ (2)', 'R%02d_%s_bulk_R25_FeNeutral' % (4 * n - 2, i)),
                    ('GB A (13)', 'R%02d_%s_GB_R25_FeNeutral' % (4 * n, i)), ('GB B (9)', 'T%02d_%s_GBclassB_atom9_R25' % (n + 6, i))):
        if rn in d.index: runs.append((i, lab, rn))
t = ['\\begin{longtable}{lllrrrrrrrrrr}\n\\caption{Crystal-field parameters of the production calculations (SA-CASSCF, $R = 25$~\\AA, cm$^{-1}$, cell frame $z\\parallel c$, $x\\parallel a$). '
     'Wybourne parameters $B^k_q$ of the one-electron AILFT operator $V=\\sum_{kq}B^k_qC^k_q$; for $q\\neq 0$ the real (Re) and imaginary (Im) parts are listed, with $B^k_{-q} = (-1)^q\\,\\overline{B^k_q}$. '
     'The boundary sites have lower symmetry than the bulk sites, so further components (up to about 600~cm$^{-1}$, e.g. $B^4_2$, $B^6_2$) are non-zero; the complete complex sets for all 85 calculations are given in the data repository (wybourne\\_Bkq\\_all\\_runs.csv). '
     '$S_k$: rank-$k$ strength; $\\zeta$: spin--orbit constant; $E_1$: first excitation of the ground multiplet.}\\label{tab:S4}\\\\\n\\toprule\n'
     'Ion & Site & & $B^2_0$ & Im$B^2_2$ & $B^4_0$ & Re$B^4_4$ & $B^6_0$ & Re$B^6_4$ & $S_2$ & $S_4$ & $S_6$ & $E_1$ \\\\\n\\midrule\n\\endfirsthead\n\\toprule\n'
     'Ion & Site & & $B^2_0$ & Im$B^2_2$ & $B^4_0$ & Re$B^4_4$ & $B^6_0$ & Re$B^6_4$ & $S_2$ & $S_4$ & $S_6$ & $E_1$ \\\\\n\\midrule\n\\endhead']
for i, lab, rn in runs:
    a, b = w.loc[rn], d.loc[rn]
    e1 = '--' if i == 'Gd' else '%.0f' % b.E1
    t.append('%s & %s & & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %s \\\\' % (i, lab, a.Re_B20, a.Im_B22, a.Re_B40, a.Re_B44, a.Re_B60, a.Re_B64, b.S2, b.S4, b.S6, e1))
t.append('\\bottomrule\n\\end{longtable}\n'); out.append('\n'.join(t))
# check that unlisted components are small
cols = [c for c in w.columns if c.startswith(('Re_B', 'Im_B')) and c not in ('Re_B20', 'Im_B22', 'Re_B40', 'Re_B44', 'Re_B60', 'Re_B64')]
mx = w.loc[[rn for _, _, rn in runs], cols].abs().max().max()
print('largest unlisted |B| among production runs: %.1f' % mx, w.loc[[rn for _, _, rn in runs], cols].abs().max().sort_values().tail(4).to_dict())

# S5 sensitivity
s = pd.read_csv(AN + '/tables/sensitivity_tests.csv')
s = s[~s.test.str.startswith('Ce CASSCF')]
tt = ['\\begin{longtable}{llllrrrrrr}\n\\caption{Convergence and sensitivity tests (cm$^{-1}$). Reference = production model unless stated; variant as given in the first column.}\\label{tab:S5}\\\\\n\\toprule\nTest & Ion & Site & & split ref & split var & \\Atwo{} ref & \\Atwo{} var & $E_1$ ref & $E_1$ var \\\\\n\\midrule\n\\endfirsthead\n\\toprule\nTest & Ion & Site & & split ref & split var & \\Atwo{} ref & \\Atwo{} var & $E_1$ ref & $E_1$ var \\\\\n\\midrule\n\\endhead']
NAME = {'radius 18->25 A': 'radius 18$\\to$25~\\AA', 'charge partition (B -1, Fe -5/14)': 'B $-1$, Fe $-5/14$', 'SC-NEVPT2': 'SC-NEVPT2',
        'single cell (net +9/+21) vs 25 A sphere': 'single cell', 'basis TZVP vs QZVP (single cell)': 'TZVP (single cell)'}
for _, x in s.iterrows():
    tt.append('%s & %s & %s & & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f \\\\' % (NAME[x.test], x.ion, x.site.replace('(4g)', '($4g$)'), x.split_ref, x.split_var, x.A20r2_ref, x.A20r2_var, x.E1_ref, x.E1_var))
tt.append('\\bottomrule\n\\end{longtable}\n'); out.append('\n'.join(tt))

# S6 Ce CASCI
ce = [('bulk $4g$ (2), 25~\\AA', 'R02_Ce_bulk_R25_FeNeutral', 'T20_Ce_bulk_atom2_R25_CASCI_check'),
      ('GB B (9), 25~\\AA', 'T07_Ce_GBclassB_atom9_R25', 'T21_Ce_GBclassB_atom9_R25_CASCI_check'),
      ('GB A (13), 18~\\AA', 'R03_Ce_GB_R18_FeNeutral', 'T22_Ce_GB_atom13_R18_CASCI_check'),
      ('bulk $4f$ (5), 25~\\AA', None, 'T18_Ce_bulk4f_atom5_R25_CASCI'), ('GB A (13), 25~\\AA', None, 'T19_Ce_GB_atom13_R25_CASCI')]
c6 = ['\\begin{table}[H]\n\\centering\n\\caption{Ce: SA-CASSCF versus CASCI on the free-ion $4f$ orbitals (cm$^{-1}$). At the $4f$-parent sites (bulk atom 5, GB atom 13 at 25~\\AA) SA-CASSCF in the field rotates one active orbital into $6s$ character; the CASCI substitute deviates strongly where both are available and is not used.}\\label{tab:S6}\n\\small\n\\begin{tabular}{lrrrrrr}\n\\toprule\nSite & split CASSCF & split CASCI & \\Atwo{} CASSCF & \\Atwo{} CASCI & $E_1$ CASSCF & $E_1$ CASCI \\\\\n\\midrule']
for lab, a, b in ce:
    B = d.loc[b]
    if a: A = d.loc[a]; c6.append('%s & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f \\\\' % (lab, A.split, B.split, A.A20r2, B.A20r2, A.E1, B.E1))
    else: c6.append('%s & failed & %.0f & failed & %.0f & failed & %.0f \\\\' % (lab, B.split, B.A20r2, B.E1))
c6.append('\\bottomrule\n\\end{tabular}\n\\end{table}\n'); out.append('\n'.join(c6))

# S7 nucleation
nu = pd.read_csv(AN + '/nucleation/nucleation_1d.csv' if os.path.exists(AN + '/nucleation/nucleation_1d.csv') else AN + '/tables/nucleation_1d.csv')
c7 = ['\\begin{table}[H]\n\\centering\n\\caption{One-dimensional nucleation-field estimates for a planar layer of thickness $t$ and anisotropy $K_{\\mathrm{layer}}$ in Nd$_2$Fe$_{14}$B ($A$ = 8~pJ\\,m$^{-1}$, $K_1$ = 4.9~MJ\\,m$^{-3}$, $\\mu_0M_s$ = 1.61~T, $\\mu_0H_A$ = 7.65~T, wall width 4.0~nm). Computed perturbations scale the rare-earth share of $K_1$ ($K_1 - 1.1$~MJ\\,m$^{-3}$) with the change of \\Atwo.}\\label{tab:S7}\n\\small\n\\begin{tabular}{lrrrr}\n\\toprule\nCase & $t$ (nm) & $K_{\\mathrm{layer}}$ (MJ\\,m$^{-3}$) & $\\mu_0H_N$ (T) & $H_N/H_A$ \\\\\n\\midrule']
for _, x in nu.iterrows():
    c7.append('%s & %.1f & %.2f & %.2f & %.3f \\\\' % (x.case.replace('%', '\\%'), x.layer_nm, x.K_layer_MJm3, x.muH_N_T, x.HN_over_HA))
c7.append('\\bottomrule\n\\end{tabular}\n\\end{table}\n'); out.append('\n'.join(c7))
# S8 capped effective-core-potential test (def2-TZVP / def2-ECP level, spin-free SA-CASSCF + AILFT)
SITE8 = [('bulk $4f$ (5)', 'bulk4f5'), ('GB A (13)', 'GBA13'), ('bulk $4g$ (2)', 'bulk4g2'), ('GB B (9)', 'GBB9')]
DKH = {('Nd', 'bulk4f5'): 'T03_Nd_bulk4f_atom5_R25', ('Nd', 'GBA13'): 'R12_Nd_GB_R25_FeNeutral', ('Nd', 'bulk4g2'): 'R10_Nd_bulk_R25_FeNeutral',
       ('Nd', 'GBB9'): 'T09_Nd_GBclassB_atom9_R25', ('Dy', 'bulk4f5'): 'T06_Dy_bulk4f_atom5_R25', ('Dy', 'GBA13'): 'R24_Dy_GB_R25_FeNeutral',
       ('Dy', 'bulk4g2'): 'R22_Dy_bulk_R25_FeNeutral', ('Dy', 'GBB9'): 'T12_Dy_GBclassB_atom9_R25'}
c8 = ['\\begin{table}[H]\n\\centering\n\\caption{Capped effective-core-potential (cECP) test (cm$^{-1}$). All Fe and rare-earth neighbors within 4.0~\\AA{} of the ion (18--19 centers) are replaced by capped ECPs with the same formal charges (Hay--Wadt for Fe, def2 small-core for the rare earth); pc = the same calculation with bare point charges. Level: def2-TZVP with the def2 small-core ECP on the central ion, spin-free SA-CASSCF + AILFT, $R = 25$~\\AA. The last column is the production DKH2/SARC2-DKH-QZVP value with bare point charges (Table~S4).}\\label{tab:S8}\n\\small\n\\begin{tabular}{llrrrrrrr}\n\\toprule\nIon & Site & split pc & split cECP & \\Atwo{} pc & \\Atwo{} cECP & $S_6$ pc & $S_6$ cECP & \\Atwo{} DKH \\\\\n\\midrule']
for ion, base in (('Nd', 1), ('Dy', 9)):
    for i, (lab, key) in enumerate(SITE8):
        pc = d.loc[[r for r in d.index if r.startswith('U%02d_' % (base + 2 * i))][0]]
        ce = d.loc[[r for r in d.index if r.startswith('U%02d_' % (base + 2 * i + 1))][0]]
        c8.append('%s & %s & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f \\\\' % (ion if i == 0 else '', lab, pc.split, ce.split, pc.A20r2, ce.A20r2, pc.S6, ce.S6, d.loc[DKH[(ion, key)]].A20r2))
c8.append('\\bottomrule\n\\end{tabular}\n\\end{table}\n'); out.append('\n'.join(c8))

# S9 mid-grain and substituted boundary-plane sites (production level)
runs9 = [('Nd', 'C (2), mid-grain', 'W01_Nd_GBC2'), ('Nd', 'D (6), mid-grain', 'W02_Nd_GBD6'),
         ('Nd', 'A (13), Zr for RE', 'V02_Nd_GBA13_Zrsub'), ('Nd', 'B (9), B for Fe', 'V03_Nd_GBB9_Bsub'), ('Nd', 'B (9), Zr for RE', 'V04_Nd_GBB9_Zrsub'),
         ('Dy', 'C (2), mid-grain', 'W03_Dy_GBC2'), ('Dy', 'D (6), mid-grain', 'W04_Dy_GBD6'),
         ('Dy', 'A (13), B for Fe', 'V05_Dy_GBA13_Bsub'), ('Dy', 'A (13), Zr for RE', 'V06_Dy_GBA13_Zrsub'),
         ('Dy', 'B (9), B for Fe', 'V07_Dy_GBB9_Bsub'), ('Dy', 'B (9), Zr for RE', 'V08_Dy_GBB9_Zrsub')]
c9 = ['\\begin{table}[H]\n\\centering\n\\caption{Principal Wybourne parameters $B^k_q$, rank invariants $S_k$, total $4f$ splitting, and ground-doublet $g_z$ at the mid-grain sites and at the substituted boundary-plane sites (cm$^{-1}$; production model, $R = 25$~\\AA, SA-CASSCF). Substitutions: B for the nearest Fe (2.77~\\AA{} at A, 2.82~\\AA{} at B) or Zr for the nearest rare-earth neighbor (3.62~\\AA), one per grain-boundary cell, unrelaxed. Nd at class A with B substitution did not converge to a $4f^3$ active space and is omitted. Derived quantities are in Table~4 of the main text.}\\label{tab:S9}\n\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n\\begin{tabular}{llrrrrrrrrrrr}\n\\toprule\nIon & GB site & $B^2_0$ & Im $B^2_2$ & $B^4_0$ & $B^4_4$ & $B^6_0$ & $B^6_4$ & $S_2$ & $S_4$ & $S_6$ & split & $g_z$ \\\\\n\\midrule']
for ion, lab, rn in runs9:
    a, b = w.loc[rn], d.loc[rn]
    c9.append('%s & %s & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.0f & %.2f \\\\' % (ion, lab, a.Re_B20, a.Im_B22, a.Re_B40, a.Re_B44, a.Re_B60, a.Re_B64, b.S2, b.S4, b.S6, b.split, b.gz))
c9.append('\\bottomrule\n\\end{tabular}\n\\end{table}\n'); out.append('\n'.join(c9))
open(OUT, 'w').write('\n'.join(out))
print('written', OUT)
