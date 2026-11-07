#!/usr/bin/env python3
"""make_figures.py - Figures 2-5 of the revised manuscript from the analysis tables.

python3 make_figures.py <analysis_dir> <out_dir>
needs: tables/cf_all_runs.csv, tables/easy_axis_lfci.csv, parsed/*.json (for multiplets and 4f densities)
Colours: validated categorical slots (bulk #2a78d6, GB class A #eb6834, GB class B #1baf7a); marker shape is
a second encoding for every series. Ranks k = 2,4,6 use one sequential blue ramp."""
import sys, os, json, math
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import cflib as c

AN, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['STIXGeneral', 'DejaVu Serif'], 'mathtext.fontset': 'stix',
                     'font.size': 9, 'axes.labelsize': 9.5, 'axes.titlesize': 9.5, 'legend.fontsize': 8, 'xtick.labelsize': 8.5,
                     'ytick.labelsize': 8.5, 'axes.linewidth': 0.8, 'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 400, 'figure.dpi': 120})
BULK, GA, GBB, NEU = '#2a78d6', '#eb6834', '#1baf7a', '#8c8b86'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
IONS = ['Ce', 'Pr', 'Nd', 'Gd', 'Tb', 'Dy']
d = pd.read_csv(os.path.join(AN, 'tables/cf_all_runs.csv')); d = d[d.block == 'CASSCF'].set_index('run')
ez = pd.read_csv(os.path.join(AN, 'tables/easy_axis_lfci.csv')).set_index('run')

BULK2 = {i: 'R%02d_%s_bulk_R25_FeNeutral' % (n, i) for i, n in zip(IONS, (2, 6, 10, 14, 18, 22))}
GB13 = {i: 'R%02d_%s_GB_R25_FeNeutral' % (n, i) for i, n in zip(IONS, (4, 8, 12, 16, 20, 24))}
BULK5 = {i: 'T%02d_%s_bulk4f_atom5_R25' % (n, i) for i, n in zip(IONS, range(1, 7))}
GB9 = {i: 'T%02d_%s_GBclassB_atom9_R25' % (n, i) for i, n in zip(IONS, range(7, 13))}
PAIR_A = {i: (BULK5[i], GB13[i]) for i in IONS if i != 'Ce'}          # Ce: 4f-parent sites not obtainable
PAIR_B = {i: (BULK2[i], GB9[i]) for i in IONS}


def label(ax, s):
    ax.text(-0.02, 1.02, s, transform=ax.transAxes, ha='right', va='bottom', fontsize=10.5, fontweight='bold', color=INK)


def style(ax):
    ax.grid(axis='y', color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.tick_params(colors=INK2, length=3); [s.set_color(INK2) for s in ax.spines.values()]


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + '.png'), bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(OUT, name + '.pdf'), bbox_inches='tight', facecolor='white'); plt.close(fig)


def ratio(pair, key):
    return [d.loc[pair[i][1], key] / d.loc[pair[i][0], key] if i in pair else np.nan for i in IONS]


x = np.arange(len(IONS))
# ------------------------------------------------------------------ Figure 2
fig = plt.figure(figsize=(7.0, 5.6))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.30)
for k, (key, ylab) in enumerate((('split', 'GB / bulk, total 4$f$ splitting'), ('A20r2', r'GB / bulk, $A_2^0\langle r^2\rangle$'))):
    ax = fig.add_subplot(gs[0, k]); style(ax)
    ax.axhspan(0.86, 1.14, color='#f0efec', zorder=0, lw=0)
    ax.axhline(1, color=INK2, lw=0.8, ls='--', zorder=1)
    ax.plot(x, ratio(PAIR_A, key), '-o', color=GA, ms=6, lw=1.6, mec='white', mew=0.8, label='class A (4$f$ parent)', zorder=3)
    ax.plot(x, ratio(PAIR_B, key), '-s', color=GBB, ms=6, lw=1.6, mec='white', mew=0.8, label='class B (4$g$ parent)', zorder=3)
    ax.set_xticks(x, IONS); ax.set_ylabel(ylab); ax.set_ylim(0.8, 1.7)
    ax.text(-0.25, 0.83, 'convergence uncertainty (±14%)', fontsize=7, color=INK2, ha='left', va='bottom')
    if k == 0: ax.legend(loc='upper left', frameon=False, handlelength=1.6)
    label(ax, '(a)' if k == 0 else '(b)')
for k, ion in enumerate(('Nd', 'Dy')):
    ax = fig.add_subplot(gs[1, k]); style(ax); ax.grid(False)
    n = {'Nd': 10, 'Dy': 16}[ion]
    for j, (run, col, xc, lab) in enumerate(((BULK5[ion], BULK, 0, 'bulk 4$f$'), (GB13[ion], GA, 1, 'GB class A'))):
        lev = np.array(json.load(open(os.path.join(AN, 'parsed', run + '.json')))['soc'][0][:n])[::2]
        for e in lev: ax.plot([xc - 0.32, xc + 0.32], [e, e], color=col, lw=2.0, solid_capstyle='round')
        ax.annotate('', xy=(xc + 0.42, lev[1]), xytext=(xc + 0.42, lev[0]), arrowprops=dict(arrowstyle='<->', color=INK, lw=0.8, shrinkA=0, shrinkB=0))
        ax.text(xc + 0.47, max(lev[1] / 2, 0.035 * lev[-1]), '$E_1$ = %.0f' % lev[1], fontsize=7.5, va='bottom' if lev[1] < 0.07 * lev[-1] else 'center', color=INK)
        ax.text(xc, -0.07, lab, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=8.5, color=INK)
    ax.set_xlim(-0.6, 1.9); ax.set_xticks([]); ax.spines['bottom'].set_visible(False)
    ax.set_ylabel('energy (cm$^{-1}$)'); ax.set_title(r'%s$^{3+}$, ground multiplet %s' % (ion, {'Nd': r'$^4I_{9/2}$', 'Dy': r'$^6H_{15/2}$'}[ion]), color=INK)
    label(ax, '(c)' if k == 0 else '(d)')
save(fig, 'Figure_2')

# ------------------------------------------------------------------ Figure 3
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.8), gridspec_kw=dict(wspace=0.32))
ax = axs[0]; style(ax)
aniso = [i for i in IONS if i != 'Gd']; xa = np.arange(len(aniso))
KR = {'Ce', 'Nd', 'Dy'}            # Kramers ions: ab initio ground-doublet g axis; Tb: LF-CI chi axis; Pr (quasi-doublet) omitted
def tilt(run, ion):
    return d.loc[run, 'g_angle'] if ion in KR else ez.loc[run, 'easy_angle_deg']
tA = [tilt(PAIR_A[i][1], i) if (i in PAIR_A and i != 'Pr') else np.nan for i in aniso]
tB = [tilt(PAIR_B[i][1], i) if i != 'Pr' else np.nan for i in aniso]
pub = {'Pr': 'PUB_Pr_GB', 'Nd': 'PUB_Nd_GB', 'Dy': 'PUB_Dy_GB'}
tP = [ez.loc[pub[i], 'easy_angle_deg'] if i in pub else np.nan for i in aniso]
ax.plot(xa - 0.12, tA, 'o', color=GA, ms=7, mec='white', mew=0.8, label='class A, converged embedding')
ax.plot(xa + 0.12, tB, 's', color=GBB, ms=7, mec='white', mew=0.8, label='class B, converged embedding')
ax.plot(xa, tP, 'D', mfc='none', mec=NEU, ms=7, mew=1.3, label='class A, single-cell embedding')
ax.set_xticks(xa, aniso); ax.set_ylim(-4, 95); ax.set_yticks([0, 15, 30, 45, 60, 75, 90])
ax.set_ylabel('easy-axis tilt from $c$ (deg)'); ax.legend(loc='center left', frameon=False, fontsize=7.3, bbox_to_anchor=(0.0, 0.62))
label(ax, '(a)')
ax = axs[1]; style(ax)
rA = ratio(PAIR_A, 'E1'); rB = ratio(PAIR_B, 'E1')
rA = [rA[IONS.index(i)] for i in aniso]; rB = [rB[IONS.index(i)] for i in aniso]
ax.axhline(1, color=INK2, lw=0.8, ls='--')
ax.plot(xa, rA, '-o', color=GA, ms=6, lw=1.6, mec='white', mew=0.8, label='class A / bulk 4$f$')
ax.plot(xa, rB, '-s', color=GBB, ms=6, lw=1.6, mec='white', mew=0.8, label='class B / bulk 4$g$')
ax.set_xticks(xa, aniso); ax.set_ylabel('GB / bulk, first excitation $E_1$'); ax.set_ylim(0.5, 2.4)
ax.legend(loc='upper left', frameon=False); label(ax, '(b)')
save(fig, 'Figure_3')

# ------------------------------------------------------------------ Figure 4
RAMP = {2: '#1c5cab', 4: '#5598e7', 6: '#a9cdf7'}
fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.7), sharey=True, gridspec_kw=dict(wspace=0.08))
for ax, pair, ttl, lab in ((axs[0], PAIR_A, 'class A vs bulk 4$f$', '(a)'), (axs[1], PAIR_B, 'class B vs bulk 4$g$', '(b)')):
    style(ax); ax.axhline(0, color=INK2, lw=0.8)
    for j, k in enumerate((2, 4, 6)):
        v = [100 * (d.loc[pair[i][1], 'S%d' % k] / d.loc[pair[i][0], 'S%d' % k] - 1) if i in pair else np.nan for i in IONS]
        ax.bar(x + (j - 1) * 0.26, v, width=0.24, color=RAMP[k], edgecolor='white', linewidth=0.6, label='$k$ = %d' % k)
    ax.set_xticks(x, IONS); ax.set_title(ttl, color=INK); label(ax, lab)
axs[0].set_ylabel(r'change of rank-$k$ strength $S_k$ (%)'); axs[0].legend(frameon=False, loc='upper left', ncol=3)
axs[0].yaxis.set_major_locator(MaxNLocator(7)); axs[0].set_ylim(-28, 55)
save(fig, 'Figure_4')

# ------------------------------------------------------------------ Figure 5 : ground-state 4f charge density
def rdm(run):
    dd, a, Vm, U = c.load(os.path.join(AN, 'parsed', run + '.json')); n = c.NEL[dd['ion']]
    E, Uc, M = c.lfci(Vm, a['zeta'], a.get('F2ff') or 0, a.get('F4ff') or 0, a.get('F6ff') or 0, n)
    nd = 2 if n % 2 else 1
    D = np.zeros((7, 7), complex)
    for i in range(7):
        for j in range(7):
            h = np.zeros((14, 14)); h[2 * i, 2 * j] = h[2 * i + 1, 2 * j + 1] = 1.0
            O = c.one_body(h, n); D[i, j] = np.mean([Uc[:, s].conj() @ O @ Uc[:, s] for s in range(nd)])
    return D


from scipy.special import sph_harm_y
def density_map(D, N=301, rmax=1.4):
    xs = np.linspace(-rmax, rmax, N); X, Z = np.meshgrid(xs, xs)
    # plane spanned by [110] (in-plane) and c
    r = np.sqrt(X ** 2 + Z ** 2) + 1e-12; th = np.arccos(Z / r); ph = np.full_like(th, np.pi / 4); ph[X < 0] = np.pi / 4 + np.pi
    Y = np.array([sph_harm_y(3, m, th, ph) for m in c.MS])            # (7,N,N)
    ang = np.einsum('iab,ij,jab->ab', Y.conj(), D, Y).real
    zeta = 3.0 / 0.55                                                # 4f radial ~ r^3 exp(-zeta r), <r> ~ 0.6 A
    rad = r ** 6 * np.exp(-2 * zeta * r)
    return xs, ang * rad


RUN_FIG5 = True
panels = [('Pr', BULK5['Pr'], GB13['Pr']), ('Nd', BULK5['Nd'], GB13['Nd']), ('Dy', BULK5['Dy'], GB13['Dy'])]
fig, axs = plt.subplots(2, 3, figsize=(7.0, 4.75), gridspec_kw=dict(hspace=0.10, wspace=0.06))
cache = {}
for jc, (ion, rb, rg) in enumerate(panels):
    maps = [density_map(rdm(rb)), density_map(rdm(rg))]
    vmax = max(m[1].max() for m in maps)
    for ir, (xs, rho) in enumerate(maps):
        ax = axs[ir, jc]
        im = ax.imshow(rho / vmax, extent=[xs[0], xs[-1], xs[0], xs[-1]], origin='lower', cmap='Blues', vmin=0, vmax=1)
        ax.contour(xs, xs, rho / vmax, levels=[0.25, 0.5, 0.75], colors=[INK2], linewidths=0.5)
        ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1])
        ax.tick_params(labelsize=7.5, colors=INK2, length=2)
        if ir == 0: ax.set_title(r'%s$^{3+}$' % ion, color=INK); ax.set_xticklabels([])
        else: ax.set_xlabel(r'[110] (Å)')
        if jc == 0: ax.set_ylabel(('bulk 4$f$ site\n' if ir == 0 else 'GB class A site\n') + r'$c$ axis (Å)')
        else: ax.set_yticklabels([])
        ax.annotate('', xy=(1.22, 1.25), xytext=(1.22, 0.75), arrowprops=dict(arrowstyle='->', color=INK, lw=0.9))
        ax.text(1.12, 1.0, '$c$', fontsize=8, ha='right', va='center', color=INK)
        label(ax, '(%s)' % 'abcdef'[ir * 3 + jc])
cb = fig.colorbar(im, ax=axs, shrink=0.75, pad=0.02); cb.set_label('ground-state 4$f$ density (normalised)', color=INK); cb.ax.tick_params(labelsize=7.5)
save(fig, 'Figure_5')
print('figures written to', OUT)
