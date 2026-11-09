#!/usr/bin/env python3
"""nucleation_1d.py - 1-D micromagnetic (Brown/Kronmueller) nucleation field of a planar layer of modified
uniaxial anisotropy embedded in Nd2Fe14B. Linearised Brown equation for small deviations phi(x) of the
magnetisation from the easy axis, field antiparallel:
      -2A phi'' + 2K(x) phi = mu0 Ms H_N phi      ->  mu0 H_N = lambda_min / Ms
(demagnetising and misalignment effects excluded, i.e. the ideal-grain limit in which uniform K gives H_N = H_A).
Parameters (300 K): A = 8 pJ/m, K1 = 4.9 MJ/m3, mu0 Ms = 1.61 T  (H_A = 2K1/Ms ~ 7.6 T, wall width ~ 4 nm).
RE-sublattice share of K1: K1 - K1(Y2Fe14B) = 4.9 - 1.1 MJ/m3, scaled linearly with the computed A20<r2> change."""
import numpy as np, csv, sys
from scipy.linalg import eigh_tridiagonal
A, K1, J = 8e-12, 4.9e6, 1.61; Ms = J / (4e-7 * np.pi); KFe = 1.1e6
def HN(t_nm, Klayer, L=200e-9, n=20001):
    x = np.linspace(-L / 2, L / 2, n); h = x[1] - x[0]
    K = np.where(np.abs(x) <= t_nm * 1e-9 / 2, Klayer, K1)
    d = 2 * A * 2 / h ** 2 + 2 * K; e = -2 * A / h ** 2 * np.ones(n - 1)
    lam = eigh_tridiagonal(d, e, select='i', select_range=(0, 0))[0][0]
    return lam / Ms                               # mu0 H_N in tesla (lambda/Ms)
HA = HN(0.0, K1)                 # uniform grain on the same grid (= 2K1/Ms up to box confinement)
rows = []
# (i) the computed boundary perturbation: change of A20<r2> applied to the RE share of K1
for lab, dA in (('class A, Nd (+38%)', 0.378), ('class B, Nd (+1.5%)', 0.015), ('mean of A and B, Nd (+19%)', 0.196),
                ('hypothetical -30% (mismatched-site comparison)', -0.30)):
    Kl = KFe + (K1 - KFe) * (1 + dA)
    for t in (0.6, 1.2):
        rows.append(dict(case=lab, layer_nm=t, K_layer_MJm3=Kl / 1e6, muH_N_T=HN(t, Kl), HN_over_HA=HN(t, Kl) / HA))
# (ii) what it takes: fully or partly softened layers
for kap in (0.0, 0.25, 0.5):
    for t in (0.6, 1.2, 2.0, 4.0, 8.0):
        Kl = kap * K1; rows.append(dict(case='soft layer K=%.2fK1' % kap, layer_nm=t, K_layer_MJm3=Kl / 1e6, muH_N_T=HN(t, Kl), HN_over_HA=HN(t, Kl) / HA))
with open(sys.argv[1] if len(sys.argv) > 1 else 'nucleation_1d.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print('mu0 H_A = %.3f T (analytic 2K1/Ms = %.3f T), wall width = %.2f nm' % (HA, 2 * K1 / Ms, np.pi * np.sqrt(A / K1) * 1e9))
for r in rows: print('%-48s t=%4.1f nm  K=%.2f MJ/m3  mu0HN=%.2f T  HN/HA=%.3f' % (r['case'], r['layer_nm'], r['K_layer_MJm3'], r['muH_N_T'], r['HN_over_HA']))
