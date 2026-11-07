#!/usr/bin/env python3
"""Magnetic easy axis from the LF-CI Van Vleck susceptibility tensor at 2 K (principal axis of largest chi),
for every run in the list; angle to the cell c axis (z). Usage: python3 easy_axis_lfci.py parsed out.csv run1 run2 ..."""
import sys, os, math, csv, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import cflib as c
src, out, runs = sys.argv[1], sys.argv[2], sys.argv[3:]
rows = []
for r in runs:
    d, a, Vm, U = c.load(os.path.join(src, r + '.json')); n = c.NEL[d['ion']]
    E, Uc, M = c.lfci(Vm, a['zeta'], a.get('F2ff') or 0, a.get('F4ff') or 0, a.get('F6ff') or 0, n)
    chi = c.chi_tensor(E, M, 2.0); w, V = np.linalg.eigh(chi); ax = V[:, 2]
    rows.append(dict(run=r, ion=d['ion'], chi2K_max=w[2], chi2K_min=w[0], anis=w[2] / max(w[0], 1e-12),
                     easy_angle_deg=math.degrees(math.acos(min(1, abs(ax[2])))),
                     E1_lfci=(E[2] if n % 2 else E[2]), E1_orca=d['soc'][0][2]))
    print(rows[-1], flush=True)
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
