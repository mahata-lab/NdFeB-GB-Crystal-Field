#!/usr/bin/env python3
"""cf_tables.py - one row per ORCA run with every descriptor used in the paper.

python3 cf_tables.py <parsed_dir> <out.csv>
Descriptors (cm-1 unless stated):
  split            spread of the seven AILFT one-electron 4f energies (total 4f ligand-field splitting)
  orb_d1, orb_d2   first and second one-electron orbital gaps (E2-E1, E3-E1)
  B20,B40,B60      Wybourne one-electron parameters along the cell c axis (z); A20r2 = B20/2 (Stevens A_2^0<r^2>)
  S2,S4,S6         rotational invariants of rank k; f_ax = |B20|/(sqrt5 S2) = axial fraction of the rank-2 part
  ax2_angle        angle (deg) between c and the unique principal axis of the rank-2 CF tensor
  E1               first excitation of the ground J multiplet from the ab initio QDPT-SOC spectrum
                   (first excited Kramers doublet for Ce, Nd, Gd, Dy; first level above the ground quasi-doublet for Pr, Tb)
  tun              ground quasi-doublet splitting (non-Kramers ions)
  Jspan            spread of the ground J multiplet (QDPT-SOC)
  gx,gy,gz, g_angle  ground Kramers-doublet g tensor from ORCA's <L>,<S> matrix elements (gauge-free) and angle of gz to c
"""
import sys, os, glob, csv, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cflib as c

TWOJ1 = dict(Ce=6, Pr=9, Nd=10, Gd=8, Tb=13, Dy=16)       # 2J+1 of the ground multiplet
KRAMERS = dict(Ce=True, Pr=False, Nd=True, Gd=True, Tb=False, Dy=True)


def gtensor_kramers(d):
    ge = c.GE
    M1 = np.array(d['Re_L_1']) + ge * np.array(d['Re_S_1'])
    M2r = np.array(d['Re_L_2']) + ge * np.array(d['Re_S_2'])
    M2i = np.array(d['Im_L_2']) + ge * np.array(d['Im_S_2'])
    g = 2 * np.array([M2r, -M2i, M1])           # rows: pseudo-spin x,y,z ; columns: lab x,y,z
    w, U = np.linalg.eigh(g.T @ g)
    return np.sqrt(np.clip(w, 0, None)), U


def row(path, block='CASSCF'):
    d, a, Vm, U = c.load(path, block)
    ion = d['ion']; r = dict(run=d['label'], ion=ion, block=block, normal=d['normal'])
    E = np.sort(np.array(a['E'])); r['split'] = E[-1] - E[0]; r['orb_d1'] = E[1] - E[0]; r['orb_d2'] = E[2] - E[0]
    B, res = c.wybourne(Vm); S = c.invariants(B)
    for k in (2, 4, 6):
        r['B%d0' % k] = B[(k, 0)].real; r['S%d' % k] = S[k]
    r['A20r2'] = B[(2, 0)].real / c.LAMBDA[2]
    r['A40r4'] = B[(4, 0)].real / c.LAMBDA[4]; r['A60r6'] = B[(6, 0)].real / c.LAMBDA[6]
    r['f_ax'] = abs(B[(2, 0)].real) / (math.sqrt(5) * S[2])
    D, w, Ua = c.rank2_axes(B)
    # unique axis = eigenvalue farthest from the mean of the other two
    dev = [abs(w[i] - np.mean(np.delete(w, i))) for i in range(3)]
    iu = int(np.argmax(dev)); ax = Ua[:, iu]
    r['ax2_angle'] = math.degrees(math.acos(min(1, abs(ax[2]))))
    r['rank2_eta'] = (w[1] - w[0]) / (w[2] - w[0]) if w[2] != w[0] else 0.0
    r['wyb_resid'] = res
    soc = d['soc'][0 if block == 'CASSCF' else 1] if d['soc'] else []
    if soc:
        n = TWOJ1[ion]; lev = np.array(soc[:n])
        r['Jspan'] = lev[-1] - lev[0]
        if KRAMERS[ion]:
            r['E1'] = lev[2]; r['E2'] = lev[4] if n > 4 else np.nan; r['tun'] = 0.0
        else:
            r['tun'] = lev[1]; r['E1'] = lev[2]; r['E2'] = lev[3]
        r['gapJ'] = soc[n] - lev[-1] if len(soc) > n else np.nan
    if d.get('kramers') and len(d['kramers']) == 6 and block == 'CASSCF' and KRAMERS[ion]:
        g, Ug = gtensor_kramers(d['kramers'])
        r['gx'], r['gy'], r['gz'] = g; r['g_angle'] = math.degrees(math.acos(min(1, abs(Ug[2, 2]))))
    return r


if __name__ == '__main__':
    src, out = sys.argv[1], sys.argv[2]
    rows = []
    for f in sorted(glob.glob(os.path.join(src, '*.json'))):
        d = json.load(open(f))
        rows.append(row(f))
        if 'NEVPT2' in d['ailft']:
            rows.append(row(f, 'NEVPT2'))
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    print('wrote', out, len(rows), 'rows')
