#!/usr/bin/env python3
"""pointcharge_sites.py - electrostatic (point-charge) crystal-field lattice sums at every rare-earth site
of the bulk and grain-boundary supercells, built with exactly the embedding used for the ORCA runs
(neutral ion-centred sphere of radius R, RE +3 / Fe -3/14 / B -3, residual spread over the outer 2 A shell).

For each site the dimensionless lattice sums
    G^k_q = sum_j q_j C^k_q(R_j) / R_j^(k+1)          (q_j in e, R_j in Angstrom)
are evaluated (k = 2, 4, 6), together with the radial decomposition of G^2_0 into shells.
In the point-charge model the electron feels B^k_q = -e^2 <r^k> (-1)^q G^k_-q, so the ratio of G^2_0 between two
sites estimates the ratio of the axial parameter A_2^0<r^2> at fixed ion.

python3 pointcharge_sites.py <bulk_CONTCAR> <GB_CONTCAR> <outdir> [R]
"""
import sys, os, csv, math
import numpy as np
from scipy.special import sph_harm_y

Q_N = {'Nd': 3.0, 'Fe': -3 / 14, 'B': -3.0}
Q_A = {'Nd': 3.0, 'Fe': -5 / 14, 'B': -1.0}


def contcar(path):
    L = open(path).read().splitlines()
    s = float(L[1]); lat = np.array([[float(x) for x in L[i].split()] for i in (2, 3, 4)]) * s
    el = L[5].split(); n = [int(x) for x in L[6].split()]
    frac = np.array([[float(x) for x in L[8 + i].split()[:3]] for i in range(sum(n))])
    els = sum(([e] * k for e, k in zip(el, n)), [])
    return lat, els, frac @ lat


def sphere(lat, el, pos, idx, q, R, shell=2.0):
    c = pos[idx]; rows = []
    n = int(np.ceil(R / min(np.linalg.norm(lat, axis=1)))) + 1
    for a in range(-n, n + 1):
        for b in range(-n, n + 1):
            for d in range(-n, n + 1):
                sh = np.array([a, b, d]) @ lat
                for e, p in zip(el, pos):
                    r = p + sh; dist = np.linalg.norm(r - c)
                    if 1e-6 < dist <= R: rows.append([q[e], *r, dist])
    rows = np.array(rows); tot = rows[:, 0].sum(); outer = rows[:, 4] > R - shell
    rows[outer, 0] -= tot / outer.sum()
    return c, rows


def Ckq(k, q, th, ph):
    return np.sqrt(4 * np.pi / (2 * k + 1)) * sph_harm_y(k, q, th, ph)


def lattice_sums(c, rows):
    v = rows[:, 1:4] - c; r = rows[:, 4]; qq = rows[:, 0]
    th = np.arccos(np.clip(v[:, 2] / r, -1, 1)); ph = np.arctan2(v[:, 1], v[:, 0])
    G = {}
    for k in (2, 4, 6):
        for q in range(-k, k + 1):
            G[(k, q)] = np.sum(qq * Ckq(k, q, th, ph) / r ** (k + 1))
    # Cartesian field-gradient-like tensor for the rank-2 axis
    T = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            T[a, b] = np.sum(qq * (3 * v[:, a] * v[:, b] - (a == b) * r ** 2) / r ** 5)
    shells = {}
    g20 = (qq * (3 * (v[:, 2] / r) ** 2 - 1) / 2 / r ** 3)
    for lo, hi in ((0, 3.5), (3.5, 6.0), (6.0, 12.0), (12.0, 99)):
        m = (r > lo) & (r <= hi); shells['G20_%g_%g' % (lo, hi)] = g20[m].sum()
    for e, val in (('RE', 3.0), ('B', None)):
        pass
    return G, T, shells


def run(bulk, gb, outdir, R=25.0):
    os.makedirs(outdir, exist_ok=True)
    out = []
    for env, path in (('bulk', bulk), ('GB', gb)):
        lat, el, pos = contcar(path)
        re_idx = [i for i, e in enumerate(el) if e == 'Nd']
        for i in re_idx:
            for tag, q in (('FeNeutral', Q_N), ('FeBalt', Q_A)):
                c, rows = sphere(lat, el, pos, i, q, R)
                G, T, sh = lattice_sums(c, rows)
                w, U = np.linalg.eigh(T); iu = int(np.argmax(np.abs(w)))
                S = {k: math.sqrt(sum(abs(G[(k, qq)]) ** 2 for qq in range(-k, k + 1)) / (2 * k + 1)) for k in (2, 4, 6)}
                # element-resolved first-shell (<3.5 A) contributions to G20
                v = rows[:, 1:4] - c; r = rows[:, 4]
                g20 = rows[:, 0] * (3 * (v[:, 2] / r) ** 2 - 1) / 2 / r ** 3
                d = dict(environment=env, atom=i + 1, charges=tag, z=round(c[2], 3),
                         G20=G[(2, 0)].real, G40=G[(4, 0)].real, G60=G[(6, 0)].real,
                         S2=S[2], S4=S[4], S6=S[6], f_ax=abs(G[(2, 0)].real) / (math.sqrt(5) * S[2]),
                         axis_angle=math.degrees(math.acos(min(1, abs(U[2, iu])))), **sh)
                out.append(d)
    keys = list(out[0].keys())
    with open(os.path.join(outdir, 'pointcharge_sites_R%d.csv' % R), 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=keys); wr.writeheader()
        for d in out: wr.writerow({k: (round(v, 5) if isinstance(v, float) else v) for k, v in d.items()})
    return out


if __name__ == '__main__':
    R = float(sys.argv[4]) if len(sys.argv) > 4 else 25.0
    run(sys.argv[1], sys.argv[2], sys.argv[3], R)
