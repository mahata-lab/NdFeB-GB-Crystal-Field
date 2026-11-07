#!/usr/bin/env python3
"""cflib.py - crystal-field analysis of ORCA AILFT results for Ln(III) 4f^n ions.

Works on the JSON files written by orca_parse.py. Provides
  * the AILFT one-electron matrix in the complex |l=3, m> basis (Condon-Shortley phases fixed from
    ORCA's own LF-SOC angular-momentum integrals, so no assumption about ORCA's real-f phase convention);
  * the exact Wybourne decomposition V = sum_kq B^k_q C^k_q (k = 0,2,4,6), rotational invariants
    S_k = sqrt(sum_q |B^k_q|^2 / (2k+1)), Stevens-type A_k^0<r^k> = B^k_0 / lambda_k0 (lambda = 2, 8, 16);
  * the principal axes of the rank-2 (quadrupolar) part;
  * a full ligand-field configuration interaction (LF-CI) of the 4f^n configuration with the AILFT
    one-electron matrix, Slater F2/F4/F6 and zeta, and the Zeeman operator mu = -muB (L + g_e S),
    from which the multiplet levels, ground-doublet g tensor, the Van Vleck susceptibility tensor chi(T)
    and the magnetic easy axis are obtained.
Energies in cm-1, chi in cm3 mol-1 (cgs-emu)."""
import json, itertools, math
from functools import lru_cache
import numpy as np

L = 3
MS = list(range(-L, L + 1))                     # m = -3..3
GE = 2.00231930436
CM2K = 1.438776877                              # cm-1 -> K
NAMUB2_KB = 0.375148                            # N_A muB^2 / k_B  [cm3 K mol-1]
NEL = dict(Ce=1, Pr=2, Nd=3, Gd=7, Tb=8, Dy=9)
LAMBDA = {2: 2.0, 4: 8.0, 6: 16.0}              # B^k_0 (Wybourne) = lambda_k A_k^0 <r^k>


# ----------------------------------------------------------------- angular-momentum algebra
@lru_cache(None)
def fact(n):
    return math.factorial(n)


def wigner3j(j1, j2, j3, m1, m2, m3):
    """Racah formula, integer arguments only (sufficient for l = 3, k <= 6)."""
    if m1 + m2 + m3 != 0 or abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    t1 = j2 - m1 - j3; t2 = j1 + m2 - j3; t3 = j1 + j2 - j3; t4 = j1 - m1; t5 = j2 + m2
    tmin = max(0, t1, t2); tmax = min(t3, t4, t5)
    tri = fact(j1 + j2 - j3) * fact(j1 - j2 + j3) * fact(-j1 + j2 + j3) / fact(j1 + j2 + j3 + 1)
    pre = math.sqrt(tri * fact(j1 + m1) * fact(j1 - m1) * fact(j2 + m2) * fact(j2 - m2) * fact(j3 + m3) * fact(j3 - m3))
    s = 0.0
    for t in range(tmin, tmax + 1):
        s += (-1) ** t / (fact(t) * fact(t - t1) * fact(t - t2) * fact(t3 - t) * fact(t4 - t) * fact(t5 - t))
    return (-1) ** (j1 - j2 - m3) * pre * s


def ck(k, m, mp):
    """c^k(l m, l m') = (-1)^m (2l+1) (l k l; 0 0 0)(l k l; -m, m-m', m')."""
    return (-1) ** m * (2 * L + 1) * wigner3j(L, k, L, 0, 0, 0) * wigner3j(L, k, L, -m, m - mp, mp)


def Ckq(k, q):
    """7x7 matrix of the Racah tensor C^k_q in the |m> basis (m = -3..3)."""
    M = np.zeros((7, 7))
    for i, m in enumerate(MS):
        for j, mp in enumerate(MS):
            M[i, j] = (-1) ** m * (2 * L + 1) * wigner3j(L, k, L, 0, 0, 0) * wigner3j(L, k, L, -m, q, mp)
    return M


def lmat():
    lz = np.diag(np.array(MS, float)).astype(complex)
    lp = np.zeros((7, 7), complex)
    for i, m in enumerate(MS[:-1]):
        lp[i + 1, i] = math.sqrt(L * (L + 1) - m * (m + 1))
    lm = lp.conj().T
    return (lp + lm) / 2, (lp - lm) / 2j, lz


LX, LY, LZ = lmat()


# ----------------------------------------------------------------- ORCA real basis -> |m>
def orca_L(Lsoc):
    """L_x, L_y, L_z (Hermitian) in ORCA's real-f basis from the printed zeta*<i|l|j> integrals."""
    X, Y, Z = [np.array(Lsoc[a]) for a in 'XYZ']
    zeta = np.abs(Z).max() / 3.0
    Ls = [1j * M / zeta for M in (X, Y, Z)]
    for sx, sy, sz in itertools.product((1, -1), repeat=3):
        Lx, Ly, Lz = sx * Ls[0], sy * Ls[1], sz * Ls[2]
        if np.allclose(Lx @ Ly - Ly @ Lx, 1j * Lz, atol=2e-3) and np.allclose(Ly @ Lz - Lz @ Ly, 1j * Lx, atol=2e-3):
            return Lx, Ly, Lz
    raise RuntimeError('angular momentum commutators not satisfied')


def real_to_m(Lsoc):
    """Unitary U (columns = |m>, m=-3..3, in ORCA real basis) with Condon-Shortley phases."""
    Lx, Ly, Lz = orca_L(Lsoc)
    w, v = np.linalg.eigh(Lz)
    U = np.zeros((7, 7), complex)
    i0 = np.argmin(np.abs(w + 3)); U[:, 0] = v[:, i0]
    Lp = Lx + 1j * Ly
    for i in range(6):
        u = Lp @ U[:, i]; U[:, i + 1] = u / np.linalg.norm(u)
    # checks
    assert np.allclose(U.conj().T @ Lz @ U, LZ, atol=2e-3)
    assert np.allclose(U.conj().T @ Lx @ U, LX, atol=2e-3)
    return U


# ----------------------------------------------------------------- CF decomposition
def wybourne(Vm):
    """Exact decomposition of a 7x7 Hermitian matrix (|m> basis) into B^k_q C^k_q, k = 0,2,4,6."""
    B = {}
    for k in (0, 2, 4, 6):
        for q in range(-k, k + 1):
            C = Ckq(k, q)
            B[(k, q)] = np.trace(C.conj().T @ Vm) / np.trace(C.conj().T @ C)
    rec = sum(B[kq] * Ckq(*kq) for kq in B)
    resid = np.abs(rec - Vm).max()
    return B, resid


def invariants(B):
    S = {k: math.sqrt(sum(abs(B[(k, q)]) ** 2 for q in range(-k, k + 1)) / (2 * k + 1)) for k in (2, 4, 6)}
    return S


def rank2_axes(B):
    """Cartesian quadrupole tensor D of the rank-2 part, V2 = sum_ab D_ab l_a l_b (traceless); returns
    eigenvalues (ascending) and eigenvectors (columns)."""
    V2 = sum(B[(2, q)] * Ckq(2, q) for q in range(-2, 3))
    Ls = (LX, LY, LZ)
    ops, idx = [], []
    for a in range(3):
        for b in range(a, 3):
            ops.append(0.5 * (Ls[a] @ Ls[b] + Ls[b] @ Ls[a])); idx.append((a, b))
    # least squares for symmetric D (6 comps) with constant absorbed
    A = np.array([o.ravel() for o in ops] + [np.eye(7).ravel()]).T
    c, *_ = np.linalg.lstsq(A, V2.ravel(), rcond=None)
    D = np.zeros((3, 3))
    for (a, b), x in zip(idx, c[:6].real):
        D[a, b] = D[b, a] = x
    D -= np.trace(D) / 3 * np.eye(3)
    w, U = np.linalg.eigh(D)
    return D, w, U


# ----------------------------------------------------------------- LF-CI of 4f^n
def spinorb_ops(Vm, zeta):
    """14x14 one-electron Hamiltonian, index p = 2*i_m + s (s=0 up, 1 down)."""
    sx = np.array([[0, 1], [1, 0]]) / 2; sy = np.array([[0, -1j], [1j, 0]]) / 2; sz = np.array([[1, 0], [0, -1]]) / 2
    h = np.kron(Vm, np.eye(2)) + zeta * (np.kron(LX, sx) + np.kron(LY, sy) + np.kron(LZ, sz))
    mu = [np.kron(La, np.eye(2)) + GE * np.kron(np.eye(7), s) for La, s in ((LX, sx), (LY, sy), (LZ, sz))]
    return h, mu


@lru_cache(None)
def twoel_table(F2, F4, F6):
    """<pq|rs> for spin orbitals (physicist notation), only nonzero entries."""
    Fk = {2: F2, 4: F4, 6: F6}
    g = {}
    for p, q, r, s in itertools.product(range(14), repeat=4):
        if p % 2 != r % 2 or q % 2 != s % 2:
            continue
        m1, m2, m3, m4 = (MS[x // 2] for x in (p, q, r, s))
        if m1 + m2 != m3 + m4:
            continue
        v = sum(ck(k, m1, m3) * ck(k, m4, m2) * Fk[k] for k in (2, 4, 6))
        if abs(v) > 1e-12:
            g[(p, q, r, s)] = v
    return g


def _ann(det, o):
    if not det >> o & 1:
        return None, 0
    return det ^ (1 << o), (-1) ** bin(det & ((1 << o) - 1)).count('1')


def _cre(det, o):
    if det >> o & 1:
        return None, 0
    return det | (1 << o), (-1) ** bin(det & ((1 << o) - 1)).count('1')


@lru_cache(None)
def dets(n):
    ds = [sum(1 << o for o in c) for c in itertools.combinations(range(14), n)]
    return ds, {d: i for i, d in enumerate(ds)}


def one_body(h, n):
    ds, ix = dets(n); N = len(ds); H = np.zeros((N, N), complex)
    nz = [(p, q, h[p, q]) for p in range(14) for q in range(14) if abs(h[p, q]) > 1e-12]
    for j, d in enumerate(ds):
        for p, q, v in nz:
            d1, s1 = _ann(d, q)
            if d1 is None: continue
            d2, s2 = _cre(d1, p)
            if d2 is None: continue
            H[ix[d2], j] += v * s1 * s2
    return H


def two_body(F2, F4, F6, n):
    ds, ix = dets(n); N = len(ds); H = np.zeros((N, N))
    g = twoel_table(F2, F4, F6)
    # antisymmetrised, p<q, r<s
    anti = {}
    for p in range(14):
        for q in range(p + 1, 14):
            for r in range(14):
                for s in range(r + 1, 14):
                    v = g.get((p, q, r, s), 0.0) - g.get((p, q, s, r), 0.0)
                    if abs(v) > 1e-12:
                        anti.setdefault((r, s), []).append((p, q, v))
    for j, d in enumerate(ds):
        occ = [o for o in range(14) if d >> o & 1]
        for r, s in itertools.combinations(occ, 2):
            if (r, s) not in anti: continue
            d1, s1 = _ann(d, r); d2, s2 = _ann(d1, s)
            for p, q, v in anti[(r, s)]:
                d3, s3 = _cre(d2, q)
                if d3 is None: continue
                d4, s4 = _cre(d3, p)
                if d4 is None: continue
                H[ix[d4], j] += v * s1 * s2 * s3 * s4
    return H


@lru_cache(None)
def _two_body_cached(F2, F4, F6, n):
    return two_body(F2, F4, F6, n)


def lfci(Vm, zeta, F2, F4, F6, n):
    h, mu = spinorb_ops(Vm, zeta)
    H = one_body(h, n) + _two_body_cached(round(F2, 1), round(F4, 1), round(F6, 1), n)
    H = (H + H.conj().T) / 2
    E, U = np.linalg.eigh(H)
    M = [U.conj().T @ one_body(m, n) @ U for m in mu]
    return E - E[0], U, M


def chi_tensor(E, M, T, nlow=None, degtol=1e-3):
    """Van Vleck susceptibility tensor (cm3/mol) at temperature T (K)."""
    Ek = E * CM2K
    if nlow is None:
        nlow = int(np.searchsorted(Ek, Ek[0] + 40 * max(T, 1)))  # thermally relevant states
        nlow = max(nlow, 2)
    w = np.exp(-(Ek[:nlow] - Ek[0]) / T); Z = w.sum()
    chi = np.zeros((3, 3))
    dE = Ek[:nlow, None] - Ek[None, :]                       # E_n - E_m
    deg = np.abs(dE) < degtol * CM2K
    inv = np.where(deg, 0.0, 1.0 / np.where(deg, 1.0, dE))
    for a in range(3):
        for b in range(a, 3):
            P = (M[a][:nlow, :] * M[b][:, :nlow].T).real        # Re <n|Ma|m><m|Mb|n>
            term = np.where(deg, P / T, -2 * P * inv)
            chi[a, b] = chi[b, a] = (w * term.sum(1)).sum() / Z
    return NAMUB2_KB * chi


def doublet_g(M, i0=0):
    """g tensor of the doublet (states i0, i0+1): G = g g^T = 2 Tr(M_a M_b) over the doublet."""
    Md = [m[i0:i0 + 2, i0:i0 + 2] for m in M]
    G = np.array([[2 * np.trace(Md[a] @ Md[b]).real for b in range(3)] for a in range(3)])
    w, U = np.linalg.eigh(G)
    return np.sqrt(np.clip(w, 0, None)), U


# ----------------------------------------------------------------- convenience
def load(path, block='CASSCF'):
    d = json.load(open(path))
    a = d['ailft'][block]
    Lsoc = d['Lsoc']
    if not all(a in Lsoc for a in 'XYZ'):
        # spin-free runs print no LF-SOC integrals; ORCA's real-f basis (f0, f+1, f-1, ...) is the same in every AILFT
        # output, so the angular-momentum matrices are taken from a reference run (only their structure is used).
        import os
        ref = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'parsed', 'R10_Nd_bulk_R25_FeNeutral.json')
        Lsoc = json.load(open(ref))['Lsoc']
    U = real_to_m(Lsoc)
    V = np.array(a['V']); V = (V + V.T) / 2
    Vm = U.conj().T @ V @ U
    Vm -= np.trace(Vm).real / 7 * np.eye(7)
    return d, a, Vm, U
