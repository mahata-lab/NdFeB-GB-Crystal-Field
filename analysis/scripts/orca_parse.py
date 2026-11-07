#!/usr/bin/env python3
"""Extract the AILFT / SOC / Kramers-pair data from an ORCA 6 SA-CASSCF(+NEVPT2)+QDPT+AILFT output.

python3 orca_parse.py <orca.out> <label> <ion> <out.json>
Stores: AILFT one-electron matrix (cm-1, ORCA real-f basis), LF eigenvalues and eigenvectors,
Slater F2,F4,F6 and zeta (cm-1) for every AILFT block (CASSCF and, if present, NEVPT2),
LF-SOC-X/Y/Z integrals, QDPT-SOC eigenvalues (first block = CASSCF, second = NEVPT2 if present),
and the Kramers-pair-1 <S>,<L> matrix elements."""
import re, sys, json
HA2CM = 219474.6313705
NUM = r'-?\d+\.\d+'

def read_matrix(L, i, n=7):
    M = [[None]*n for _ in range(n)]
    while i < len(L):
        cols = L[i].split()
        if not cols or not all(c.isdigit() for c in cols): break
        idx = [int(c) for c in cols]; i += 1
        for r in range(n):
            p = L[i].split(); i += 1
            for k, c in enumerate(idx): M[int(p[0])][c] = float(p[1+k])
    return M, i

def parse(path):
    txt = open(path, errors='replace').read(); L = txt.splitlines()
    out = dict(normal='ORCA TERMINATED NORMALLY' in txt, ailft={}, soc=[], kramers=None, Lsoc={})
    for i, l in enumerate(L):
        s = l.strip()
        m = re.match(r'AILFT MATRIX ELEMENTS \((\w+)\)', s)
        if m:
            lab = m.group(1); j = i+1
            while not L[j].strip().startswith('Orbital'): j += 1
            V = [[float(x)*HA2CM for x in L[j+1+r].split()[1:8]] for r in range(7)]
            blk = dict(V=V)
            k = j
            while 'one electron eigenfunctions' not in L[k]:
                for key in ('F2ff', 'F4ff', 'F6ff'):
                    if L[k].startswith(key) and key not in blk:
                        blk[key] = float(L[k].split('=')[-1].split()[0])
                k += 1
            k += 3; E = []; C = []
            for r in range(7):
                p = L[k+r].split(); E.append(float(p[2])); C.append([float(x) for x in p[3:10]])
            blk['E'] = E; blk['C'] = C
            for kk in range(k, min(k+600, len(L))):
                mm = re.search(r'SOC constant zeta\s*=\s*' + NUM + r' eV =\s*(' + NUM + ')', L[kk])
                if mm: blk['zeta'] = float(mm.group(1)); break
            if lab not in out['ailft']: out['ailft'][lab] = blk
        if s.startswith('LF-SOC-') and s.endswith('integrals (cm-1)'):
            ax = s[7]
            if ax not in out['Lsoc']:
                M, _ = read_matrix(L, i+1); out['Lsoc'][ax] = M
        if l.startswith('Lowest eigenvalue of the SOC matrix'):
            j = i+1
            while not L[j].startswith('Eigenvalues:'): j += 1
            j += 1; ev = []
            while j < len(L):
                mm = re.match(r'\s*(\d+):\s+(' + NUM + ')', L[j])
                if not mm: break
                ev.append(float(mm.group(2))); j += 1
            out['soc'].append(ev[:200])
        if s == 'KRAMERS PAIR 1' and out['kramers'] is None:
            d = {}
            for j in range(i+1, i+12):
                mm = re.match(r'\s*Matrix elements (Re|Im)<1\|([SL])\|([12])>\s+(' + NUM + r')\s+(' + NUM + r')\s+(' + NUM + ')', L[j])
                if mm: d['%s_%s_%s' % mm.groups()[:3]] = [float(mm.group(x)) for x in (4, 5, 6)]
            out['kramers'] = d
    return out

if __name__ == '__main__':
    p, label, ion, dst = sys.argv[1:5]
    o = parse(p); o.update(label=label, ion=ion, source=p)
    json.dump(o, open(dst, 'w'))
    print(label, 'normal=%s' % o['normal'], 'ailft=%s' % list(o['ailft']), 'soc_blocks=%d' % len(o['soc']))
