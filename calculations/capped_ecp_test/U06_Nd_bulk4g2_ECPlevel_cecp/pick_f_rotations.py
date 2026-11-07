#!/usr/bin/env python3
"""pick_f_rotations.py step0a.out nel norb out_rotate.txt
Reads the Loewdin reduced orbital populations per MO of the ROHF run, takes the 7 lowest MOs with >= 50% f character on
atom 0 (the lanthanide), and writes an ORCA %scf rotate block that swaps them into the CAS window
(window = inactive..inactive+norb-1, inactive = (N_el - nel)/2)."""
import re, sys
out, nel, norb, dst = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
L = open(out, errors='replace').read().splitlines()
nelec = None
for l in L:
    m = re.search(r'Number of Electrons\s+NEL\s+\.\.\.\.\s+(\d+)', l)
    if m: nelec = int(m.group(1))
start = max(i for i, l in enumerate(L) if 'REDUCED ORBITAL POPULATIONS PER MO' in l and 'LOEWDIN' in l)
fpct = {}; i = start
while i < len(L) and not L[i].startswith('*****') and 'MAYER' not in L[i] and i < start + 20000:
    cols = L[i].split()
    if cols and all(c.isdigit() for c in cols) and len(cols) <= 10:
        idx = [int(c) for c in cols]; i += 3
        if set(L[i].strip()) <= {'-', ' '}: i += 1
        while i < len(L) and L[i].strip():
            m = re.match(r'\s*(\d+)\s+(\w+)\s+([spdfgh])\S*\s+(.*)$', L[i])
            if m and m.group(1) == '0' and m.group(3) == 'f':
                vals = [float(v) for v in m.group(4).split()]
                for k, v in zip(idx, vals): fpct[k] = fpct.get(k, 0) + v
            i += 1
    i += 1
inact = (nelec - nel) // 2
window = list(range(inact, inact + norb))
best = [k for k in sorted(fpct) if fpct[k] >= 50.0][:norb]   # 4f = the lowest-lying MOs of mainly f character (5f-like virtuals lie far above)
print('N_el', nelec, 'window', window, 'top-f MOs', sorted(best), ['%d:%.0f%%' % (k, fpct[k]) for k in sorted(best)])
out_of = [k for k in best if k not in window]; free = [k for k in window if k not in best]
rot = ['   {%d, %d, 90}' % (a, b) for a, b in zip(sorted(out_of), sorted(free))]
open(dst, 'w').write('%scf\n  rotate\n' + '\n'.join(rot) + '\n  end\nend\n' if rot else '')
print('rotations:', rot if rot else 'none needed')
