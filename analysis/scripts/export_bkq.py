#!/usr/bin/env python3
"""Full Wybourne one-electron crystal-field parameters B^k_q (cm-1, cell frame: z = c, x = a) for every run,
plus Stevens-normalised A_k^0<r^k>. These are the transferable single-ion parameters (one-electron 4f operator
V = sum_kq B^k_q C^k_q). Complex B^k_q for q<0 follow from B^k_{-q} = (-1)^q conj(B^k_q).
python3 export_bkq.py parsed out.csv"""
import sys, os, glob, csv, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import cflib as c
rows = []
for f in sorted(glob.glob(os.path.join(sys.argv[1], '*.json'))):
    d = json.load(open(f))
    for block in d['ailft']:
        _, a, Vm, U = c.load(f, block); B, res = c.wybourne(Vm)
        r = dict(run=d['label'], ion=d['ion'], block=block, zeta=a.get('zeta'), F2=a.get('F2ff'), F4=a.get('F4ff'), F6=a.get('F6ff'), resid=res)
        for k in (2, 4, 6):
            for q in range(0, k + 1):
                r['Re_B%d%d' % (k, q)] = B[(k, q)].real
                if q: r['Im_B%d%d' % (k, q)] = B[(k, q)].imag
        rows.append(r)
with open(sys.argv[2], 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader()
    for r in rows: w.writerow({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
print(len(rows), 'rows')
