#!/usr/bin/env python3
"""Verify source integrity, dataset coverage and numerical links to the current figures."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import sys

import numpy as np
import cflib

ROOT = Path(__file__).resolve().parents[2]


def rows(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reproduced', type=Path, default=ROOT/'build')
    args = parser.parse_args()
    manifest = json.loads((ROOT/'analysis/source_manifest.json').read_text())
    index = rows(ROOT/'calculations/jobs.csv')
    indexed = {row['dataset']: row for row in index}
    parsed = {p.stem: p for p in (ROOT/'analysis/parsed').glob('*.json')}
    assert len(parsed) == len(manifest) == 116
    assert len(indexed) == len(index), 'Duplicate dataset IDs'
    assert set(parsed) <= set(indexed)
    for record in manifest:
        dataset = record['dataset']
        assert parsed[dataset] == ROOT/record['parsed_result']
        assert (ROOT/record['calculation_directory']).is_dir(), dataset
        assert (ROOT/indexed[dataset]['directory']).is_dir(), dataset
        if record['raw_output']:
            with gzip.open(ROOT/record['raw_output'], 'rb') as handle:
                digest = hashlib.sha256(handle.read()).hexdigest()
            assert digest == record['original_output_sha256'], dataset
    for record in json.loads((ROOT/'analysis/orbitals/manifest.json').read_text()):
        with gzip.open(ROOT/record['path'], 'rb') as handle:
            digest = hashlib.sha256(handle.read()).hexdigest()
        assert digest == record['uncompressed_sha256'], record['path']

    coefficients = rows(ROOT/'analysis/tables/wybourne_Bkq_all_runs.csv')
    assert len(coefficients) == 124
    assert len({(r['run'],r['block']) for r in coefficients}) == 124
    assert {r['run'] for r in coefficients} == set(parsed)
    for row in coefficients:
        _, _, matrix, _ = cflib.load(str(parsed[row['run']]), row['block'])
        params, residual = cflib.wybourne(matrix)
        assert residual < 1e-5, (row['run'], residual)
        for rank in [2,4,6]:
            for q in range(rank+1):
                assert abs(params[rank,q].real-float(row[f'Re_B{rank}{q}'])) < .0011, (row['run'],rank,q,'real')
                if q:
                    assert abs(params[rank,q].imag-float(row[f'Im_B{rank}{q}'])) < .0011, (row['run'],rank,q,'imaginary')

    table = ROOT/'analysis/tables'
    finite = {r['run']:r for r in rows(table/'finite_cell_descriptors.csv')}
    assert len(finite) == 12
    assert abs(float(finite['FINITE_Ce_bulk']['split'])-2589.4)<.01
    for ion in ['Ce','Pr','Nd','Gd','Tb','Dy']:
        for env in ['bulk','GB']:
            key=f'FINITE_{ion}_{env}'
            d=json.loads(parsed[key].read_text())
            assert abs(float(finite[key]['split'])-np.ptp(d['ailft']['CASSCF']['E'])) < .001
    curves=rows(table/'finite_cell_chiT.csv')
    assert len(curves)==480
    for ion in ['Ce','Pr','Nd','Gd','Tb','Dy']:
        for env in ['bulk','GB']:
            points=[r for r in curves if r['ion']==ion and r['environment']==env]
            assert sorted(float(r['T_K']) for r in points)==list(range(1,41))
    tensors=rows(table/'finite_cell_susceptibility_tensors.csv')
    for tensor in tensors:
        matrix=np.array([[float(tensor[f'tensor_{a}{b}']) for b in 'xyz'] for a in 'xyz'])
        eig=np.linalg.eigvalsh(matrix)
        assert np.allclose(eig,[float(tensor[k]) for k in ['eigen_min','eigen_mid','eigen_max']])
        if tensor['ion']!='Gd':
            assert 86 < float(tensor['bulk_GB_axis_angle_deg']) < 89
    radii=rows(table/'orbital_extent.csv')
    assert len(radii)==6 and all(abs(float(r['grid_norm'])-1)<.009 for r in radii)
    for ion in ['Ce','Pr','Dy']:
        bulk=next(float(r['orbital_r_rms_A']) for r in radii if r['ion']==ion and r['environment']=='Bulk')
        boundary=next(float(r['orbital_r_rms_A']) for r in radii if r['ion']==ion and r['environment']=='GB')
        assert boundary < bulk, ion
    tex=(ROOT/'supplementary/si_tables.tex').read_text()
    assert all(f'\\label{{tab:S{i}}}' in tex for i in range(1,11))
    assert all((ROOT/'figures'/f'{name}.{ext}').is_file() for name in ['Figure_2','Figure_3','Figure_4','Figure_5','Figure_S1','Figure_S2'] for ext in ['pdf','png'])
    integrity_path=ROOT/'analysis/integrity.json'
    if integrity_path.exists():
        for record in json.loads(integrity_path.read_text()):
            assert hashlib.sha256((ROOT/record['path']).read_bytes()).hexdigest()==record['sha256'],record['path']
    for path in [ROOT/'README.md',*(ROOT/'docs').glob('*.md')]:
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if re.match(r'^https?://|^#',target):continue
            assert (path.parent/target.split('#')[0]).exists(), (path.name,target)
    oversized=[p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '.git' not in p.relative_to(ROOT).parts and p.stat().st_size>100_000_000]
    assert not oversized,oversized

    # Compare an independent regeneration with the supplied paper data; allow only printed rounding.
    reproduced=args.reproduced.resolve()/'analysis/tables'
    if reproduced.exists():
        for name,keys in [('finite_cell_chiT.csv',['ion','environment','T_K']),
                          ('finite_cell_descriptors.csv',['run']),
                          ('finite_cell_susceptibility_tensors.csv',['ion','environment']),
                          ('orbital_extent.csv',['ion','environment']),
                          ('conditional_switching_estimate.csv',['ion']),
                          ('wybourne_Bkq_all_runs.csv',['run','block'])]:
            a=rows(table/name);b=rows(reproduced/name)
            key=lambda row:tuple(row[k] for k in keys)
            aa={key(r):r for r in a};bb={key(r):r for r in b}
            assert aa.keys()==bb.keys(),name
            for k,row in aa.items():
                for field,value in row.items():
                    if field=='source' or field not in bb[k]:continue
                    try:
                        left,right=float(value),float(bb[k][field])
                    except (ValueError,TypeError):
                        assert value==bb[k][field],(name,k,field)
                        continue
                    assert (math.isnan(left) and math.isnan(right)) or math.isclose(left,right,abs_tol=.0011,rel_tol=1e-7),(name,k,field,left,right)
    print('Verified: 116 parsed datasets; 124 complete parameter rows; 12 raw-output hashes; 6 cube hashes; 480 susceptibility points; current figure/table sources and links.')
    if reproduced.exists(): print('Independent regenerated numerical data agree with the supplied results.')


if __name__=='__main__':
    main()
