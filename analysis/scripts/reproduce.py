#!/usr/bin/env python3
"""Rebuild the current paper's data, figures and supplementary tables without ORCA."""
from pathlib import Path
import argparse
import csv
import gzip
import json
import math
import re
import subprocess
import sys

import numpy as np

import cf_tables
import orca_parse
from zero_field_response import zero_field_tensors

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
IONS = ['Ce', 'Pr', 'Nd', 'Gd', 'Tb', 'Dy']


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def call(script, *args):
    subprocess.run([sys.executable, '-X', 'utf8', str(HERE/script), *map(str, args)], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'build')
    args = parser.parse_args()
    out = args.output.resolve()
    # Never regenerate into the tracked source tree itself.
    if out == ROOT or any(out.is_relative_to(ROOT/p) for p in ['analysis', 'figures', 'supplementary', 'calculations', 'structures', 'archive', 'docs', 'embedding']):
        parser.error('Choose a separate output directory, such as build/.')
    an = out/'analysis'
    tables = an/'tables'
    tables.mkdir(parents=True, exist_ok=True)
    (an/'parsed').mkdir(exist_ok=True)
    call('cf_tables.py', ROOT/'analysis/parsed', tables/'cf_all_runs.csv')
    call('export_bkq.py', ROOT/'analysis/parsed', tables/'wybourne_Bkq_all_runs.csv')

    finite, tensors, curves, reconcile = [], [], [], []
    for ion in IONS:
        axes = {}
        for env in ['bulk', 'GB']:
            label = f'FINITE_{ion}_{env}'
            record = json.loads((ROOT/'analysis/parsed'/f'{label}.json').read_text())
            source = ROOT/record['source']
            with gzip.open(source, 'rt', errors='replace') as handle:
                text = handle.read()
            # Re-extract AILFT/SOC data from the actual raw output, then check the supplied JSON.
            fresh = orca_parse.parse(source)
            assert fresh == {k: record[k] for k in fresh}, f'Raw/JSON mismatch: {label}'
            fresh.update(label=label, ion=ion, source=record['source'])
            path = an/'parsed'/f'{label}.json'
            path.write_text(json.dumps(fresh)+'\n', encoding='utf-8')
            row = cf_tables.row(path)
            finite.append(row)
            prior = cf_tables.row(ROOT/'analysis/parsed'/f'PUB_{ion}_{env}.json')
            reconcile.append(dict(ion=ion, environment=env,
                magnetic_property_orbital_spread=row['split'],
                separate_production_orbital_spread=prior['split'],
                difference_cm1=row['split']-prior['split']))
            zero_tensors = zero_field_tensors(text)
            matrix = zero_tensors[2.]
            eigenvalues, eigenvectors = np.linalg.eigh(matrix)
            axes[env] = eigenvectors[:, -1]
            tensor = dict(ion=ion, environment=env, temperature_K=2., tensor_units='cm3 K/mol (chiT)',
                eigen_min=eigenvalues[0], eigen_mid=eigenvalues[1], eigen_max=eigenvalues[2],
                axis_x=axes[env][0], axis_y=axes[env][1], axis_z=axes[env][2], source=record['source'])
            for i, a in enumerate('xyz'):
                for j, b in enumerate('xyz'):
                    tensor[f'tensor_{a}{b}'] = matrix[i, j]
            tensors.append(tensor)
            local = [dict(ion=ion, environment=env, T_K=temp, chiT=float(np.trace(matrix)/3))
                     for temp,matrix in sorted(zero_tensors.items()) if 1<=temp<=40]
            assert len(local) == 40, (label, len(local))
            curves.extend(local)
        angle = float(np.degrees(np.arccos(np.clip(abs(axes['bulk']@axes['GB']), 0, 1))))
        for item in tensors[-2:]:
            item['bulk_GB_axis_angle_deg'] = angle
    write_csv(tables/'finite_cell_descriptors.csv', finite)
    write_csv(tables/'finite_cell_susceptibility_tensors.csv', tensors)
    write_csv(tables/'finite_cell_chiT.csv', curves)
    write_csv(tables/'finite_cell_run_reconciliation.csv', reconcile)

    # Orbital-level panels read these JSON files; all other descriptors come from the complete table.
    call('make_figures.py', an, out/'figures')
    call('make_supplementary.py', an, out/'supplementary')
    # Retain canonical names requested with the supplementary document.
    import shutil
    shutil.copy2(tables/'wybourne_Bkq_all_runs.csv', out/'supplementary/wybourne_Bkq_all_runs.csv')
    shutil.copy2(ROOT/'structures/bulk_Nd8Fe56B4_CONTCAR', out/'supplementary/bulk_CONTCAR')
    shutil.copy2(ROOT/'structures/GB90twist_Nd16Fe112B8_CONTCAR', out/'supplementary/GB_CONTCAR')
    report = dict(parsed_datasets=len(list((ROOT/'analysis/parsed').glob('*.json'))),
        susceptibility_points=len(curves), verified_raw_magnetic_outputs=12,
        figures=['Figure_1','Figure_2','Figure_3','Figure_4','Figure_5','Figure_S1','Figure_S2'],
        output=str(out), orca_invoked=False)
    (out/'reproduction-report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(f'Completed current analysis: {len(curves)} susceptibility points, 6 full cubes, 7 figure topics. Output: {out}', flush=True)


if __name__ == '__main__':
    main()
