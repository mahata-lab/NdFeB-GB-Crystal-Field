#!/usr/bin/env python3
"""build_inputs.py - regenerate every ORCA job folder of this repository from the two relaxed structures.

    python3 build_inputs.py OUTDIR            # write all job folders into OUTDIR/<group>/<job>/
    python3 build_inputs.py OUTDIR R10 U04    # only jobs whose name starts with R10 or U04
    python3 build_inputs.py --check           # regenerate into a temporary folder and compare with ../calculations

Each job folder contains
    sphere.pc        point charges (ORCA format: count, then q x y z in Angstrom)
    step*.inp        ORCA 6.1 inputs (see "Protocols" below)
    run.sh           driver; needs ORCA_BIN = full path of the orca binary
    provenance.json  how the embedding was built (centre, radius, charges, neutralisation, nearest charge)

Embedding (all jobs except the single-cell tests, which are kept as published in ../calculations/single_cell_*):
    the periodic images of the relaxed cell within a sphere of radius R around the rare-earth ion become point charges,
    RE +3, Fe -3/14, B -3 (formula-unit neutral; alternative partition 'FeBalt': B -1, Fe -5/14); the residual charge of
    the truncated sphere is spread equally over all charges with R - 2 A < r <= R, so the array is exactly neutral.
Protocols:
    DKH     production: DKH2, SARC2-DKH-QZVP, RIJCOSX/AutoAux; step0 free ion -> step1 ion in field -> step2 SA-CASSCF
            + QDPT-SOC + AILFT + g tensor (ActConstraints 1, actorbs forbs). Options: SC-NEVPT2 (T13-T16), CASCI on the
            free-ion orbitals (T18-T22).
    ECP     capped-ECP test (U01-U16): def2-TZVP with the def2 small-core ECP on the ion, spin-free SA-CASSCF + AILFT;
            arm 'cecp' replaces every Fe and rare-earth neighbour within 4.0 A by a capped ECP with the same charge
            (Fe: Hay-Wadt, RE: def2-ECP) and starts from a ROHF in the same environment (see run.sh, pick_f_rotations.py).
Site numbering (one-based, as in ../structures/*_CONTCAR): bulk atoms 1-4 = 4g (x,-x,0), 5-8 = 4f (x,x,0);
GB classes A {5,8,13,16}, B {1,4,9,12}, C {2,3,10,11}, D {6,7,14,15}.
"""
import os, sys, json, re, shutil, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
BULK_CONTCAR = os.path.join(ROOT, 'structures', 'bulk_Nd8Fe56B4_CONTCAR')
GB_CONTCAR = os.path.join(ROOT, 'structures', 'GB90twist_Nd16Fe112B8_CONTCAR')

Q_PROD = {'Nd': 3.0, 'Fe': -3 / 14, 'B': -3.0, 'Zr': 4.0}
Q_ALT = {'Nd': 3.0, 'Fe': -5 / 14, 'B': -1.0}
ION = {'Ce': dict(nel=1, m1='2', r1='7',  m2='2',   r2='7',      mult=2),
       'Pr': dict(nel=2, m1='3', r1='21', m2='3,1', r2='21,28',  mult=3),
       'Nd': dict(nel=3, m1='4', r1='35', m2='4,2', r2='35,112', mult=4),
       'Gd': dict(nel=7, m1='8', r1='1',  m2='8,6', r2='1,48',   mult=8),
       'Tb': dict(nel=8, m1='7', r1='7',  m2='7,5', r2='7,140',  mult=7),
       'Dy': dict(nel=9, m1='6', r1='21', m2='6,4', r2='21,224', mult=6)}
CAP_R = 4.0


# ----------------------------------------------------------------------------------------------------- structure
def contcar(path):
    L = open(path).read().split('\n'); s = float(L[1])
    lat = np.array([[float(x) for x in L[i].split()] for i in (2, 3, 4)]) * s
    syms = L[5].split(); cnt = [int(x) for x in L[6].split()]; i = 7
    if L[i].strip()[0] in 'Ss': i += 1
    mode = L[i].strip()[0].lower(); i += 1
    pos = np.array([[float(x) for x in L[i + k].split()[:3]] for k in range(sum(cnt))])
    if mode == 'd': pos = pos @ lat
    return lat, sum([[a] * c for a, c in zip(syms, cnt)], []), pos


def sphere(lat, el, pos, idx, q, R, shell=2.0):
    """Point charges of all periodic images within R of atom idx; residual charge spread over the outer shell."""
    c = pos[idx]; rows = []
    n = int(np.ceil(R / min(np.linalg.norm(lat, axis=1)))) + 1
    for a in range(-n, n + 1):
        for b in range(-n, n + 1):
            for d in range(-n, n + 1):
                sh = np.array([a, b, d]) @ lat
                for e, p in zip(el, pos):
                    r = p + sh; dist = np.linalg.norm(r - c)
                    if 1e-6 < dist <= R: rows.append([e, q[e], r[0], r[1], r[2], dist])
    tot = sum(x[1] for x in rows); outer = [x for x in rows if x[5] > R - shell]
    for x in outer: x[1] -= tot / len(outer)
    return c, rows, tot, len(outer)


def nearest(lat, el, pos, idx, elem):
    """Nearest atom of type elem to atom idx (minimum image); returns (index, distance)."""
    best = None
    for j, (e, p) in enumerate(zip(el, pos)):
        if e != elem or j == idx: continue
        for a in (-1, 0, 1):
            for b in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    dd = np.linalg.norm(p + np.array([a, b, dz]) @ lat - pos[idx])
                    if best is None or dd < best[1]: best = (j, dd)
    return best


def w(path, text):
    with open(path, 'w', newline='\n') as fh: fh.write(text)


def write_pc(path, rows):
    with open(path, 'w', newline='\n') as fh:
        fh.write('%d\n' % len(rows))
        for e, q, x, y, z, d in rows: fh.write('%10.6f  %12.6f %12.6f %12.6f\n' % (q, x, y, z))


# ----------------------------------------------------------------------------------------------------- run scripts
RUN_DKH = """#!/bin/bash
# step0 free ion -> step1 ion in the point-charge field -> step2 SA-CASSCF + SOC + AILFT (result: orca.out)
# step0 is identical for every job of a given ion; if step0.gbw is supplied it is reused.
set -e
: "${ORCA_BIN:?set ORCA_BIN to the full path of the orca binary}"
if [ ! -f step0.gbw ]; then
  "$ORCA_BIN" step0.inp > step0.out 2>&1 || true
  grep -q "ORCA TERMINATED NORMALLY" step0.out || { echo STEP0_FAIL; tail -25 step0.out; exit 1; }
fi
echo STEP0_OK
"$ORCA_BIN" step1.inp > step1.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" step1.out || { echo STEP1_FAIL; tail -25 step1.out; exit 1; }
echo STEP1_OK
"$ORCA_BIN" step2.inp > orca.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" orca.out || { echo STEP2_FAIL; tail -25 orca.out; exit 1; }
grep -q "active space contains 7 f Orbitals : OK" orca.out && echo "STEP2_OK (7 pure 4f active orbitals)" || echo "STEP2_WARN: check active space"
"""
RUN_CASCI = """#!/bin/bash
# Ce CASCI variant: step0 free ion, then step2 = SA-CASCI (NoIter) on the free-ion 4f orbitals + SOC + AILFT (orca.out)
set -e
: "${ORCA_BIN:?set ORCA_BIN to the full path of the orca binary}"
if [ ! -f step0.gbw ]; then
  "$ORCA_BIN" step0.inp > step0.out 2>&1 || true
  grep -q "ORCA TERMINATED NORMALLY" step0.out || { echo STEP0_FAIL; tail -25 step0.out; exit 1; }
fi
"$ORCA_BIN" step2.inp > orca.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" orca.out && echo STEP2_OK || { echo STEP2_FAIL; tail -25 orca.out; exit 1; }
"""
RUN_GUESS = """#!/bin/bash
# step1 starts from supplied orbitals (see NOTE.md), then step2 (orca.out)
set -e
: "${ORCA_BIN:?set ORCA_BIN to the full path of the orca binary}"
"$ORCA_BIN" step1.inp > step1.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" step1.out || { echo STEP1_FAIL; tail -25 step1.out; exit 1; }
"$ORCA_BIN" step2.inp > orca.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" orca.out && echo STEP2_OK || { echo STEP2_FAIL; tail -25 orca.out; exit 1; }
"""
RUN_ECP_PC = """#!/bin/bash
# ECP-level spin-free test, bare point charges: step0 free ion -> step1 in field -> step2 SA-CASSCF + AILFT (orca.out)
set -e
: "${ORCA_BIN:?set ORCA_BIN to the full path of the orca binary}"
for s in step0 step1; do
  "$ORCA_BIN" $s.inp > $s.out 2>&1 || true
  grep -q "ORCA TERMINATED NORMALLY" $s.out || { echo ${s}_FAIL; tail -25 $s.out; exit 1; }
done
"$ORCA_BIN" step2.inp > orca.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" orca.out && echo STEP2_OK || { echo STEP2_FAIL; tail -25 orca.out; exit 1; }
"""
RUN_ECP_CAP = """#!/bin/bash
# ECP-level spin-free test, capped ECPs on the Fe and rare-earth neighbours within 4.0 A (orca.out).
# step0a: ROHF on the same atom list (ion + uncharged caps) with Loewdin MO populations;
# pick_f_rotations.py rotates the seven lowest MOs of mainly 4f character into the CAS window (if needed);
# step0: SA-CASSCF from those orbitals; step1/step2 in the field as for the point-charge arm.
set -e
: "${ORCA_BIN:?set ORCA_BIN to the full path of the orca binary}"
HERE="$(cd "$(dirname "$0")" && pwd)"; PICK="$HERE/pick_f_rotations.py"; [ -f "$PICK" ] || PICK="$HERE/../../../embedding/pick_f_rotations.py"
"$ORCA_BIN" step0a.inp > step0a.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" step0a.out || { echo step0a_FAIL; tail -20 step0a.out; exit 1; }
python3 "$PICK" step0a.out NEL 7 rotate.txt
if [ -s rotate.txt ]; then
  awk 'NR==FNR{r=r $0 "\\n"; next} /^%maxcore/{print; printf "%s", r; next} {print}' rotate.txt step0_rohf.inp > step0.inp
else
  cp step0_rohf.inp step0.inp
fi
for s in step0 step1; do
  "$ORCA_BIN" $s.inp > $s.out 2>&1 || true
  grep -q "ORCA TERMINATED NORMALLY" $s.out || { echo ${s}_FAIL; tail -25 $s.out; exit 1; }
done
"$ORCA_BIN" step2.inp > orca.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" orca.out && echo STEP2_OK || { echo STEP2_FAIL; tail -25 orca.out; exit 1; }
"""


# ----------------------------------------------------------------------------------------------------- DKH inputs
def dkh_inputs(folder, ion, c, nevpt2=False, casci=False, guess=None):
    I = ION[ion]; pc = 'sphere.pc'
    xyz = "* xyz 3 %d\n%s   %.8f %.8f %.8f\n*\n" % (I['mult'], ion, *c)
    s0 = ("# STEP0 : free %s3+ ion, no point charges -> clean atomic 4f\n! DKH2 SARC2-DKH-QZVP AutoAux CASSCF\n%%maxcore 14000\n%%pal nprocs 16 end\n%%scf maxiter 1000 end\n"
          "%%casscf\n  nel  %d\n  norb 7\n  mult %s\n  nroots %s\n  maxiter 300\n  orbstep SuperCI\n  switchstep NR\n  shiftup 1.0\n  shiftdn 1.0\n  minshift 0.2\n  actorbs forbs\nend\n"
          "* xyz 3 %d\n%s   0.0 0.0 0.0\n*\n") % (ion, I['nel'], I['m1'], I['r1'], I['mult'], ion)
    head1 = guess[0] if guess else "# STEP1 : %s3+ in neutral spherical point-charge field, reads step0.gbw" % ion
    s1 = ("%s\n! DKH2 RIJCOSX SARC2-DKH-QZVP AutoAux CASSCF MORead\n%%moinp \"%s\"\n%%maxcore 14000\n%%pal nprocs 16 end\n%%pointcharges \"%s\"\n%%scf maxiter 800 end\n"
          "%%casscf\n  ActConstraints 1\n  nel  %d\n  norb 7\n  mult %s\n  nroots %s\n  maxiter 800\n  orbstep    SuperCI\n  switchstep DIIS\n  shiftup 1.0\n  shiftdn 1.0\n  minshift 0.2\n  actorbs forbs\nend\n"
          ) % (head1, guess[1] if guess else 'step0.gbw', pc, I['nel'], I['m1'], I['r1']) + xyz
    if casci:
        head2, kw, moinp = "# STEP2 : SA-CASCI (NoIter, free-ion 4f orbitals) + SOC + AILFT + g-tensor", " NoIter", "step0.gbw"
    elif nevpt2:
        head2, kw, moinp = "# STEP2 : SA-CASSCF + SC-NEVPT2 + SOC + AILFT + g-tensor, reads step1.gbw", "", "step1.gbw"
    else:
        head2, kw, moinp = "# STEP2 : SA-CASSCF + SOC + AILFT + g-tensor, reads step1.gbw", "", "step1.gbw"
    s2 = ("%s\n! DKH2 RIJCOSX SARC2-DKH-QZVP AutoAux CASSCF MORead%s\n%%moinp \"%s\"\n%%maxcore 14000\n%%pal nprocs 16 end\n%%pointcharges \"%s\"\n%%casscf\n%s"
          "  ActConstraints 1\n  nel  %d\n  norb 7\n  mult %s\n  nroots %s\n  maxiter 200\n  orbstep    SuperCI\n  switchstep DIIS\n  shiftup 0.5\n  shiftdn 0.5\n  minshift 0.1\n"
          "  rel\n    dosoc true\n    gtensor true\n    printlevel 3\n  end\n  actorbs forbs\nend\n"
          ) % (head2, kw, moinp, pc, "  PTMethod SC_NEVPT2\n" if nevpt2 else "", I['nel'], I['m2'], I['r2']) + xyz
    files = {'step0.inp': s0, 'step2.inp': s2}
    if not casci: files['step1.inp'] = s1
    for n, t in files.items(): w(os.path.join(folder, n), t)
    w(os.path.join(folder, 'run.sh'), RUN_CASCI if casci else (RUN_GUESS if guess else RUN_DKH))


# ----------------------------------------------------------------------------------------------------- ECP inputs
def ecp_inputs(folder, ion, c, caps):
    I = ION[ion]
    capblk = ''.join('%s> %10.6f  %12.6f %12.6f %12.6f  NewECP "%s" end\n' % (e, q, x, y, z, 'HayWadt' if e == 'Fe' else 'def2-ECP')
                     for e, q, x, y, z, d in caps)
    xyz = "* xyz 3 %d\n%s   %.8f %.8f %.8f\n%s*\n" % (I['mult'], ion, *c, capblk)
    common = "! RIJCOSX def2-TZVP AutoAux CASSCF%s\n%s%%maxcore 3000\n%%pal nprocs 12 end\n"
    cas = "%%casscf\n%s  nel  %d\n  norb 7\n  mult %s\n  nroots %s\n  maxiter %d\n  orbstep SuperCI\n  switchstep DIIS\n  shiftup %s\n  shiftdn %s\n  minshift %s\n  actorbs forbs\nend\n"
    s1 = ("# STEP1 : %s3+ in field\n" % ion + common % (' MORead', '%moinp "step0.gbw"\n') + '%pointcharges "sphere.pc"\n'
          + cas % ('  ActConstraints 1\n', I['nel'], I['m1'], I['r1'], 800, '1.0', '1.0', '0.2') + xyz)
    s2 = ("# STEP2 : spin-free SA-CASSCF + AILFT\n" + common % (' MORead', '%moinp "step1.gbw"\n') + '%pointcharges "sphere.pc"\n'
          + cas % ('  ActConstraints 1\n', I['nel'], I['m2'], I['r2'], 200, '0.5', '0.5', '0.1') + xyz)
    files = {'step1.inp': s1, 'step2.inp': s2}
    if not caps:
        files['step0.inp'] = ("# STEP0 : free %s3+ ion, def2-TZVP/def2-ECP\n" % ion + common % ('', '')
                              + cas % ('', I['nel'], I['m1'], I['r1'], 300, '1.0', '1.0', '0.2') + "* xyz 3 %d\n%s   0.0 0.0 0.0\n*\n" % (I['mult'], ion))
        run = RUN_ECP_PC
    else:   # step0 must contain the same atom list as step1/2 (MORead); caps carry no charge there, no point charges
        xyz0 = re.sub(r'^(\w+>)\s+-?\d+\.\d+', lambda m: m.group(1) + '   0.000000', xyz, flags=re.M)
        files['step0a.inp'] = ("# STEP0a : ROHF high-spin orbitals on the same atom list (ion + uncharged capped-ECP centres), guess for step0\n"
                               "! ROHF RIJCOSX def2-TZVP AutoAux SlowConv\n%maxcore 3000\n%pal nprocs 12 end\n%scf maxiter 500 end\n"
                               "%output\n  Print[ P_OrbPopMO_L ] 1\n  Print[ P_ReducedOrbPopMO_L ] 1\nend\n" + xyz0)
        files['step0_rohf.inp'] = ("# STEP0 : guess = ion + uncharged capped-ECP centres (same atom list as step1/2, no point charges); %s3+ ion, def2-TZVP/def2-ECP\n" % ion
                                   + common % (' MORead', '%moinp "step0a.gbw"\n')
                                   + cas % ('  ActConstraints 1\n', I['nel'], I['m1'], I['r1'], 300, '1.0', '1.0', '0.2') + xyz0)
        run = RUN_ECP_CAP.replace('NEL', str(I['nel']))
    for n, t in files.items(): w(os.path.join(folder, n), t)
    w(os.path.join(folder, 'run.sh'), run)


# ----------------------------------------------------------------------------------------------------- catalogue
def catalogue():
    """(group, job, builder kwargs) for every regenerable job of the paper."""
    J = []; ions = ['Ce', 'Pr', 'Nd', 'Gd', 'Tb', 'Dy']
    k = 0
    for ion in ions:                                                      # R01-R24: radius test, bulk 4g / GB class A
        for env, atom in (('bulk', 2), ('GB', 13)):
            for R in (18, 25):
                k += 1; J.append(('converged_sphere', 'R%02d_%s_%s_R%d_FeNeutral' % (k, ion, env, R), dict(kind='dkh', ion=ion, env=env, atom=atom, R=R)))
    for ion in ('Nd', 'Dy'):                                              # R25-R28: alternative charge partition
        for env, atom in (('bulk', 2), ('GB', 13)):
            k += 1; J.append(('converged_sphere', 'R%02d_%s_%s_R25_FeBalt' % (k, ion, env), dict(kind='dkh', ion=ion, env=env, atom=atom, R=25, alt=True)))
    k = 0
    for env, atom, lab in (('bulk', 5, 'bulk4f_atom5'), ('GB', 9, 'GBclassB_atom9')):   # T01-T12: matched sublattice sites
        for ion in ions:
            k += 1; J.append(('matched_sites', 'T%02d_%s_%s_R25' % (k, ion, lab), dict(kind='dkh', ion=ion, env=env, atom=atom, R=25)))
    for k, (ion, env, atom, lab) in enumerate((('Nd', 'bulk', 2, 'bulk_atom2'), ('Nd', 'GB', 13, 'GB_atom13'),
                                               ('Dy', 'bulk', 2, 'bulk_atom2'), ('Dy', 'GB', 13, 'GB_atom13')), 13):
        J.append(('nevpt2_production', 'T%02d_%s_%s_R25_NEVPT2' % (k, ion, lab), dict(kind='dkh', ion=ion, env=env, atom=atom, R=25, nevpt2=True)))
    J.append(('ce_active_space_tests', 'T17_Ce_GB_atom13_R25_retry', dict(kind='dkh', ion='Ce', env='GB', atom=13, R=25,
              guess=("# STEP1 : Ce3+ GB atom 13 in 25 A sphere; guess = converged 4f orbitals of R03 (same site, 18 A sphere)", "guess_R18_step2.gbw"))))
    for job, env, atom, R in (('T18_Ce_bulk4f_atom5_R25_CASCI', 'bulk', 5, 25), ('T19_Ce_GB_atom13_R25_CASCI', 'GB', 13, 25),
                              ('T20_Ce_bulk_atom2_R25_CASCI_check', 'bulk', 2, 25), ('T21_Ce_GBclassB_atom9_R25_CASCI_check', 'GB', 9, 25),
                              ('T22_Ce_GB_atom13_R18_CASCI_check', 'GB', 13, 18)):
        J.append(('ce_active_space_tests', job, dict(kind='dkh', ion='Ce', env=env, atom=atom, R=R, casci=True)))
    k = 0
    for ion in ('Nd', 'Dy'):                                              # U01-U16: capped-ECP test
        for site, env, atom in (('bulk4f5', 'bulk', 5), ('GBA13', 'GB', 13), ('bulk4g2', 'bulk', 2), ('GBB9', 'GB', 9)):
            for arm in ('pc', 'cecp'):
                k += 1; J.append(('capped_ecp_test', 'U%02d_%s_%s_ECPlevel_%s' % (k, ion, site, arm), dict(kind='ecp', ion=ion, env=env, atom=atom, R=25, site=site, cap=arm == 'cecp')))
    k = 0
    for ion in ('Nd', 'Dy'):                                              # V01-V08: B / Zr substitution at the boundary plane
        for site, atom in (('GBA13', 13), ('GBB9', 9)):
            for sub in ('Bsub', 'Zrsub'):
                k += 1; J.append(('boundary_chemistry', 'V%02d_%s_%s_%s' % (k, ion, site, sub), dict(kind='dkh', ion=ion, env='GB', atom=atom, R=25, site=site, sub=sub)))
    k = 0
    for ion in ('Nd', 'Dy'):                                              # W01-W04: mid-grain classes C and D
        for site, atom in (('GBC2', 2), ('GBD6', 6)):
            k += 1; J.append(('midgrain_sites', 'W%02d_%s_%s' % (k, ion, site), dict(kind='dkh', ion=ion, env='GB', atom=atom, R=25, site=site)))
    return J


def build(outdir, sel=None):
    S = {'bulk': contcar(BULK_CONTCAR), 'GB': contcar(GB_CONTCAR)}
    done = []
    for group, job, a in catalogue():
        if sel and not any(job.startswith(s) for s in sel): continue
        folder = os.path.join(outdir, group, job); os.makedirs(folder, exist_ok=True)
        lat, el, pos = S[a['env']]; el = list(el); idx = a['atom'] - 1
        q = dict(Q_ALT if a.get('alt') else Q_PROD); prov = dict(job=job, ion=a['ion'], environment=a['env'], atom_1based=a['atom'], R_A=a['R'])
        if a.get('sub'):
            j, dist = nearest(lat, el, pos, idx, 'Fe' if a['sub'] == 'Bsub' else 'Nd')
            prov.update(substitution=a['sub'], substituted_atom_1based=j + 1, substituted_element=el[j],
                        new_element='B' if a['sub'] == 'Bsub' else 'Zr', distance_A=round(dist, 4), note='unrelaxed substitution, one per periodic GB cell')
            el[j] = prov['new_element']
        c, rows, raw, nout = sphere(lat, el, pos, idx, q, a['R'])
        caps = [x for x in rows if a.get('cap') and x[5] <= CAP_R and x[0] in ('Fe', 'Nd')]
        pcs = [x for x in rows if not any(x is y for y in caps)]
        write_pc(os.path.join(folder, 'sphere.pc'), pcs)
        if a['kind'] == 'ecp':
            ecp_inputs(folder, a['ion'], c, caps)
            prov.update(level='def2-TZVP/def2-ECP, spin-free SA-CASSCF + AILFT', arm='cecp' if a.get('cap') else 'pc', n_capped=len(caps),
                        capped=[[x[0], round(x[5], 3)] for x in caps])
            if caps: shutil.copy(os.path.join(HERE, 'pick_f_rotations.py'), folder)
        else:
            dkh_inputs(folder, a['ion'], c, nevpt2=a.get('nevpt2', False), casci=a.get('casci', False), guess=a.get('guess'))
            prov.update(level='DKH2/SARC2-DKH-QZVP, SA-CASSCF%s + QDPT-SOC + AILFT' % (' + SC-NEVPT2' if a.get('nevpt2') else (' (CASCI, NoIter)' if a.get('casci') else '')))
        d = np.array([x[5] for x in pcs])
        prov.update(center_A=[round(float(v), 8) for v in c], charges=dict(RE=q['Nd'], Fe=round(q['Fe'], 6), B=q['B'], **({'Zr': 4.0} if a.get('sub') == 'Zrsub' else {})),
                    n_pc=len(pcs), net_charge=round(sum(x[1] for x in rows), 6), residual_before_neutralisation=round(raw, 3),
                    outer_shell_charges=nout, outer_shell_correction_e=round(-raw / nout, 4), min_pc_dist_A=round(float(d.min()), 4))
        json.dump(prov, open(os.path.join(folder, 'provenance.json'), 'w', newline='\n'), indent=1)
        done.append(os.path.join(group, job))
    return done


# ----------------------------------------------------------------------------------------------------- check
def check(ref=os.path.join(ROOT, 'calculations')):
    tmp = tempfile.mkdtemp(prefix='nfb_check_'); jobs = build(tmp); bad = 0
    skip_line = re.compile(r'^(#|%pal|%maxcore)')
    for j in jobs:
        a, b = os.path.join(tmp, j), os.path.join(ref, j)
        if not os.path.isdir(b): print('missing in repository:', j); bad += 1; continue
        msg = []
        P, Pr = np.loadtxt(os.path.join(a, 'sphere.pc'), skiprows=1), np.loadtxt(os.path.join(b, 'sphere.pc'), skiprows=1)
        if P.shape != Pr.shape or np.abs(P - Pr).max() > 2e-6: msg.append('sphere.pc differs')
        for f in sorted(os.listdir(a)):
            if not f.endswith('.inp'): continue
            if not os.path.exists(os.path.join(b, f)): msg.append(f + ' missing'); continue
            ta = [l for l in open(os.path.join(a, f)).read().splitlines() if not skip_line.match(l)]
            tb = [l for l in open(os.path.join(b, f)).read().splitlines() if not skip_line.match(l)]
            if ta != tb: msg.append(f + ' differs')
        print('%-48s %s' % (j, 'identical' if not msg else '; '.join(msg))); bad += bool(msg)
    shutil.rmtree(tmp)
    print('\n%d of %d job folders reproduced exactly (comment, %%pal and %%maxcore lines ignored)' % (len(jobs) - bad, len(jobs)))
    return bad


if __name__ == '__main__':
    if '--check' in sys.argv: sys.exit(1 if check() else 0)
    if len(sys.argv) < 2: print(__doc__); sys.exit(1)
    out = build(sys.argv[1], sys.argv[2:] or None); print('wrote %d job folders into %s' % (len(out), sys.argv[1]))
