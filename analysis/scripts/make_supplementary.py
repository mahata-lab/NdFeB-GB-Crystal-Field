from pathlib import Path
import csv,json,re,math,sys,argparse

ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description='Rebuild current numerical supplementary tables from a reproduced analysis directory.')
parser.add_argument('analysis_dir',type=Path)
parser.add_argument('output_dir',type=Path)
args=parser.parse_args()
AN=args.analysis_dir.resolve()
OUT=args.output_dir.resolve()
OUT.mkdir(parents=True,exist_ok=True)
RALL=list(csv.DictReader(open(AN/'tables/cf_all_runs.csv')))
R={r['run']:r for r in RALL if r['block']=='CASSCF'}
R.update({r['run']:r for r in csv.DictReader(open(AN/'tables/finite_cell_descriptors.csv'))})
N={r['run']:r for r in RALL if r['block']=='NEVPT2'}
IONS=['Ce','Pr','Nd','Gd','Tb','Dy']
num=lambda r,k:float(r[k])
keep=[(ROOT/'analysis/templates/structural_tables.tex').read_text()]
tex=(ROOT/'supplementary/4f-supplementary.tex').read_text()

alltables='\n\n'.join(keep)+'\n\n'
def table(caption,label,cols,header,lines):
    return '\\begin{table}[H]\n\\centering\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n\\caption{'+caption+'}\\label{'+label+'}\n\\begin{tabular}{'+cols+'}\n\\toprule\n'+header+' \\\\\n\\midrule\n'+'\n'.join(lines)+'\n\\bottomrule\n\\end{tabular}\n\\end{table}\n\n'

lines=[]
for ion in IONS:
    b=R['FINITE_'+ion+'_bulk'];g=R['FINITE_'+ion+'_GB']
    lines.append(f"{ion} & {num(b,'tun'):.2f} & {num(g,'tun'):.2f} & {num(b,'E1'):.2f} & {num(g,'E1'):.2f} & {num(b,'Jspan'):.2f} & {num(g,'Jspan'):.2f} \\\\")
alltables+=table(r'Many-electron QDPT-SOC descriptors of the finite-cell calculations (cm$^{-1}$). $\delta_{01}=E_1-E_0$ is zero for the Kramers ions; for Pr and Tb it is the lowest-pair splitting. $E_{\mathrm{next}}=E_2-E_0$ is the next Kramers doublet or the level above the lowest non-Kramers pair. $W_J$ is the span of the assigned ground-$J$ manifold. These quantities are distinct from the AILFT orbital gaps in main Table~1.','tab:S4','lrrrrrr',r'Ion & $\delta_{01}^{b}$ & $\delta_{01}^{g}$ & $E_{\mathrm{next}}^b$ & $E_{\mathrm{next}}^g$ & $W_J^b$ & $W_J^g$',lines)
lines=[]
for ion in IONS:
    for env in ['bulk','GB']:
        r=R['FINITE_'+ion+'_'+env]
        lines.append(f"{ion} & {env} & {num(r,'A20r2'):.2f} & {num(r,'S2'):.2f} & {num(r,'S4'):.2f} & {num(r,'S6'):.2f} \\\\")
alltables+=table(r'Exact axial parameter and rank invariants of the finite-cell AILFT operator (cm$^{-1}$). Bulk = atom~2, interface = atom~13. The full complex Wybourne sets for these and the other archived calculations accompany the PDF as \texttt{wybourne\_Bkq\_all\_runs.csv}.','tab:S5','llrrrr',r'Ion & Site & $\Atwo$ & $S_2$ & $S_4$ & $S_6$',lines)

lines=[]
for ion,bs,gs in [('Nd',9,11),('Dy',21,23)]:
    for env,n,site in [('bulk',bs,'bulk 2'),('GB',gs,'GB 13')]:
        variants=[('finite cell',R[f'FINITE_{ion}_{env}'])]
        for radius,nr in [(18,n),(25,n+1)]:variants.append((f'sphere {radius} A',R[f'R{nr:02d}_{ion}_{env}_R{radius}_FeNeutral']))
        nn={'Nd':{'bulk':25,'GB':26},'Dy':{'bulk':27,'GB':28}}[ion][env]
        variants.append(('alternative charges',R[f'R{nn:02d}_{ion}_{env}_R25_FeBalt']))
        nk={'Nd':{'bulk':1,'GB':2},'Dy':{'bulk':3,'GB':4}}[ion][env]
        original_nev=next(r for k,r in N.items() if k.startswith(f'N{nk:02d}_'))
        variants.append(('finite cell + NEVPT2',original_nev))
        tk={'Nd':{'bulk':13,'GB':14},'Dy':{'bulk':15,'GB':16}}[ion][env]
        variants.append(('sphere 25 A + NEVPT2',next(r for k,r in N.items() if k.startswith(f'T{tk:02d}_'))))
        for lab,r in variants:lines.append(f"{ion} & {site} & {lab} & {num(r,'split'):.1f} & {num(r,'A20r2'):.1f} & {num(r,'E1'):.2f} \\\\")
alltables+=table(r'Selected sensitivity calculations at fixed centers (cm$^{-1}$). Sphere models use Fe $-3/14$, B $-3$ unless marked alternative (Fe $-5/14$, B $-1$). All use DKH2/SARC2-DKH-QZVP. NEVPT2 rows report the corrected block, not its accompanying CASSCF block. $E_{\mathrm{next}}$ is the first excited Kramers doublet for these two ions. The finite-cell versus sphere step changes multiple embedding choices; only the 18/25~\AA{} sphere comparison isolates radius. These two centers are not a matched parent-sublattice pair.','tab:S6','lllrrr',r'Ion & Center & Model & Orbital spread & $\Atwo$ & $E_{\mathrm{next}}$',lines)

lines=[]
for cls in ['A','B']:
    for j,ion in enumerate(IONS):
        if cls=='A' and ion=='Ce':continue
        br=f'T{j+1:02d}_{ion}_bulk4f_atom5_R25' if cls=='A' else f'R{[2,6,10,14,18,22][j]:02d}_{ion}_bulk_R25_FeNeutral'
        gr=f'R{[4,8,12,16,20,24][j]:02d}_{ion}_GB_R25_FeNeutral' if cls=='A' else f'T{j+7:02d}_{ion}_GBclassB_atom9_R25'
        b=R[br];g=R[gr]
        ang='--'
        if ion in ['Ce','Nd','Dy']:ang=f"{num(g,'g_angle'):.1f}"
        lines.append(f"{ion} & {cls} & {num(b,'split'):.1f}$\\to${num(g,'split'):.1f} & {num(b,'A20r2'):.1f}$\\to${num(g,'A20r2'):.1f} & {num(b,'E1'):.2f}$\\to${num(g,'E1'):.2f} & {ang} \\\\")
alltables+=table(r'Matched pairs in the 25~\AA{} sphere model (energies in cm$^{-1}$). A: bulk $4f$ atom~5 $\to$ interface class A atom~13. B: bulk $4g$ atom~2 $\to$ class B atom~9. $E_{\mathrm{next}}$ is defined as in Table~S4. The final column is the interface axis tilt from $c$ in degrees: direct ground-doublet magnetic matrix elements for Ce/Nd/Dy. Pr, Tb and Gd directions are omitted (Section~S4). Ce class A is unavailable.','tab:S7','llrrrr',r'Ion & Pair & Orbital spread (b$\to$g) & $\Atwo$ (b$\to$g) & $E_{\mathrm{next}}$ (b$\to$g) & Tilt',lines)

SITE8=[('bulk $4f$ (5)','bulk4f5'),('GB A (13)','GBA13'),('bulk $4g$ (2)','bulk4g2'),('GB B (9)','GBB9')]
DKH={('Nd','bulk4f5'):'T03_Nd_bulk4f_atom5_R25',('Nd','GBA13'):'R12_Nd_GB_R25_FeNeutral',('Nd','bulk4g2'):'R10_Nd_bulk_R25_FeNeutral',('Nd','GBB9'):'T09_Nd_GBclassB_atom9_R25',('Dy','bulk4f5'):'T06_Dy_bulk4f_atom5_R25',('Dy','GBA13'):'R24_Dy_GB_R25_FeNeutral',('Dy','bulk4g2'):'R22_Dy_bulk_R25_FeNeutral',('Dy','GBB9'):'T12_Dy_GBclassB_atom9_R25'}
caplines=[]
for ion,base in [('Nd',1),('Dy',9)]:
    for i,(lab,key) in enumerate(SITE8):
        pc=next(r for k,r in R.items() if k.startswith(f'U{base+2*i:02d}_'))
        ce=next(r for k,r in R.items() if k.startswith(f'U{base+2*i+1:02d}_'))
        vals=[num(pc,'split'),num(ce,'split'),num(pc,'A20r2'),num(ce,'A20r2'),num(pc,'S6'),num(ce,'S6'),num(R[DKH[(ion,key)]],'A20r2')]
        caplines.append((ion if i==0 else '')+' & '+lab+' & '+' & '.join(f'{v:.0f}' for v in vals)+r' \\')
alltables+=table(r'Capped effective-core-potential (cECP) test (cm$^{-1}$). All Fe and rare-earth neighbors within 4.0~\AA{} of the ion (18--19 centers) are replaced by capped ECPs with the same formal charges (Hay--Wadt for Fe, def2 small-core for the rare earth); pc = the same calculation with bare point charges. Level: def2-TZVP with the def2 small-core ECP on the central ion, spin-free SA-CASSCF + AILFT, $R = 25$~\AA. The last column is the sphere-model DKH2/SARC2-DKH-QZVP value with bare point charges.','tab:S8','llrrrrrrr',r'Ion & Site & split pc & split cECP & \Atwo{} pc & \Atwo{} cECP & $S_6$ pc & $S_6$ cECP & \Atwo{} DKH',caplines)

lines=[];swrows=[]
for r in csv.DictReader(open(AN/'tables/finite_cell_susceptibility_tensors.csv')):
    if r['environment']!='bulk' or r['ion']=='Gd':continue
    angle=float(r['bulk_GB_axis_angle_deg']);a=math.radians(angle)
    value=(abs(math.cos(a))**(2/3)+abs(math.sin(a))**(2/3))**(-1.5)
    lines.append(f"{r['ion']} & {angle:.2f} & {value:.3f} \\\\")
    swrows.append(dict(ion=r['ion'],angle_deg=angle,conditional_Hsw_over_HA=value))
alltables+=table(r'Conditional Stoner--Wohlfarth loss-of-metastability field for the finite-cell susceptibility-axis angles, assuming the applied field is along the bulk response axis. The listed value is dimensionless; it assumes an isolated uniaxial region with the same $K$ and $M_s$ at each angle. It is not a prediction of nucleation or macroscopic coercivity.','tab:S9','lrr',r'Ion & $\theta$ (degrees) & $H_{\mathrm{sw}}/H_A$',lines)
with open(AN/'tables/conditional_switching_estimate.csv','w',newline='') as h:
    w=csv.DictWriter(h,fieldnames=list(swrows[0]));w.writeheader();w.writerows(swrows)

lines=[]
for row in csv.DictReader(open(AN/'tables/orbital_extent.csv')):
    lines.append(f"{row['ion']} & {row['environment']} & {float(row['grid_norm']):.5f} & {float(row['orbital_r_rms_A']):.4f} \\\\")
alltables+=table(r'Full-cube integrals for the individual AILFT orbital in main Fig.~5. $N$ is the numerical one-electron norm before normalization; $r_{\mathrm{rms}}$ uses the full three-dimensional density and ion-centered distances. Cubes have $80^3$ points with spacing about 0.094~\AA. These are eigenorbital metrics of the finite-cell models, not total correlated-shell radii.','tab:S10','llrr',r'Ion & Environment & Grid norm $N$ & $r_{\mathrm{rms}}$ (\AA)',lines)

# Retain the manuscript's existing S8 typography while regenerating its numerical rows.
alltables=alltables.replace('\\centering\n\\footnotesize\n\\setlength{\\tabcolsep}{4pt}\n\\caption{Capped effective-core-potential', '\\centering\n\\caption{Capped effective-core-potential')
alltables=alltables.replace('\\label{tab:S8}\n', '\\label{tab:S8}\n\\small\n')

(OUT/'4f-supplementary.tex').write_text(tex,encoding='utf-8')
(OUT/'si_tables.tex').write_text(alltables,encoding='utf-8')
print('Wrote supplement with 2 figures and 10 tables. Conditional switching range:',min(r['conditional_Hsw_over_HA'] for r in swrows),max(r['conditional_Hsw_over_HA'] for r in swrows))
