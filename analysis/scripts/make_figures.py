from pathlib import Path
import csv, json, math, shutil, gzip, sys, argparse
import numpy as np
from PIL import Image
from contour_overlay import draw_contours, CONTOUR_LEVELS
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor, Color, white
from reportlab.lib.utils import ImageReader
import pypdfium2 as pdfium

ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description='Generate current figures from a reproduced analysis directory and included orbital cubes.')
parser.add_argument('analysis_dir',type=Path)
parser.add_argument('output_dir',type=Path)
args=parser.parse_args()
AN=args.analysis_dir.resolve()
OUT=args.output_dir.resolve()
OUT.mkdir(parents=True,exist_ok=True)
(AN/'figures').mkdir(parents=True,exist_ok=True)
(AN/'tables').mkdir(parents=True,exist_ok=True)
regular=Path('C:/Windows/Fonts/cambria.ttc')
bold=Path('C:/Windows/Fonts/cambriab.ttf')
if not regular.exists():
    from matplotlib.font_manager import FontProperties, findfont
    regular=Path(findfont(FontProperties(family='DejaVu Serif')))
    bold=Path(findfont(FontProperties(family='DejaVu Serif',weight='bold')))

pdfmetrics.registerFont(TTFont('A',str(regular)))
pdfmetrics.registerFont(TTFont('AB',str(bold)))
GRAY='#60666D';RED='#C7352A';BLUE='#2573A6';GREEN='#16806A';INK='#202124';LIGHT='#D9DDE1'
IONS=['Ce','Pr','Nd','Gd','Tb','Dy']
R={r['run']:r for r in csv.DictReader(open(AN/'tables/cf_all_runs.csv')) if r['block']=='CASSCF'}
R.update({r['run']:r for r in csv.DictReader(open(AN/'tables/finite_cell_descriptors.csv'))})
CHI=list(csv.DictReader(open(AN/'tables/finite_cell_chiT.csv')))
MAG=list(csv.DictReader(open(AN/'tables/finite_cell_susceptibility_tensors.csv')))
f=lambda run,key:float(R[run][key])

def txt(c,x,y,s,size=8,font='A',color=INK,align='left'):
    c.setFillColor(HexColor(color));c.setFont(font,size)
    getattr(c,{'left':'drawString','center':'drawCentredString','right':'drawRightString'}[align])(x,y,s)
def start(name,w,h):
    c=canvas.Canvas(str(AN/'figures'/f'{name}.pdf'),pagesize=(w,h));c.setTitle(name);return c
def finish(c,name):
    c.save();p=AN/'figures'/f'{name}.pdf';doc=pdfium.PdfDocument(p)
    doc[0].render(scale=3).to_pil().save(AN/'figures'/f'{name}.png')
    for ext in ['pdf','png']:shutil.copy2(AN/'figures'/f'{name}.{ext}',OUT/f'{name}.{ext}')
def axes(c,x,y,w,h,xlim,ylim,xticks,yticks,xlab='',ylab='',tag='',title=''):
    c.setStrokeColor(HexColor(GRAY));c.setLineWidth(.6);c.line(x,y,x+w,y);c.line(x,y,x,y+h)
    xx=lambda z:x+(z-xlim[0])/(xlim[1]-xlim[0])*w
    yy=lambda z:y+(z-ylim[0])/(ylim[1]-ylim[0])*h
    for v,lab in xticks:
        c.line(xx(v),y,xx(v),y-3);txt(c,xx(v),y-13,str(lab),7.5,align='center')
    for v,lab in yticks:
        c.line(x-3,yy(v),x,yy(v));txt(c,x-6,yy(v)-2.5,str(lab),7.5,align='right')
    if xlab:txt(c,x+w/2,y-29,xlab,8.3,align='center')
    if ylab:
        c.saveState();c.translate(x-35,y+h/2);c.rotate(90);txt(c,0,0,ylab,8.3,align='center');c.restoreState()
    if tag:txt(c,x+5,y+h-14,tag,9.5,'AB')
    if title:txt(c,x+w/2,y+h+9,title,9,'AB',align='center')
    return xx,yy
def line(c,pts,col=GRAY,dash=None,width=1.1):
    c.setStrokeColor(HexColor(col));c.setLineWidth(width);c.setDash(dash or [])
    p=c.beginPath();p.moveTo(*pts[0])
    for a in pts[1:]:p.lineTo(*a)
    c.drawPath(p);c.setDash([])
def marker(c,x,y,col,open=False,r=2):
    c.setStrokeColor(HexColor(col));c.setFillColor(white if open else HexColor(col));c.setLineWidth(.8);c.circle(x,y,r,stroke=1,fill=1)
def legend(c,x,y,s,col,open=False):
    line(c,[(x,y),(x+14,y)],col);marker(c,x+7,y,col,open);txt(c,x+20,y-2.5,s,7.5)

# Keep the original structure illustration without changing its geometry.
shutil.copy2(ROOT/'figures/Figure_1.png',OUT/'Figure_1.png')

# Original AILFT observable, re-plotted from the archived energies rather than fallback constants.
c=start('Figure_2',504,399)
bulk=[f('FINITE_'+i+'_bulk','split') for i in IONS];gb=[f('FINITE_'+i+'_GB','split') for i in IONS]
xx,yy=axes(c,49,238,177,130,(-.6,5.6),(0,3400),list(enumerate(IONS)),[(i,str(i)) for i in [0,1000,2000,3000]],ylab='4f orbital-energy spread (cm⁻¹)',tag='(a)')
for k,(a,b) in enumerate(zip(bulk,gb)):
    for dx,v,col in [(-.17,a,GRAY),(.17,b,RED)]:
        c.setFillColor(HexColor(col));c.rect(xx(k+dx)-3.7,yy(0),7.4,yy(v)-yy(0),fill=1,stroke=0)
    txt(c,xx(k),yy(max(a,b))+6,f'{100*(b/a-1):+.0f}%',7,align='center')
legend(c,80,382,'bulk',GRAY);legend(c,141,382,'GB',RED)
xx,yy=axes(c,307,238,177,130,(-.6,5.6),(.55,1.4),list(enumerate(IONS)),[(v,str(v)) for v in [.6,.8,1,1.2,1.4]],ylab='GB / bulk orbital-energy spread',tag='(b)')
line(c,[(xx(-.6),yy(1)),(xx(5.6),yy(1))],GRAY,[3,3],.6)
vals=[b/a for a,b in zip(bulk,gb)];line(c,[(xx(k),yy(v)) for k,v in enumerate(vals)],INK)
for k,v in enumerate(vals):marker(c,xx(k),yy(v),RED if v<1 else BLUE)
for k,ion in enumerate(['Ce','Dy']):
    x=49+258*k; mx=2200 if ion=='Ce' else 700
    ticks=[0,500,1000,1500,2000] if ion=='Ce' else [0,200,400,600]
    xx,yy=axes(c,x,49,177,128,(-.6,1.7),(0,mx),[(0,'bulk'),(1,'GB')],[(v,str(v)) for v in ticks],ylab='Orbital energy above lowest (cm⁻¹)',tag='(c)' if k==0 else '(d)',title=ion+'³⁺  one-electron levels')
    for j,env,col in [(0,'bulk',GRAY),(1,'GB',RED)]:
        d=json.loads((AN/'parsed'/f'FINITE_{ion}_{env}.json').read_text())
        E=d['ailft']['CASSCF']['E'][:3]
        for ev in E:line(c,[(xx(j)-19,yy(ev)),(xx(j)+19,yy(ev))],col,width=1.8)
        txt(c,xx(j)+24,yy(E[1])+3,f'δ = {E[1]:.0f}',6.8,color=col)
        txt(c,xx(j),yy(E[2])+5,f'Δ = {E[2]-E[1]:.0f}',7,color=col,align='center')
finish(c,'Figure_2')

# Magnetic directions from the actual zero-field susceptibility tensors of the same finite model.
c=start('Figure_3',504,215)
txt(c,17,195,'(a)',9.5,'AB');txt(c,117,195,'Local single-ion response',9,'AB',align='center')
c.setFillColor(HexColor('#F2F4F6'));c.setStrokeColor(HexColor(LIGHT));c.ellipse(23,36,210,88,stroke=1,fill=1)
line(c,[(116,61),(116,171)],GRAY,width=2)
line(c,[(116,171),(112,162)],GRAY,width=1.6);line(c,[(116,171),(120,162)],GRAY,width=1.6)
line(c,[(116,61),(200,66)],RED,width=2)
line(c,[(200,66),(192,70)],RED,width=1.6);line(c,[(200,66),(191,62)],RED,width=1.6)
txt(c,123,157,'bulk ≈ c',8,color=GRAY);txt(c,154,80,'GB',8,color=RED)
txt(c,127,108,'≈ 87°',10,'AB');txt(c,126,31,'basal plane',8,color=GRAY)
txt(c,117,11,'Finite-cell embedding',8,align='center')
xx,yy=axes(c,285,49,201,127,(-.6,5.6),(0,110),list(enumerate(IONS)),[(v,str(v)) for v in [0,30,60,90]],ylab='Bulk–GB axis angle at 2 K (deg)',tag='(b)')
for k,ion in enumerate(IONS):
    if ion=='Gd':txt(c,xx(k),yy(13),'—',10,align='center');continue
    v=float(next(r for r in MAG if r['ion']==ion)['bulk_GB_axis_angle_deg'])
    c.setFillColor(HexColor(RED if k<3 else BLUE));c.rect(xx(k)-8,yy(0),16,yy(v)-yy(0),fill=1,stroke=0)
    txt(c,xx(k),yy(v)+4,f'{v:.1f}°',7,align='center')
txt(c,383,11,'Gd: no orbital easy-axis assignment',7,align='center')
finish(c,'Figure_3')

# Six original susceptibility panels, now uniformly lettered at upper left.
c=start('Figure_4',504,323)
curie=[.8035714286,1.6,1.6363636364,7.875,11.8125,14.1666666667]
limits=[(.36,.86),(.1,1.75),(1.18,1.72),(7.38,7.94),(9.8,12.1),(11.9,14.45)]
ticks=[[.4,.6,.8],[.25,.75,1.25,1.75],[1.2,1.4,1.6],[7.4,7.6,7.8],[10,11,12],[12,13,14]]
for k,ion in enumerate(IONS):
    x=46+(k%3)*167;y=191-(k//3)*150
    xx,yy=axes(c,x,y,115,104,(0,40),limits[k],[(v,str(v)) for v in [0,10,20,30,40]],[(v,str(v)) for v in ticks[k]],xlab='Temperature (K)' if k>=3 else '',ylab='χM T (cm³ K mol⁻¹)' if k%3==0 else '',title=ion+'³⁺')
    txt(c,x+4,y+113,f'({"abcdef"[k]})',9.5,'AB')
    line(c,[(xx(0),yy(curie[k])),(xx(40),yy(curie[k]))],'#9AA0A6',[3,3],.7)
    txt(c,xx(39),yy(curie[k])+3,'free ion',6.5,color='#8A9098',align='right')
    for env,col,op in [('bulk',GRAY,True),('GB',RED,False)]:
        data=[r for r in CHI if r['ion']==ion and r['environment']==env]
        pts=[(xx(float(r['T_K'])),yy(float(r['chiT']))) for r in data];line(c,pts,col,None,.65)
        for px,py in pts:marker(c,px,py,col,op,1.6)
    if k==0:legend(c,x+50,y+64,'bulk',GRAY,True);legend(c,x+50,y+50,'GB',RED)
finish(c,'Figure_4')

# Figure 5: bulk/interface comparison on one absolute density scale, from the archived cubes.
c=start('Figure_5',550,405);cube_meta=[];maps={};bohr=.5291772109
for ion in ['Ce','Pr','Dy']:
    for env in ['Bulk','GB']:
        path=ROOT/'analysis/orbitals'/f'{ion}_{env}_mo27a.cube.gz'
        with gzip.open(path,'rt') as h:
            h.readline();h.readline();a=h.readline().split();nat=abs(int(a[0]));origin=np.array(a[1:4],float)
            grid=[h.readline().split() for _ in range(3)];n=np.array([int(q[0]) for q in grid]);step=np.array([q[1:4] for q in grid],float)
            atoms=[h.readline().split() for _ in range(nat)];center=np.array(atoms[0][2:5],float)
            if int(a[0])<0:h.readline()
            data=np.fromstring(h.read(),sep=' ').reshape(n)
        coords=[origin[j]+np.arange(n[j])*step[j,j]-center[j] for j in range(3)]
        rho=data**2;voxel=abs(np.linalg.det(step));norm=float(rho.sum()*voxel)
        r2grid=sum(coords[j].reshape(tuple(n[j] if i==j else 1 for i in range(3)))**2 for j in range(3))
        rms=float(np.sqrt((rho*r2grid).sum()*voxel/norm))*bohr
        xyz=[z*bohr for z in coords]
        keepx=np.where(abs(xyz[0])<=1.08)[0];keepy=np.where(abs(xyz[1])<=1.08)[0]
        # Integrate the full density along c; a centred slice can be a nodal plane.
        projection=(rho.sum(axis=2)*abs(step[2,2])/norm/bohr**2).T
        density=projection[np.ix_(keepy,keepx)]
        maps[(ion,env)]=(xyz[0][keepx],xyz[1][keepy],density,rms)
        cube_meta.append(dict(ion=ion,environment=env,source=path.relative_to(ROOT).as_posix(),plane='Full orbital density projected along c onto a,b',grid_norm=norm,orbital_r_rms_A=rms,voxel_bohr3=voxel))
vmax=math.ceil(max(v[2].max() for v in maps.values())*2)/2
# Sequential white-to-blue palette; a single square-root display mapping for every panel.
PALETTE=np.array([[255,255,255],[222,235,247],[198,219,239],[158,202,225],
                  [107,174,214],[66,146,198],[33,113,181],[8,81,156],[8,48,107]])
def density_rgb(density):
    coordinate=np.sqrt(np.clip(np.asarray(density)/vmax,0,1))
    return np.stack([np.interp(coordinate,np.linspace(0,1,len(PALETTE)),PALETTE[:,k]) for k in range(3)],axis=-1)
assert all(v[2].min()>=0 and v[2].max()<=vmax for v in maps.values())
txt(c,264,386,'Projected lowest-orbital density — shared square-root color scale',10,align='center')
for j,ion in enumerate(['Ce','Pr','Dy']):
    for row,env in enumerate(['Bulk','GB']):
        xs,ys,density,rms=maps[(ion,env)]
        # No separate panel normalization; only the shared density-to-color transfer changes.
        rgb=density_rgb(density).clip(0,255).astype('uint8')
        im=Image.fromarray(rgb[::-1]).resize((550,550),Image.Resampling.BICUBIC)
        im.save(AN/'figures'/f'{ion}_{env}_orbital_density.png')
        x=48+j*163;y=237-row*180;w=118;h=118
        c.drawImage(ImageReader(im),x,y,w,h)
        draw_contours(c,density,x,y,w,h)
        ticks=[(v,f'{v:g}') for v in [-1,-.5,0,.5,1]]
        xx,yy=axes(c,x,y,w,h,(xs[0],xs[-1]),(ys[0],ys[-1]),ticks,ticks,xlab='a direction (Å)' if row==1 else '',ylab=('Bulk' if row==0 else 'Boundary')+'   b (Å)' if j==0 else '',tag=f'({"abcdef"[row*3+j]})',title=f'{ion}³⁺  {env} MO27')
        # Labeled 0.5-Angstrom intervals plus minor ticks every 0.25 Angstrom.
        c.setStrokeColor(HexColor(GRAY));c.setLineWidth(.5)
        for v in [-.75,-.25,.25,.75]:
            c.line(xx(v),y,xx(v),y-1.8)
            c.line(x-1.8,yy(v),x,yy(v))
        txt(c,x+w/2,y-(29 if row==0 else 42),'orbital rᵣₘₛ = %.2f Å'%rms,8,align='center')
# Continuous colorbar; its nonlinear tick positions are labelled with physical density.
gradient=density_rgb(vmax*np.linspace(1,0,800)**2).astype('uint8')[:,None,:]
c.drawImage(ImageReader(Image.fromarray(np.repeat(gradient,12,axis=1))),505,65,8,278)
for value in [0,.25,1,2,4,vmax]:
    py=65+278*np.sqrt(value/vmax)
    c.setStrokeColor(HexColor(GRAY));c.setLineWidth(.5);c.line(513,py,516,py)
    txt(c,519,py-2.5,f'{value:g}',7)
c.saveState();c.translate(545,204);c.rotate(90);txt(c,0,0,'Projected density (Å⁻²)',8,align='center');c.restoreState()
(AN/'tables/orbital_slice_provenance.json').write_text(json.dumps(dict(common_density_vmin_Aminus2=0,common_density_vmax_Aminus2=vmax,color_mapping='square root: u=(P/Pmax)^(1/2)',gamma=0.5,independent_panel_normalization=False,contour_levels_Aminus2=list(CONTOUR_LEVELS),contour_display_interpolation='bicubic',orbitals=cube_meta),indent=2),encoding='utf-8')
with open(AN/'tables/orbital_extent.csv','w',newline='') as h:
    ww=csv.DictWriter(h,fieldnames=['ion','environment','grid_norm','orbital_r_rms_A','source']);ww.writeheader();ww.writerows([{k:r[k] for k in ww.fieldnames} for r in cube_meta])
finish(c,'Figure_5')

# S1: show original and extended embeddings explicitly; do not label a multi-change comparison a radius test.
c=start('Figure_S1',504,396)
for k,key,label,ylim,ticks in [(0,'A20r2','A₂⁰⟨r²⟩ at GB atom 13 (cm⁻¹)',(-1100,1300),[-1000,0,1000]),(1,'split','4f orbital-energy spread at GB (cm⁻¹)',(500,2100),[500,1000,1500,2000])]:
    x=49+258*k;xx,yy=axes(c,x,237,177,128,(-.25,2.25),ylim,[(0,'finite cell'),(1,'18 Å'),(2,'25 Å')],[(v,str(v)) for v in ticks],ylab=label,tag=f'({"ab"[k]})')
    if ylim[0]<0:line(c,[(xx(-.25),yy(0)),(xx(2.25),yy(0))],GRAY,[3,3],.6)
    for ion,ns,col,op in [('Nd',[None,11,12],RED,False),('Dy',[None,23,24],BLUE,True)]:
        vals=[f('FINITE_'+ion+'_GB',key)]+[f(f'R{n:02d}_{ion}_GB_R{18 if j==1 else 25}_FeNeutral',key) for j,n in enumerate(ns) if n]
        line(c,[(xx(j),yy(v)) for j,v in enumerate(vals)],col)
        for j,v in enumerate(vals):marker(c,xx(j),yy(v),col,op)
    if k==0:legend(c,x+10,379,'Nd',RED);legend(c,x+83,379,'Dy',BLUE,True)
for k,key,lab in [(0,'split','Matched GB / bulk orbital spread'),(1,'A20r2','Matched GB / bulk A₂⁰⟨r²⟩')]:
    x=49+258*k;xx,yy=axes(c,x,51,177,126,(-.6,5.6),(.8,1.7),list(enumerate(IONS)),[(v,str(v)) for v in [.8,1,1.2,1.4,1.6]],ylab=lab,tag=f'({"cd"[k]})')
    line(c,[(xx(-.6),yy(1)),(xx(5.6),yy(1))],GRAY,[3,3],.6)
    for cls,col,op in [('A',RED,False),('B',GREEN,True)]:
        pts=[]
        for j,ion in enumerate(IONS):
            if cls=='A' and ion=='Ce':continue
            br=f'T{j+1:02d}_{ion}_bulk4f_atom5_R25' if cls=='A' else f'R{[2,6,10,14,18,22][j]:02d}_{ion}_bulk_R25_FeNeutral'
            gr=f'R{[4,8,12,16,20,24][j]:02d}_{ion}_GB_R25_FeNeutral' if cls=='A' else f'T{j+7:02d}_{ion}_GBclassB_atom9_R25'
            pts.append((xx(j),yy(f(gr,key)/f(br,key))))
        line(c,pts,col)
        for p in pts:marker(c,*p,col,op)
    if k==0:legend(c,x,192,'A / bulk 4f',RED);legend(c,x+94,192,'B / bulk 4g',GREEN,True)
finish(c,'Figure_S1')

# S2: exact rank decomposition, with no claimed radial-density measurement.
c=start('Figure_S2',504,219)
for panel,cls,ttl in [(0,'original','Finite-cell site comparison'),(1,'A','Sphere: A / bulk 4f'),(2,'B','Sphere: B / bulk 4g')]:
    x=43+169*panel;xx,yy=axes(c,x,51,119,126,(-.6,5.6),(-105,60),list(enumerate(IONS)),[(v,str(v)) for v in [-100,-50,0,50]],ylab='Change in rank strength Sₖ (%)' if panel==0 else '',tag=f'({"abc"[panel]})',title=ttl)
    line(c,[(xx(-.6),yy(0)),(xx(5.6),yy(0))],GRAY,width=.5)
    for j,ion in enumerate(IONS):
        if cls=='A' and ion=='Ce':continue
        if cls=='original':br='FINITE_'+ion+'_bulk';gr='FINITE_'+ion+'_GB'
        elif cls=='A':br=f'T{j+1:02d}_{ion}_bulk4f_atom5_R25';gr=f'R{[4,8,12,16,20,24][j]:02d}_{ion}_GB_R25_FeNeutral'
        else:br=f'R{[2,6,10,14,18,22][j]:02d}_{ion}_bulk_R25_FeNeutral';gr=f'T{j+7:02d}_{ion}_GBclassB_atom9_R25'
        for n,(rank,col) in enumerate([(2,'#215A95'),(4,'#6299C5'),(6,'#B8D4E8')]):
            val=100*(f(gr,f'S{rank}')/f(br,f'S{rank}')-1);c.setFillColor(HexColor(col));c.rect(xx(j+(n-1)*.24)-1.9,min(yy(0),yy(val)),3.8,abs(yy(val)-yy(0)),fill=1,stroke=0)
for j,(rank,col) in enumerate([(2,'#215A95'),(4,'#6299C5'),(6,'#B8D4E8')]):
    c.setFillColor(HexColor(col));c.rect(154+j*75,201,9,7,fill=1,stroke=0);txt(c,168+j*75,201,f'k = {rank}',8)
finish(c,'Figure_S2')
print('Built Figures 2–5 and S1–S2 from archived data; retained Figure 1.')
