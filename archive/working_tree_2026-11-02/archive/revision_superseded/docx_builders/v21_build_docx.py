"""Build v2.1 docx with paragraph-level tracked changes from base (accepted v2) + changes.py + figures."""
import copy, datetime, re, shutil, sys
from lxml import etree
sys.path.insert(0, '/home/claude/v3/tools'); from changes_b import REPLACE, INSERT_AFTER, NEW_REFS, NEW_TABLE_ROW
import os, pandas as pd
_ts='/mnt/user-data/outputs/revision/outputs/a8/a8_trends_twosat.csv'
if os.path.exists(_ts):
    _d=pd.read_csv(_ts).set_index('series').loc['ACC_T']; _txt=f" and by {_d.pct_per_decade:.1f} % per decade ({_d.pct_lo:.1f}–{_d.pct_hi:.1f}; p = {_d.p:.3f}) in the two-satellite product"
else: _txt=''
REPLACE[78]=REPLACE[78].replace('__TWOSAT_T__',_txt)
from PIL import Image
W="http://schemas.openxmlformats.org/wordprocessingml/2006/main"; R="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
q=lambda t:f"{{{W}}}{t}"
AUTH="Claude (v2.1 draft)"; DATE=datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
nid=[9000]
def att(el):
    nid[0]+=1; el.set(q("id"),str(nid[0])); el.set(q("author"),AUTH); el.set(q("date"),DATE)
SRC='/home/claude/v3/base'; DST='/home/claude/v3/out'
shutil.rmtree(DST, ignore_errors=True); shutil.copytree(SRC, DST)
tree=etree.parse(f'{DST}/word/document.xml'); body=tree.getroot().find(q('body'))
ps=[p for p in body.iter(q('p'))]
def mark_deleted(p):
    ppr=p.find(q('pPr'))
    if ppr is None:
        ppr=etree.SubElement(p,q('pPr')); p.insert(0,ppr)
    rpr=ppr.find(q('rPr'))
    if rpr is None:
        rpr=etree.SubElement(ppr,q('rPr'))
    d=etree.Element(q('del')); att(d); rpr.insert(0,d)
    def wrap_runs(container):
        for r in list(container):
            if r.tag==q('pPr'): continue
            if r.tag==q('hyperlink'): wrap_runs(r); continue
            for t in r.iter(q('t')): t.tag=q('delText')
            wrap=etree.Element(q('del')); att(wrap); container.replace(r,wrap); wrap.append(r)
    wrap_runs(p)
def new_par(template, text, style=None):
    p=etree.Element(q('p'))
    tp=template.find(q('pPr'))
    ppr=copy.deepcopy(tp) if tp is not None else etree.Element(q('pPr'))
    for el in list(ppr):
        if el.tag in (q('rPr'),): ppr.remove(el)
    if style:
        st=ppr.find(q('pStyle'))
        if st is None: st=etree.Element(q('pStyle')); ppr.insert(0,st)
        st.set(q('val'),style)
    rpr=etree.SubElement(ppr,q('rPr')); i=etree.SubElement(rpr,q('ins')); att(i)
    p.append(ppr)
    ins=etree.SubElement(p,q('ins')); att(ins)
    # caption style: bold "Figure N." + 9pt as in v2 captions
    m=re.match(r'(Figure \d+\.\s)(.*)',text,flags=re.S)
    def run(txt,bold=False,small=False):
        r=etree.SubElement(ins,q('r')); rp=etree.SubElement(r,q('rPr'))
        if bold: etree.SubElement(rp,q('b'))
        if small:
            for tg in ('sz','szCs'):
                e=etree.SubElement(rp,q(tg)); e.set(q('val'),'18')
        t=etree.SubElement(r,q('t')); t.text=txt; t.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
    if m: run(m.group(1),True,True); run(m.group(2),False,True)
    else: run(text)
    return p
# figure drawing template & image replacement
draw_pars={}
for i,p in enumerate(ps):
    if p.find('.//'+q('drawing')) is not None: draw_pars[i]=p
figmap={60:'Figure1',67:'Figure2',74:'Figure3'}
rels=open(f'{DST}/word/_rels/document.xml.rels',encoding='utf8').read()
EMU=914400; maxw=5715000
def set_drawing(p, png):
    im=Image.open(png); w,h=im.size
    cx=maxw; cy=int(maxw*h/w)
    if cy>8400000: cy=8400000; cx=int(cy*w/h)
    for ext in p.iter('{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}extent'):
        ext.set('cx',str(cx)); ext.set('cy',str(cy))
    for ext in p.iter('{http://schemas.openxmlformats.org/drawingml/2006/main}ext'):
        ext.set('cx',str(cx)); ext.set('cy',str(cy))
    for sr in list(p.iter('{http://schemas.openxmlformats.org/drawingml/2006/main}srcRect')):
        sr.getparent().remove(sr)
for i,name in figmap.items():
    p=draw_pars[i]; blip=p.find('.//{http://schemas.openxmlformats.org/drawingml/2006/main}blip')
    rid=blip.get(f'{{{R}}}embed'); target=re.search(rf'Id="{rid}"[^>]*Target="([^"]+)"',rels) or re.search(rf'Target="([^"]+)"[^>]*Id="{rid}"',rels)
    tgt=target.group(1)
    src=f'/mnt/user-data/outputs/revision/outputs/figures_v21/{name}.png'
    out=f'{DST}/word/{tgt}'
    Image.open(src).convert('RGB').save(out, quality=92) if tgt.lower().endswith(('.jpg','.jpeg')) else shutil.copy(src,out)
    set_drawing(p, src)
# Figures 4 and 5: new media + rels
FIGSRC='/mnt/user-data/outputs/revision/outputs/figures_v21'
for n in (4,5):
    shutil.copy(f'{FIGSRC}/Figure{n}.png', f'{DST}/word/media/image_fig{n}.png')
    rels=rels.replace('</Relationships>',f'<Relationship Id="rIdFig{n}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/image_fig{n}.png"/></Relationships>')
open(f'{DST}/word/_rels/document.xml.rels','w',encoding='utf8').write(rels)
ct=open(f'{DST}/[Content_Types].xml',encoding='utf8').read()
if 'Extension="png"' not in ct: ct=ct.replace('</Types>','<Default Extension="png" ContentType="image/png"/></Types>'); open(f'{DST}/[Content_Types].xml','w',encoding='utf8').write(ct)
def fig_par(template_draw, n):
    p=copy.deepcopy(template_draw)
    blip=p.find('.//{http://schemas.openxmlformats.org/drawingml/2006/main}blip'); blip.set(f'{{{R}}}embed',f'rIdFig{n}')
    for dp in p.iter('{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}docPr'): dp.set('id',str(40+n)); dp.set('name',f'Figure{n}.png')
    for cp in p.iter('{http://schemas.openxmlformats.org/drawingml/2006/picture}cNvPr'): cp.set('id',str(40+n)); cp.set('name',f'Figure{n}.png')
    set_drawing(p,f'{FIGSRC}/Figure{n}.png')
    # mark inserted
    ppr=p.find(q('pPr')); rpr=ppr.find(q('rPr'))
    if rpr is None: rpr=etree.SubElement(ppr,q('rPr'))
    i=etree.Element(q('ins')); att(i); rpr.insert(0,i)
    for r in list(p):
        if r.tag==q('pPr'): continue
        wrap=etree.Element(q('ins')); att(wrap); p.replace(r,wrap); wrap.append(r)
    return p
# Table 1: Product 5 row (tracked insertion), cloned from the row of Product 4
row4=ps[28].getparent().getparent(); assert row4.tag==q('tr')
row5=copy.deepcopy(row4)
cells=[tc for tc in row5 if tc.tag==q('tc')]
texts=[NEW_TABLE_ROW[0],NEW_TABLE_ROW[1],NEW_TABLE_ROW[2]]+[NEW_TABLE_ROW[3]]
for tc,txt in zip(cells,texts):
    pars=[p for p in tc if p.tag==q('p')]
    lines=txt if isinstance(txt,list) else [txt]
    for extra in pars[len(lines):]: tc.remove(extra)
    for p,line in zip(pars,lines):
        for h in list(p.iter(q('hyperlink'))): h.getparent().remove(h)   # plain text instead of the Product-4 hyperlink
        runs=[r for r in p if r.tag==q('r')]
        for r in runs[1:]: p.remove(r)
        if runs:
            r=runs[0]
            for t in list(r.iter(q('t'))): r.remove(t)
        else:
            r=etree.SubElement(p,q('r')); rp=etree.SubElement(r,q('rPr'))
            for tg in ('sz','szCs'): e=etree.SubElement(rp,q(tg)); e.set(q('val'),'18')
        t=etree.SubElement(r,q('t')); t.text=line; t.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
        wrap=etree.Element(q('ins')); att(wrap); p.replace(r,wrap); wrap.append(r)
        ppr=p.find(q('pPr')); rpr=ppr.find(q('rPr')) if ppr is not None else None
        if ppr is None: ppr=etree.Element(q('pPr')); p.insert(0,ppr)
        if rpr is None: rpr=etree.SubElement(ppr,q('rPr'))
        i=etree.Element(q('ins')); att(i); rpr.insert(0,i)
trpr=row5.find(q('trPr'))
if trpr is None: trpr=etree.Element(q('trPr')); row5.insert(0,trpr)
i=etree.SubElement(trpr,q('ins')); att(i)
row4.addnext(row5)
# apply replacements
for idx,new in REPLACE.items():
    p=ps[idx]; mark_deleted(p)
    if new is not None:
        p.addnext(new_par(p,new))
# insertions
for idx,items in INSERT_AFTER.items():
    anchor=ps[idx]
    # if anchor was replaced, insert after its replacement
    nxt=anchor.getnext()
    if idx in REPLACE and REPLACE[idx] is not None: anchor=nxt
    for style,text in items:
        if text in ('__FIGURE4__','__FIGURE5__'): np_=fig_par(draw_pars[74], int(text[8]))
        else: np_=new_par(ps[idx] if style=='Normal' else ps[43], text, style)
        anchor.addnext(np_); anchor=np_
# new references alphabetically
refs=[(i,''.join(x.text or '' for x in p.iter(q('t')))) for i,p in enumerate(ps) if 102<=i<=139 and i not in REPLACE]
for ref in NEW_REFS:
    key=ref.split(',')[0].lower(); prev=ps[101]
    for i,txt in refs:
        if txt.split(',')[0].lower()<key: prev=ps[i]
    # if prev was replaced, use its replacement
    anchor=prev if prev.getparent() is not None else None
    if prev.find(q('pPr')) is not None and prev.find(q('pPr')).find(q('rPr')) is not None and prev.find(q('pPr')).find(q('rPr')).find(q('del')) is not None:
        prev=prev.getnext()
    prev.addnext(new_par(ps[102],ref))
tree.write(f'{DST}/word/document.xml',xml_declaration=True,encoding='UTF-8',standalone=True)
print('built')
