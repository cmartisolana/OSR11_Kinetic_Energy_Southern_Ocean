"""Assemble the clean v2.2 docx: title block, abstract and table from the accepted v2.1 (clean/), body from
manuscript_v22.md via pandoc (reference styles = clean docx), equations and figure drawings re-inserted, references
from the accepted v2.1 with the edits listed in REF_EDITS."""
import copy, re, shutil, subprocess, sys
from pathlib import Path
from lxml import etree
from PIL import Image
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"; R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
q = lambda t: f"{{{W}}}{t}"
ROOT = Path("/home/claude/v3"); CLEAN = ROOT / "clean"; OUT = ROOT / "v22"; MD = ROOT / "reviews" / "manuscript_v22.md"
FIGS = Path("/mnt/user-data/outputs/revision/outputs/figures_v22")
shutil.rmtree(OUT, ignore_errors=True); shutil.copytree(CLEAN, OUT)

# ---- 1. body via pandoc
md = MD.read_text(encoding="utf8")
abstract = re.search(r"\*\*Abstract\.\*\*\s*(.+?)\n\n", md, re.S).group(1).strip()
body_md = md[md.index("# 1 Introduction"):]
body_md = re.sub(r"<!--.*?-->", "", body_md, flags=re.S)
(ROOT / "v22_body.md").write_text(body_md, encoding="utf8")
subprocess.run(["pandoc", str(ROOT / "v22_body.md"), "-o", str(ROOT / "v22_body.docx"), f"--reference-doc={ROOT}/v21b_clean.docx"], check=True)
subprocess.run(["rm", "-rf", str(ROOT / "v22_body_x")]); subprocess.run(["mkdir", "-p", str(ROOT / "v22_body_x")])
subprocess.run(["unzip", "-q", "-o", str(ROOT / "v22_body.docx"), "-d", str(ROOT / "v22_body_x")], check=True)
bt = etree.parse(str(ROOT / "v22_body_x/word/document.xml")); bbody = bt.getroot().find(q("body"))
body_pars = [p for p in bbody if p.tag == q("p")]

# ---- 2. clean document pieces
ct = etree.parse(str(OUT / "word/document.xml")); cbody = ct.getroot().find(q("body"))
cp = [p for p in cbody.iter(q("p"))]
txt = lambda p: "".join(x.text or "" for x in p.iter(q("t")))
i_abs = next(i for i, p in enumerate(cp) if txt(p).startswith("Abstract."))
i_tab = next(i for i, p in enumerate(cp) if txt(p).startswith("Product Table"))
i_tabcap = next(i for i, p in enumerate(cp) if txt(p).startswith("Table 1."))
i_intro = next(i for i, p in enumerate(cp) if txt(p).startswith("1 Introduction"))
i_refs = next(i for i, p in enumerate(cp) if txt(p).strip() == "References")
eqs = [p for p in cp if p.find(f".//{{{M}}}oMath") is not None and txt(p).strip() in ("(1)", "(2)")]
draw = next(p for p in cp if p.find(".//" + q("drawing")) is not None)
sectPr = cbody.find(q("sectPr"))
# top-level children of body: title block up to abstract, the table element, etc.
children = list(cbody)
def top_index_of(p):
    el = p
    while el.getparent() is not cbody: el = el.getparent()
    return children.index(el)
k_abs = top_index_of(cp[i_abs]); k_tab = top_index_of(cp[i_tab]); k_tabcap = top_index_of(cp[i_tabcap]); k_intro = top_index_of(cp[i_intro]); k_refs = top_index_of(cp[i_refs])
title_block = children[:k_abs]
table_block = children[k_tab:k_tabcap + 1]
ref_block = children[k_refs:]
ref_block = [e for e in ref_block if e.tag != q("sectPr")]

# ---- 3. helpers
def normal_par(text, template=cp[i_intro + 1], bold_prefix=None, size=None):
    p = etree.Element(q("p")); ppr = etree.SubElement(p, q("pPr"))
    tp = template.find(q("pPr"))
    if tp is not None:
        for el in tp:
            if el.tag != q("rPr"): ppr.append(copy.deepcopy(el))
    def run(t, b=False):
        r = etree.SubElement(p, q("r")); rp = etree.SubElement(r, q("rPr"))
        if b: etree.SubElement(rp, q("b"))
        if size:
            for tg in ("sz", "szCs"): e = etree.SubElement(rp, q(tg)); e.set(q("val"), str(size))
        t_ = etree.SubElement(r, q("t")); t_.text = t; t_.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    if bold_prefix: run(bold_prefix, True); run(text[len(bold_prefix):])
    else: run(text)
    return p
def fig_par(n):
    p = copy.deepcopy(draw)
    blip = p.find(".//{http://schemas.openxmlformats.org/drawingml/2006/main}blip"); blip.set(f"{{{R}}}embed", f"rIdV22Fig{n}")
    for dp in p.iter("{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}docPr"): dp.set("id", str(60 + n)); dp.set("name", f"Figure{n}.png")
    for cpr in p.iter("{http://schemas.openxmlformats.org/drawingml/2006/picture}cNvPr"): cpr.set("id", str(60 + n)); cpr.set("name", f"Figure{n}.png")
    im = Image.open(FIGS / f"Figure{n}.png"); w, h = im.size; maxw = 5715000; cx = maxw; cy = int(maxw * h / w)
    if cy > 8400000: cy = 8400000; cx = int(cy * w / h)
    for ext in p.iter("{http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing}extent"): ext.set("cx", str(cx)); ext.set("cy", str(cy))
    for ext in p.iter("{http://schemas.openxmlformats.org/drawingml/2006/main}ext"): ext.set("cx", str(cx)); ext.set("cy", str(cy))
    for sr in list(p.iter("{http://schemas.openxmlformats.org/drawingml/2006/main}srcRect")): sr.getparent().remove(sr)
    return p

# ---- 4. convert pandoc paragraphs
new_body = []
for p in body_pars:
    ppr = p.find(q("pPr")); st = ppr.find(q("pStyle")) if ppr is not None else None
    sid = st.get(q("val")) if st is not None else None
    text = txt(p).strip()
    if text.startswith("[EQUATION 1]"): new_body.append(copy.deepcopy(eqs[0])); continue
    if text.startswith("[EQUATION 2]"): new_body.append(copy.deepcopy(eqs[1])); continue
    m = re.match(r"Figure (\d)\. ", text)
    if m:
        n = int(m.group(1)); new_body.append(fig_par(n)); new_body.append(normal_par(text, bold_prefix=f"Figure {n}.", size=18)); continue
    if sid and sid.startswith("Heading"):
        # keep pandoc heading paragraph but normalise style id to the document's
        st.set(q("val"), sid.replace("Heading", "Heading")); new_body.append(copy.deepcopy(p)); continue
    if sid in ("BodyText", "FirstParagraph", "Compact") or sid is None:
        if ppr is not None and st is not None: ppr.remove(st)
        # drop pandoc run styles
        for rs in list(p.iter(q("rStyle"))): rs.getparent().remove(rs)
        new_body.append(copy.deepcopy(p)); continue
    new_body.append(copy.deepcopy(p))

# ---- 5. references edits
REF_DROP = ("Kendall, M. G.", "Mann, H. B.", "Yue, S. and Wang", "Meredith, M. P., Naveira Garabato", "Sallée, J.-B.", "Chelton, D. B., deSzoeke")
refs = [e for e in ref_block[1:] if e.tag == q("p")]
keep = []
for p in refs:
    t = txt(p)
    if any(t.startswith(d) for d in REF_DROP): continue
    for tnode in p.iter(q("t")):
        if tnode.text and "[check against the product page]" in tnode.text: tnode.text = tnode.text.replace(" [check against the product page]", "").replace("[check against the product page]", "")
    keep.append(p)
def sortkey(p):
    t = txt(p).strip(); m = re.match(r"(.+?):", t); authors = m.group(1) if m else t
    first = authors.split(",")[0].strip().lower(); year = re.search(r"(\d{4})[a-z]?\.?\s*$", t)
    return (first, authors.lower(), year.group(1) if year else "")
keep.sort(key=sortkey)

# ---- 6. assemble
for el in list(cbody): cbody.remove(el)
for el in title_block: cbody.append(el)
cbody.append(normal_par("Abstract. " + abstract, template=cp[i_abs], bold_prefix="Abstract."))
for el in table_block: cbody.append(el)
for el in new_body: cbody.append(el)
cbody.append(ref_block[0])
for p in keep: cbody.append(p)
cbody.append(sectPr)
ct.write(str(OUT / "word/document.xml"), xml_declaration=True, encoding="UTF-8", standalone=True)

# ---- 7. media + rels
rels = (OUT / "word/_rels/document.xml.rels").read_text(encoding="utf8")
for n in (1, 2, 3, 4):
    shutil.copy(FIGS / f"Figure{n}.png", OUT / f"word/media/v22_fig{n}.png")
    rels = rels.replace("</Relationships>", f'<Relationship Id="rIdV22Fig{n}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/v22_fig{n}.png"/></Relationships>')
(OUT / "word/_rels/document.xml.rels").write_text(rels, encoding="utf8")
subprocess.run(f"cd {OUT} && rm -f ../OSR11_Marti-Solana_Main_Manuscript_v2.2.docx && zip -q -X -r ../OSR11_Marti-Solana_Main_Manuscript_v2.2.docx '[Content_Types].xml' _rels docProps word customXml", shell=True, check=True)
print("built; refs kept", len(keep))
