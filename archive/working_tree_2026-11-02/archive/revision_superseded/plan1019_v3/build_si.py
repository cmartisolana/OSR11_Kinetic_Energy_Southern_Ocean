"""Build OSR11_Marti-Solana_etal_Supporting_Information_v3.0.docx from the v2.5 SI Word file (styles, footer, section
settings) and SI_v3.0.md, with figures new_S1..new_S10.png. References: the entries of the v2.5 main and SI lists (and the
two added in v3) that are cited in the SI text.
Usage: python build_si.py <unpacked SI v2.5> <unpacked main v2.5> <SI_v3.0.md> <figS dir> <out unpacked dir>"""
import sys, re, shutil, html
from pathlib import Path
from PIL import Image

SRC, MAIN, MD, FIG, DST = (Path(a) for a in sys.argv[1:6])
if DST.exists(): shutil.rmtree(DST)
shutil.copytree(SRC, DST)
X = (SRC / "word" / "document.xml").read_text(encoding="utf8")
head, rest = X.split("<w:body>", 1); body, tail = rest.rsplit("</w:body>", 1)
P = re.findall(r'<w:p[ >].*?</w:p>|<w:tbl>.*?</w:tbl>|<w:sectPr.*?</w:sectPr>', body, flags=re.S)
XM = (MAIN / "word" / "document.xml").read_text(encoding="utf8")
PM = re.findall(r'<w:p[ >].*?</w:p>|<w:tbl>.*?</w:tbl>|<w:sectPr.*?</w:sectPr>', XM.split("<w:body>", 1)[1], flags=re.S)


def accept(xml):
    xml = re.sub(r'<w:del [^>]*>.*?</w:del>', '', xml, flags=re.S)
    return re.sub(r'<w:ins [^>]*>(.*?)</w:ins>', r'\1', xml, flags=re.S)


def txt(p): return "".join(re.findall(r'<w:t[^>]*>([^<]*)</w:t>', p))
def esc(s): return html.escape(s, quote=False)


NOTE = re.compile(r'(\[(?:PENDING|AUTHOR CHECK|AUTHOR DECISION|CHECK|verify|TO RUN|NUMBER NEEDED)[^\]]*\])')
SUBS = {"T_hogg": ("T", "hogg"), "T_perm": ("T", "perm"), "T_all": ("T", "all"), "h_f": ("h", "f"), "n_eff": ("n", "eff"), "φhp": ("φ", "hp")}
SUB = re.compile("(" + "|".join(map(re.escape, SUBS)) + ")")


def run(t, rpr=""):
    return f'<w:r>{"<w:rPr>" + rpr + "</w:rPr>" if rpr else ""}<w:t xml:space="preserve">{esc(t)}</w:t></w:r>'


def runs(text, rpr=""):
    out = []
    for k, seg in enumerate(re.split(r'\*\*(.+?)\*\*', text)):
        r0 = rpr + ("<w:b/>" if k % 2 else "")
        for part in NOTE.split(seg):
            if not part: continue
            if NOTE.fullmatch(part): out.append(run(part, r0 + '<w:highlight w:val="yellow"/>')); continue
            for q in SUB.split(part):
                if not q: continue
                if q in SUBS:
                    b, s_ = SUBS[q]; out += [run(b, r0), run(s_, r0 + '<w:vertAlign w:val="subscript"/>')]
                else: out.append(run(q, r0))
    return "".join(out)


def para(text, ppr="", rpr=""): return f'<w:p>{"<w:pPr>" + ppr + "</w:pPr>" if ppr else ""}{runs(text, rpr)}</w:p>'
def heading(text, lvl=1): return f'<w:p><w:pPr><w:pStyle w:val="Heading{lvl}"/></w:pPr>{runs(text)}</w:p>'


rels = (DST / "word" / "_rels" / "document.xml.rels").read_text(encoding="utf8")
newrels = []; docpr = [200]


def figure_par(n):
    src = FIG / f"new_S{n}.png"; im = Image.open(src)
    w = 2400; h = int(im.height * w / im.width); name = f"image_v3_S{n}.png"
    im.convert("RGB").resize((w, h), Image.LANCZOS).save(DST / "word" / "media" / name, optimize=True)
    rid = f"rIdV3S{n}"; newrels.append(f'<Relationship Id="{rid}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="media/{name}"/>')
    cx = 6035040; cy = int(cx * h / w)
    if cy > 7600000: cy = 7600000; cx = int(cy * w / h)
    docpr[0] += 1; i = docpr[0]
    return (f'<w:p><w:pPr><w:keepNext/><w:jc w:val="center"/></w:pPr><w:r><w:drawing><wp:inline><wp:extent cx="{cx}" cy="{cy}"/>'
            f'<wp:docPr id="{i}" name="FigureS{n}"/><wp:cNvGraphicFramePr><a:graphicFrameLocks noChangeAspect="1"/></wp:cNvGraphicFramePr>'
            f'<a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:pic><pic:nvPicPr><pic:cNvPr id="{i}" name="FigureS{n}.png"/><pic:cNvPicPr/></pic:nvPicPr>'
            f'<pic:blipFill><a:blip r:embed="{rid}"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{cx}" cy="{cy}"/></a:xfrm>'
            f'<a:prstGeom prst="rect"/></pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>')


def placeholder_fig(text):
    bd = '<w:pBdr>' + "".join(f'<w:{s} w:val="single" w:sz="4" w:space="8" w:color="999999"/>' for s in ("top", "left", "bottom", "right")) + '</w:pBdr>'
    return f'<w:p><w:pPr><w:keepNext/>{bd}<w:spacing w:before="240" w:after="240"/><w:jc w:val="center"/></w:pPr>{runs(text)}</w:p>'


SZ = '<w:sz w:val="18"/><w:szCs w:val="18"/>'
BORD = ('<w:tblBorders>' + "".join(f'<w:{s} w:val="single" w:sz="4" w:space="0" w:color="000000"/>'
                                   for s in ("top", "left", "bottom", "right", "insideH", "insideV")) + '</w:tblBorders>')


def md_table(lines):
    rows = [[c.strip() for c in l.strip().strip("|").split("|")] for l in lines if l.strip().startswith("|") and not re.match(r'^\|\s*:?-', l.strip())]
    nc = max(len(r) for r in rows); total = 9900
    first = 3000 if nc > 2 else 4950; widths = [first] + [int((total - first) / (nc - 1))] * (nc - 1)
    grid = "".join(f'<w:gridCol w:w="{w}"/>' for w in widths)
    tbl = f'<w:tbl><w:tblPr><w:tblW w:w="{sum(widths)}" w:type="dxa"/>{BORD}<w:tblLayout w:type="fixed"/></w:tblPr><w:tblGrid>{grid}</w:tblGrid>'
    for k, r in enumerate(rows):
        r = r + [""] * (nc - len(r))
        tbl += '<w:tr>' + ("<w:trPr><w:tblHeader/></w:trPr>" if k == 0 else "")
        for c, w in zip(r, widths):
            tbl += (f'<w:tc><w:tcPr><w:tcW w:w="{w}" w:type="dxa"/><w:vAlign w:val="center"/></w:tcPr><w:p><w:pPr><w:spacing w:before="20" w:after="20"/>'
                    f'<w:jc w:val="{"center" if k == 0 else "left"}"/></w:pPr>{runs(c, ("<w:b/>" if k == 0 else "") + SZ)}</w:p></w:tc>')
        tbl += '</w:tr>'
    return tbl + '</w:tbl>'


md = re.sub(r'<!--.*?-->', '', MD.read_text(encoding="utf8"), flags=re.S)
lines = md.splitlines()
title = lines[0].lstrip("# ").strip()
blocks = [b.strip("\n") for b in re.split(r'\n\s*\n', "\n".join(lines[1:])) if b.strip()]
out = [f'<w:p><w:pPr/>{run(title, "<w:b/><w:sz w:val=\"28\"/>")}</w:p>']
i = 0
while i < len(blocks):
    b = blocks[i]
    if b.strip() == "---": pass
    elif b.startswith("## "): out.append(heading(b[3:], 1))
    elif b.startswith("### "): out.append(heading(b[4:], 2))
    elif re.match(r'\*\*Figure S(\d+)\.\*\*', b):
        n = int(re.match(r'\*\*Figure S(\d+)\.\*\*', b).group(1))
        if (FIG / f"new_S{n}.png").exists(): out.append(figure_par(n))
        else: out.append(placeholder_fig(f"[PENDING: Figure S{n} — Antarctic-slope index]"))
        out.append(para(b, '<w:spacing w:after="200"/>', SZ))
    elif b.lstrip().startswith("|"):
        out.append(md_table(b.splitlines()))
        out.append('<w:p><w:pPr><w:spacing w:after="120"/></w:pPr></w:p>')
    elif all(l.lstrip().startswith(("- ", "* ")) for l in b.splitlines()):
        for l in b.splitlines():
            out.append(para("– " + l.lstrip()[2:], '<w:ind w:left="360" w:hanging="240"/>'))
    else:
        bl = b.splitlines()
        if bl[0].startswith("**Table S") and any(l.lstrip().startswith("|") for l in bl[1:]):
            k = next(j for j, l in enumerate(bl) if l.lstrip().startswith("|"))
            out.append(para(" ".join(bl[:k]), '<w:keepNext/>', SZ)); out.append(md_table(bl[k:]))
            out.append('<w:p><w:pPr><w:spacing w:after="120"/></w:pPr></w:p>')
        elif bl[0].startswith("**Table S"):
            out.append(para(" ".join(bl), '<w:keepNext/>', SZ))
        else:
            out.append(para(" ".join(l.strip() for l in bl)))
    i += 1

# references: v2.5 main list + v2.5 SI list + added ones, those cited in the SI text
corpus = re.sub(r'\s+', ' ', md)
cand = {}
for p in [accept(p) for p in PM if re.match(r'^[A-ZÁÉÍÓÚ][^,]{1,40}, [A-Z]\.', txt(p)) and re.search(r'\d{4}\.\s*$', txt(p).strip())] + \
         [p for p in P if re.match(r'^[A-Z][^,]{1,40}, [A-Z]\.|^International Altimetry Team', txt(p))]:
    cand.setdefault(txt(p)[:80], p)
NEW = ['Hughes, C. W., Williams, J., Coward, A. C., and de Cuevas, B. A.: Antarctic circumpolar transport and the southern mode: a model investigation of interannual to decadal timescales, Ocean Sci., 10, 215–225, https://doi.org/10.5194/os-10-215-2014, 2014. [verify]']
for v in NEW: cand.setdefault(v[:80], para(v))
kept = []
for k, p in cand.items():
    t = txt(p) or k
    first = re.match(r'([^,:]+)[,:]', t).group(1).strip(); yr = re.findall(r'(\d{4})[a-z]?\.?\s*(?:\[verify\])?\s*$', t.strip())
    yr = yr[-1] if yr else ""
    sur = first.split()[-1] if not first.startswith("International") else "International Altimetry Team"
    if first.startswith("Thompson"): sur = "Thompson"
    if re.search(re.escape(sur) + r'.{0,60}?' + yr, corpus): kept.append((t, p))
kept.sort(key=lambda e: e[0].lower())
out.append(heading("References", 1)); out += [p for _, p in kept]
print("SI references:", len(kept)); print([t[:40] for t, _ in kept])
(DST / "word" / "document.xml").write_text(head + "<w:body>" + "".join(out) + P[-1] + "</w:body>" + tail, encoding="utf8")
rels = rels.replace("</Relationships>", "".join(newrels) + "</Relationships>")
(DST / "word" / "_rels" / "document.xml.rels").write_text(rels, encoding="utf8")
# hyperlinks inside copied main-text references need their relationships
mrels = (MAIN / "word" / "_rels" / "document.xml.rels").read_text(encoding="utf8")
used = set(re.findall(r'r:id="(rId\d+)"', "".join(p for _, p in kept)))
have = set(re.findall(r'Id="(rId\d+)"', rels))
add = []
for rid in used:
    m = re.search(rf'<Relationship Id="{rid}"[^>]*/>', mrels)
    if m and ('hyperlink' in m.group(0)):
        new_id = rid.replace("rId", "rIdM")
        add.append(m.group(0).replace(f'Id="{rid}"', f'Id="{new_id}"'))
doc = (DST / "word" / "document.xml").read_text(encoding="utf8")
# only main-list references carry main rIds; SI-list ones already match the SI rels
for _, p in kept:
    if p in [accept(q) for q in PM]:
        doc = doc.replace(p, re.sub(r'r:id="rId(\d+)"', r'r:id="rIdM\1"', p))
(DST / "word" / "document.xml").write_text(doc, encoding="utf8")
rels = rels.replace("</Relationships>", "".join(add) + "</Relationships>")
(DST / "word" / "_rels" / "document.xml.rels").write_text(rels, encoding="utf8")
print("blocks", len(out))
