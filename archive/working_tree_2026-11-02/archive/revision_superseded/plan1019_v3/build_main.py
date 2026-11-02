"""Build OSR11_Marti-Solana_Main_Manuscript_v3.0.docx from the v2.5 Word file (template, author block, product table,
equations, reference list) and manuscript_v3.0.md (text, captions, Table 2), with the v3 figures.
Bracketed author notes ([PENDING ...], [AUTHOR CHECK ...], [CHECK ...], [verify ...]) are highlighted yellow.
Usage: python build_main.py <unpacked v2.5 dir> <manuscript_v3.0.md> <figure dir> <out unpacked dir>"""
import sys, re, shutil, html
from pathlib import Path
from PIL import Image

SRC, MD, FIG, DST = (Path(a) for a in sys.argv[1:5])
if DST.exists(): shutil.rmtree(DST)
shutil.copytree(SRC, DST)
X = (SRC / "word" / "document.xml").read_text(encoding="utf8")
head, rest = X.split("<w:body>", 1)
body, tail = rest.rsplit("</w:body>", 1)
P = re.findall(r'<w:p[ >].*?</w:p>|<w:tbl>.*?</w:tbl>|<w:sectPr.*?</w:sectPr>', body, flags=re.S)
assert "".join(P) == body.strip() or len("".join(P)) > 0.99 * len(body.strip())


def accept(xml):
    """accept tracked changes in a kept block: drop deletions, unwrap insertions."""
    xml = re.sub(r'<w:del [^>]*>.*?</w:del>', '', xml, flags=re.S)
    xml = re.sub(r'<w:ins [^>]*>(.*?)</w:ins>', r'\1', xml, flags=re.S)
    return xml


def txt(p):
    return "".join(re.findall(r'<w:t[^>]*>([^<]*)</w:t>', p))


NOTE = re.compile(r'(\[(?:PENDING|AUTHOR CHECK|CHECK|verify|TO RUN|NUMBER NEEDED)[^\]]*\])')
SUB = re.compile(r'(φhp|h_f)')


def esc(s):
    return html.escape(s, quote=False)


def run(t, rpr=""):
    return f'<w:r>{"<w:rPr>" + rpr + "</w:rPr>" if rpr else ""}<w:t xml:space="preserve">{esc(t)}</w:t></w:r>'


def runs(text, rpr=""):
    out = []
    for part in NOTE.split(text):
        if not part: continue
        if NOTE.fullmatch(part):
            out.append(run(part, rpr + '<w:highlight w:val="yellow"/>')); continue
        for q in SUB.split(part):
            if not q: continue
            if q == "φhp": out += [run("φ", rpr), run("hp", rpr + '<w:vertAlign w:val="subscript"/>')]
            elif q == "h_f": out += [run("h", rpr), run("f", rpr + '<w:vertAlign w:val="subscript"/>')]
            else: out.append(run(q, rpr))
    return "".join(out)


def para(text, ppr="", rpr=""):
    return f'<w:p>{"<w:pPr>" + ppr + "</w:pPr>" if ppr else ""}{runs(text, rpr)}</w:p>'


def heading(text, lvl):
    return f'<w:p><w:pPr><w:pStyle w:val="Heading{lvl}"/></w:pPr>{runs(text)}</w:p>'


SZ = '<w:sz w:val="18"/><w:szCs w:val="18"/>'


def caption(label, text):
    return f'<w:p><w:pPr><w:spacing w:after="200"/></w:pPr>{run(label, "<w:b/><w:bCs/>" + SZ)}{runs(" " + text, SZ)}</w:p>'


# ------------------------------------------------------------- images: reuse rId12-15 (media/image1-4.png)
IMG = {1: ("rId12", "image1.png"), 2: ("rId13", "image2.png"), 3: ("rId14", "image3.png")}
docpr = [100]


def figure_par(n):
    rid, name = IMG[n]
    src = FIG / f"Figure{n}.png"; im = Image.open(src)
    w = 2400; h = int(im.height * w / im.width)
    im.convert("RGB").resize((w, h), Image.LANCZOS).save(DST / "word" / "media" / name, optimize=True)
    cx = 5715000; cy = int(cx * h / w)
    if cy > 7600000:   # keep within the page height (about 8.3 in)
        cy = 7600000; cx = int(cy * w / h)
    docpr[0] += 1; i = docpr[0]
    return (f'<w:p><w:pPr><w:keepNext/><w:jc w:val="center"/></w:pPr><w:r><w:rPr><w:noProof/></w:rPr><w:drawing>'
            f'<wp:inline distT="0" distB="0" distL="0" distR="0"><wp:extent cx="{cx}" cy="{cy}"/><wp:effectExtent l="0" t="0" r="0" b="0"/>'
            f'<wp:docPr id="{i}" name="Figure{n}.png"/><wp:cNvGraphicFramePr/><a:graphic xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">'
            f'<a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:pic xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
            f'<pic:nvPicPr><pic:cNvPr id="{i}" name="Figure{n}.png"/><pic:cNvPicPr preferRelativeResize="0"/></pic:nvPicPr><pic:blipFill><a:blip r:embed="{rid}"/>'
            f'<a:stretch><a:fillRect/></a:stretch></pic:blipFill><pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{cx}" cy="{cy}"/></a:xfrm>'
            f'<a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:ln/></pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline></w:drawing></w:r></w:p>')


def placeholder_fig(text):
    bd = '<w:pBdr><w:top w:val="single" w:sz="4" w:space="8" w:color="999999"/><w:left w:val="single" w:sz="4" w:space="8" w:color="999999"/><w:bottom w:val="single" w:sz="4" w:space="8" w:color="999999"/><w:right w:val="single" w:sz="4" w:space="8" w:color="999999"/></w:pBdr>'
    return f'<w:p><w:pPr><w:keepNext/>{bd}<w:spacing w:before="240" w:after="240"/><w:jc w:val="center"/></w:pPr>{runs(text)}</w:p>'


# ------------------------------------------------------------- equations
EQ1, EQ2 = accept(P[19]), accept(P[24])
CM = '<w:rPr><w:rFonts w:ascii="Cambria Math" w:eastAsia="Cambria Math" w:hAnsi="Cambria Math" w:cs="Cambria Math"/><w:color w:val="000000"/></w:rPr>'
def mr(t): return f'<m:r>{CM}<m:t xml:space="preserve">{esc(t)}</m:t></m:r>'
def msub(b, s): return f'<m:sSub><m:e>{mr(b)}</m:e><m:sub>{mr(s)}</m:sub></m:sSub>'
ppr_eq = re.search(r'<w:pPr>.*?</w:pPr>', EQ1, flags=re.S).group(0)
EQ3 = (f'<w:p>{ppr_eq}<m:oMath>{mr("T=−3.69")}<m:d><m:e>{mr("h−")}{msub("h", "f")}</m:e></m:d>{mr("−1.11 ")}{msub("h", "f")}</m:oMath>'
       f'<w:r>{CM}<w:tab/><w:t>(3)</w:t></w:r></w:p>')

# ------------------------------------------------------------- tables
TBLPR = re.search(r'<w:tblPr>.*?</w:tblPr>', P[8], flags=re.S).group(0)


def md_table(lines, widths):
    rows = [[c.strip() for c in l.strip().strip("|").split("|")] for l in lines if not re.match(r'^\|\s*-', l)]
    grid = "".join(f'<w:gridCol w:w="{w}"/>' for w in widths)
    tbl = f'<w:tbl>{TBLPR.replace("10095", str(sum(widths)))}<w:tblGrid>{grid}</w:tblGrid>'
    for k, r in enumerate(rows):
        tbl += '<w:tr>' + ("<w:trPr><w:tblHeader/></w:trPr>" if k == 0 else "")
        for c, w in zip(r, widths):
            rpr = ("<w:b/><w:bCs/>" if k == 0 else "") + SZ
            tbl += (f'<w:tc><w:tcPr><w:tcW w:w="{w}" w:type="dxa"/><w:vAlign w:val="center"/></w:tcPr>'
                    f'<w:p><w:pPr><w:spacing w:before="20" w:after="20"/>{"<w:jc w:val=\"center\"/>" if k == 0 else "<w:jc w:val=\"left\"/>"}</w:pPr>{runs(c, rpr)}</w:p></w:tc>')
        tbl += '</w:tr>'
    return tbl + '</w:tbl>'


# ------------------------------------------------------------- parse the markdown
md = MD.read_text(encoding="utf8")
md = re.sub(r'<!--.*?-->', '', md, flags=re.S)
title = re.match(r'\*\*(.+?)\*\*', md).group(1)
abstract = re.search(r'\*\*Abstract\.\*\*\s*(.+)', md).group(1).strip()
main = md[md.index("# 1 Introduction"):md.index("# References added or changed")]
blocks = [b.strip("\n") for b in re.split(r'\n\s*\n', main) if b.strip()]

out = []
# front matter from v2.5: chapter line, title (replaced), authors, affiliations, correspondence
out.append(accept(P[0]))
t1 = accept(P[1]); t1 = re.sub(r'(<w:t[^>]*>)[^<]*(</w:t>)', lambda m: m.group(1) + esc(title) + m.group(2), t1, count=1)
t1 = re.sub(r'(</w:t></w:r>)(<w:r>.*?</w:r>)+(?=</w:p>)', r'\1', t1, flags=re.S)
out.append(t1)
out += [accept(p) for p in P[2:6]]
out.append(f'<w:p>{run("Abstract.", "<w:b/>")}{runs(" " + abstract)}</w:p>')
# product table (OSR template), Product 3 description updated
tb = accept(P[8]).replace("Drake Passage transport time series - in-situ observations",
                          "Drake Passage ADCP velocity, EKE and transport (LM Gould) - in-situ observations")
out += [accept(P[7]), tb, accept(P[9]), accept(P[10])]

i = 0
while i < len(blocks):
    b = blocks[i]
    if b.startswith("## "): out.append(heading(b[3:], 2))
    elif b.startswith("# "): out.append(heading(b[2:], 1))
    elif b == "[EQUATION 1]": out.append(EQ1)
    elif b == "[EQUATION 2]": out.append(EQ2)
    elif "[EQUATION 3]" in b:
        out.append(EQ3)
        rest_ = b.split("[EQUATION 3]", 1)[1].strip()
        if rest_: out.append(para(rest_, '<w:jc w:val="right"/>', SZ))
    elif b.startswith("**Figure"):
        m = re.match(r'\*\*Figure (\d)\.\*\*\s*(.*)', b, flags=re.S); n = int(m.group(1))
        if n in IMG: out.append(figure_par(n))
        else: out.append(placeholder_fig(f"[PENDING: Figure {n} — Antarctic-slope sea-level index and model-calibrated transport; requires the southern altimetry extraction and the bathymetry to 80°S]"))
        out.append(caption(f"Figure {n}.", m.group(2)))
    elif b.startswith("**Table 2.**"):
        cap = b[len("**Table 2.**"):].strip()
        tl = blocks[i + 1].splitlines(); i += 1
        out.append(md_table(tl, [2700, 2000, 2700, 2600]))
        out.append(caption("Table 2.", cap))
    else:
        out.append(para(b.replace("\n", " ")))
    i += 1

# ------------------------------------------------------------- references: keep cited ones, add new ones
refs = [accept(p) for p in P[80:128]]
corpus = re.sub(r'\s+', ' ', main + " " + txt(tb))
NEW = {
    "Hughes": ('Hughes, C. W., Williams, J., Coward, A. C., and de Cuevas, B. A.: Antarctic circumpolar transport and the southern mode: '
               'a model investigation of interannual to decadal timescales, Ocean Sci., 10, 215–225, https://doi.org/10.5194/os-10-215-2014, 2014. [verify]'),
    "Marshall, D": ('Marshall, D. P., Ambaum, M. H. P., Maddison, J. R., Munday, D. R., and Novak, L.: Eddy saturation and frictional control '
                    'of the Antarctic Circumpolar Current, Geophys. Res. Lett., 44, 286–292, https://doi.org/10.1002/2016GL071702, 2017. [verify]'),
}
kept, dropped = [], []
for r in refs:
    t = txt(r); first = re.match(r'([^,]+),', t).group(1).strip(); yr = re.findall(r'(\d{4})[a-z]?\.?\s*$', t.strip()) or re.findall(r', (\d{4})[a-z]?\.', t)
    yr = yr[-1] if yr else ""
    two = re.match(r'([^,]+), [^,]+?, and ([^,:]+?),', t)
    pat = re.escape(first.split()[-1] if first.split()[-1] != "Naveira" else first) + r'[^()]{0,80}?' + yr
    hit = re.search(pat, corpus) or re.search(re.escape(first) + r'.{0,120}?' + yr, corpus)
    (kept if hit else dropped).append((t[:60], r))
def key(t): return t.lower()
entries = [(t, r) for t, r in kept] + [(v, para(v)) for v in NEW.values()]
entries.sort(key=lambda e: key(re.sub(r'[ÁÀÉÍÓÚáéíóúñçü]', lambda m: "aaeiouaeiouncu"["ÁÀÉÍÓÚáéíóúñçü".index(m.group(0))], e[0])))
out.append(heading("References", 1))
out += [r for _, r in entries]
print("references kept", len(kept), "dropped:", [d for d, _ in dropped])

sect = P[-1]
newbody = "".join(out) + sect
(DST / "word" / "document.xml").write_text(head + "<w:body>" + newbody + "</w:body>" + tail, encoding="utf8")
rp = DST / "word" / "_rels" / "document.xml.rels"; rr = rp.read_text(encoding="utf8")
rr = re.sub(r'<Relationship Id="rId15"[^>]*/>', '', rr); rp.write_text(rr, encoding="utf8")
(DST / "word" / "media" / "image4.png").unlink(missing_ok=True)
print("paragraph blocks", len(out))
