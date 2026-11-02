"""Accept all tracked changes in a document.xml (lxml), keeping OMML intact."""
import sys
from lxml import etree
W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
q = lambda t: f"{{{W}}}{t}"
src, dst = sys.argv[1], sys.argv[2]
tree = etree.parse(src); root = tree.getroot()
# 1. remove deleted content wrappers; unwrap insertions
for el in list(root.iter(q("del"))):
    par = el.getparent()
    if par is not None and par.tag != q("rPr"):
        par.remove(el)
for el in list(root.iter(q("ins"))):
    par = el.getparent()
    if par is not None and par.tag != q("rPr"):
        i = par.index(el)
        for ch in list(el):
            par.insert(i, ch); i += 1
        par.remove(el)
for el in list(root.iter(q("rPrChange"))) + list(root.iter(q("pPrChange"))):
    el.getparent().remove(el)
# 2. paragraph-mark markers
for p in list(root.iter(q("p"))):
    ppr = p.find(q("pPr"))
    rpr = ppr.find(q("rPr")) if ppr is not None else None
    if rpr is None:
        continue
    ins = rpr.find(q("ins"))
    if ins is not None:
        rpr.remove(ins)
    dl = rpr.find(q("del"))
    if dl is not None:
        rpr.remove(dl)
        nxt = p.getnext()
        while nxt is not None and nxt.tag != q("p"):
            nxt = nxt.getnext()
        if nxt is None:
            continue
        content = [c for c in p if c.tag != q("pPr")]
        nppr = nxt.find(q("pPr"))
        pos = 1 if nppr is not None else 0
        for c in content:
            nxt.insert(pos, c); pos += 1
        p.getparent().remove(p)
tree.write(dst, xml_declaration=True, encoding="UTF-8", standalone=True)
