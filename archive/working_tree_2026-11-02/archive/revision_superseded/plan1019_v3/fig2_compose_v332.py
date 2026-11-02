"""Compose Fig. 2 (v3.32): panels a-d from the v2.0 Drake figure PDF (unchanged), old panels e-f removed, new panels
e-f (fig2ef_v332.py) overlaid. Usage: python fig2_compose_v332.py <v2.0 Figure1.pdf> <fig2ef.pdf> <out.pdf>"""
import sys, pikepdf
src, ef, dst = sys.argv[1:4]
pdf = pikepdf.open(src); pg = pdf.pages[0]
ops = list(pikepdf.parse_content_stream(pg))
# 15024-15035: white background patch of old panel e (end of the panel-d block); 15037-end: old panels e and f
assert str(ops[15024].operator) == 'gs' and str(ops[15035].operator) == 'f' and str(ops[15036].operator) == 'Q'
assert str(ops[15037].operator) == 'q' and abs(float(ops[15030].operands[1]) - 23.170957) < 1e-3
keep = ops[:15024] + ops[15036:15037]
pg.Contents = pdf.make_stream(pikepdf.unparse_content_stream(keep))
new = pikepdf.open(ef)
pg.add_overlay(new.pages[0], pikepdf.Rectangle(*[float(v) for v in pg.MediaBox]))
pdf.save(dst); print("ok")
