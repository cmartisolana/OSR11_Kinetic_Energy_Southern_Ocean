"""Append a new row of panels (fig2gh_v341.py, same width) below the Drake figure PDF (v3.35 composite), unchanged above.
Usage: python fig2_append_v341.py <Figure2_drake.pdf> <fig2gh.pdf> <out.pdf>"""
import sys, pikepdf
src, row, dst = sys.argv[1:4]
A = pikepdf.open(src); B = pikepdf.open(row)
W, H = [float(v) for v in A.pages[0].MediaBox[2:]]; h = float(B.pages[0].MediaBox[3])
out = pikepdf.new(); out.add_blank_page(page_size=(W, H + h)); pg = out.pages[0]
pg.add_overlay(out.copy_foreign(A.pages[0]).as_form_xobject() if False else A.pages[0], pikepdf.Rectangle(0, h, W, H + h))
pg.add_overlay(B.pages[0], pikepdf.Rectangle(0, 0, W, h))
out.save(dst); print("ok", W, H + h)
