"""SI figure (v3.32): the former Fig. 2e-f of the v2.0 Drake figure (comparison along the LM Gould transect with
Gutierrez-Villanueva et al. 2023), cropped from the v2.0 PDF with the panel letters e, f relabelled a, b.
Usage: python figS_gv_transect_v332.py <v2.0 Figure1.pdf> <out.pdf>"""
import sys, pikepdf
src, dst = sys.argv[1:3]
pdf = pikepdf.open(src); pg = pdf.pages[0]
ops = list(pikepdf.parse_content_stream(pg))
def relabel(i_cm, x, ch):
    assert str(ops[i_cm].operator) == 'cm' and abs(float(ops[i_cm].operands[4]) - x) < 0.01
    j = next(k for k in range(i_cm, i_cm + 6) if str(ops[k].operator) in ('TJ', 'Tj'))
    ops[j] = pikepdf.ContentStreamInstruction([pikepdf.String(ch.encode('utf-16-be'))], pikepdf.Operator('Tj'))
cm_f = next(i for i in range(15509, 15531) if str(ops[i].operator) == 'cm')
relabel(16163, 13.676328, "a"); relabel(cm_f, 246.584366, "b")
pg.Contents = pdf.make_stream(pikepdf.unparse_content_stream(ops))
pg.MediaBox = pikepdf.Array([5, 0, 466, 177]); pg.CropBox = pg.MediaBox
pdf.save(dst); print("ok")
