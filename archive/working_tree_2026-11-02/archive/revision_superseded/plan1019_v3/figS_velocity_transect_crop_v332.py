"""SI figure: the original Fig. 2e comparison (ADCP vs altimetry, cross-track velocity trends) cropped from the
v2.0 Drake Passage figure PDF, with the panel letter removed."""
import sys, pikepdf
src, dst = sys.argv[1:3]
pdf = pikepdf.open(src); pg = pdf.pages[0]
ops = list(pikepdf.parse_content_stream(pg))
i = 16163; assert str(ops[i].operator) == 'cm' and abs(float(ops[i].operands[4]) - 13.676) < 0.01
j = next(k for k in range(i, i + 6) if str(ops[k].operator) in ('TJ', 'Tj'))
ops.pop(j)
pg.Contents = pdf.make_stream(pikepdf.unparse_content_stream(ops))
pg.MediaBox = pikepdf.Array([8, 0, 245.6, 177]); pg.CropBox = pg.MediaBox
pdf.save(dst); print("ok")
