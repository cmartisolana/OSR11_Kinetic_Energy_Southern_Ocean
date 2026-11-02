"""Edit Figure 2 (Drake Passage) vector PDF in place: remove the in situ ADCP series from panel (e),
open markers on the 2005-2019 altimetry line, legend compacted, panel (e) title changed. Everything else untouched."""
import sys, pikepdf
from pikepdf import Operator as Op, String, Name
src, dst = sys.argv[1:3]
pdf = pikepdf.open(src); pg = pdf.pages[0]
ops = list(pikepdf.parse_content_stream(pg))
D = lambda i: str(ops[i].operator)
ORANGE = (0.902, 0.624, 0.0)
def col(i): return tuple(round(float(x), 3) for x in ops[i].operands)
# sanity checks on the block layout found by inspection
assert D(16005) == 'RG' and col(16005) == ORANGE and D(15998) == 'q' and D(16041) == 'Q'      # ADCP line
assert D(16179) == 'q' and D(16195) == 'Q' and col(16184) == ORANGE                         # ADCP filled markers
assert D(16196) == 'q' and D(16269) == 'Q' and col(16203) == ORANGE                         # ADCP open markers
assert D(16270) == 'q' and D(16292) == 'Q'                                                  # 2005-2019 altimetry filled markers
assert D(16474) == 'RG' and col(16474) == ORANGE                                            # legend entry 1
assert D(16175) == 'TJ'
new = {}
drop = set(range(15998, 16042)) | set(range(16179, 16270)) | set(range(16474, 16480)) | set(range(16485, 16493))
# filled -> open markers for the 2005-2019 altimetry line (same as the open-marker style: white fill, 0.6 w, M5)
for i in range(16270, 16293):
    if D(i) == 'rg': new[i] = [pikepdf.ContentStreamInstruction([1], Op('g'))]
    elif D(i) == 'RG': new[i] = [ops[i], pikepdf.ContentStreamInstruction([0.6], Op('w'))]
    elif D(i) == 'Do':
        assert str(ops[i].operands[0]) == '/M4'; new[i] = [pikepdf.ContentStreamInstruction([Name('/M5')], Op('Do'))]
# legend: move entries 2-3 up one row and raise the frame bottom by the same amount
dy = 150.644525 - 142.013666
def shift(i, idxs):
    o = list(ops[i].operands)
    for k in idxs: o[k] = float(o[k]) + dy
    new[i] = [pikepdf.ContentStreamInstruction(o, ops[i].operator)]
for i in (16497, 16498, 16499, 16518, 16519, 16520): shift(i, [1])
for i in (16507, 16528): shift(i, [5])
for i in (16459, 16460, 16466): shift(i, [1])
for i in (16461, 16467): shift(i, [1, 3, 5])
# panel (e) title (F1 subset: every character used here is in the embedded subset)
title = "Cross-track velocity trend along the LM Gould transect"
cmap = pg.Resources.Font['/F1'].ToUnicode.read_bytes().decode()
assert all(f"<{ord(c):04x}>" in cmap for c in set(title)), [c for c in set(title) if f"<{ord(c):04x}>" not in cmap]
new[16175] = [pikepdf.ContentStreamInstruction([String(title.encode('utf-16-be'))], Op('Tj'))]
out = []
for i, ins in enumerate(ops):
    if i in drop: continue
    out.extend(new.get(i, [ins]))
pg.Contents = pdf.make_stream(pikepdf.unparse_content_stream(out))
pdf.save(dst)
print("ok", len(ops), "->", len(out))
