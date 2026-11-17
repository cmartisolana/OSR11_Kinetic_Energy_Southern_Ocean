"""Compare two sets of rendered figures (PNG), with a tolerance for the rasterisation of text.

Purpose : check re-rendered figures against the figures of the manuscript (rendered/png). Fonts and their rasterisation differ
          from one system to another, so labels can move by a pixel; data marks do not. A pixel therefore counts as different only
          when its colour is outside the range of colours found within +-R pixels in the other image (in both directions).
Inputs  : two folders of PNG files with the same names
Outputs : one line per figure: number of pixels that differ exactly, and number beyond the tolerance with their bounding box
Run     : python scripts/compare_figures.py figures_out/png rendered/png [R=3]
"""
import sys
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import maximum_filter, minimum_filter

Image.MAX_IMAGE_PIXELS = None
A, B = Path(sys.argv[1]), Path(sys.argv[2]); R = int(sys.argv[3]) if len(sys.argv) > 3 else 3; T = 24


def far(x, y):
    k = (2 * R + 1, 2 * R + 1, 1)
    return ((x < minimum_filter(y, size=k) - T) | (x > maximum_filter(y, size=k) + T)).any(2)


for f in sorted(A.glob("*.png")):
    g = B / f.name
    if not g.exists():
        print(f"{f.name}: not in {B}"); continue
    x = np.asarray(Image.open(f).convert("RGB"), dtype=np.int16); y = np.asarray(Image.open(g).convert("RGB"), dtype=np.int16)
    if x.shape != y.shape:
        print(f"{f.name}: sizes differ {x.shape[:2]} vs {y.shape[:2]}"); continue
    exact = int((np.abs(x - y).max(2) > 0).sum()); d = far(x, y) | far(y, x); n = int(d.sum())
    msg = f"{f.name}: {'identical' if exact == 0 else f'{exact} pixels differ ({100 * exact / d.size:.2f} %)'}"
    if n:
        ys, xs = np.where(d); msg += f"; {n} beyond +-{R} px (x {xs.min()}-{xs.max()}, y {ys.min()}-{ys.max()}: check that these are labels)"
    elif exact:
        msg += f"; none beyond +-{R} px"
    print(msg)
