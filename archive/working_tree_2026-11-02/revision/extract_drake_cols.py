"""Extract Drake-region (90-35W) columns of adt_mean (and ke if present) from staged monthly files, then delete them."""
import glob, os, numpy as np, xarray as xr
from pathlib import Path
SRC='/mnt/user-data/uploads/OSR11/repository/revision/data/monthly_allsat'
OUT=Path('/home/claude/x1drake/cols'); OUT.mkdir(exist_ok=True)
for f in sorted(glob.glob(f'{SRC}/*/*.nc')):
    y=Path(f).parent.name; m=Path(f).stem
    with xr.open_dataset(f) as d:
        sub=d.sel(longitude=slice(-90,-35))
        np.savez_compressed(OUT/f'{y}{m}.npz', lat=sub.latitude.values, lon=sub.longitude.values,
            adt=sub.adt_mean.values.astype('float32'), **({'ke':sub.ke.values.astype('float32')} if 'ke' in sub else {}))
    os.remove(f)
print(len(list(OUT.glob('*.npz'))),'months extracted')
