"""Subset the AVISO META (4.0 or 3.2) DT allsat eddy atlas to the Southern Ocean and 2001 onwards.

Input : META4_DT_allsat_<anticyclonic|cyclonic>_*.nc (one file per polarity; META3.2 long/short files also work).
        Observation-based layout: one record per eddy observation (track, observation_number, time, latitude,
        longitude, amplitude, speed_radius, speed_average, effective_radius, ...). Variables are packed
        integers with scale_factor/add_offset; they are copied packed, with their attributes.
Output: <out_dir>/eddies_so_<polarity>_<kind>.nc with records at latitude <= -35 and time >= 2001-01-01.
Usage: python eddies_subset_aviso.py <atlas_dir> <out_dir>     (about 60 s and 1 GB of memory per file)
"""
import sys, glob, os
import numpy as np
from netCDF4 import Dataset

KEEP = ["track", "observation_number", "observation_flag", "time", "latitude", "longitude",
        "amplitude", "speed_radius", "speed_average", "effective_radius", "effective_area",
        "latitude_max", "longitude_max", "num_point_e", "num_point_s", "segment"]
atlas_dir, out_dir = sys.argv[1], sys.argv[2]
os.makedirs(out_dir, exist_ok=True)


def unpacked(var, raw):
    sf = var.getncattr("scale_factor") if "scale_factor" in var.ncattrs() else 1.0
    ao = var.getncattr("add_offset") if "add_offset" in var.ncattrs() else 0.0
    return raw.astype(np.float64) * sf + ao


for f in sorted(glob.glob(os.path.join(atlas_dir, "META*_DT_allsat_*.nc"))):
    base = os.path.basename(f)
    pol = "anticyclonic" if "nticyclonic" in base else "cyclonic"
    kind = "long" if "_long_" in base else "short" if "_short_" in base else "all"
    with Dataset(f) as ds:
        ds.set_auto_maskandscale(False)
        lat = unpacked(ds.variables["latitude"], ds.variables["latitude"][:])
        tvar = ds.variables["time"]
        assert tvar.units.startswith("days since 1950-01-01"), tvar.units
        tdays = unpacked(tvar, tvar[:])
        sel = np.nonzero((lat <= -35.0) & (tdays >= 18628.0))[0]      # 18628 = 2001-01-01
        del lat, tdays
        print(base, "records", ds.dimensions["obs"].size, "-> kept", sel.size, flush=True)
        out = os.path.join(out_dir, f"eddies_so_{pol}_{kind}.nc")
        with Dataset(out, "w") as o:
            o.set_auto_maskandscale(False)
            o.createDimension("obs", sel.size)
            o.source = base; o.subset = "latitude <= -35, time >= 2001-01-01"
            for v in KEEP:
                if v not in ds.variables:
                    continue
                src = ds.variables[v]
                fv = src.getncattr("_FillValue") if "_FillValue" in src.ncattrs() else None
                dst = o.createVariable(v, src.dtype, ("obs",), zlib=True, complevel=4, fill_value=fv)
                dst.set_auto_maskandscale(False)     # write the packed integers as they are
                for a in src.ncattrs():
                    if a != "_FillValue":
                        dst.setncattr(a, src.getncattr(a))
                dst[:] = src[:][sel]
                print("  ", v, flush=True)
print("done")
