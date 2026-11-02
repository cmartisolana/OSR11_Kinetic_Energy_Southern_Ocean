"""dk32a: extract the Southern Ocean (south of 35 S) from the Roemmich-Gilson Argo climatology (netCDF-3 classic) with
numpy only (minimal netCDF-3 reader; no netCDF library needed). Writes per variable: coordinates, the 2004-2018 mean,
the masks, and the monthly anomalies split into one .npy per year (float32), small enough to transfer.

Purpose : convert the Roemmich-Gilson Argo climatology (Product 5) to small .npy files south of 35 S.
Inputs  : RG_ArgoClim_Temperature_2019.nc and RG_ArgoClim_Salinity_2019.nc (https://sio-argo.ucsd.edu/RG_Climatology.html).
Outputs : $ARGO_RG_DIR/so_T/*.npy and $ARGO_RG_DIR/so_S/*.npy (default $OSR11_DATA/argo_rg); the RG monthly extension files
          (RG_ArgoClim_20*.nc, 2019-2025) are placed unchanged in $ARGO_RG_DIR/ and read there by 04a and 04e.
Run     : python scripts/00_extract/00d_argo_rg_extract.py <RG_ArgoClim_Temperature_2019.nc> $OSR11_DATA/argo_rg/so_T
          python scripts/00_extract/00d_argo_rg_extract.py <RG_ArgoClim_Salinity_2019.nc>    $OSR11_DATA/argo_rg/so_S
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import sys, struct, os
import numpy as np
F, OUT = sys.argv[1], sys.argv[2]; LATMAX = float(sys.argv[3]) if len(sys.argv) > 3 else -35.0
os.makedirs(OUT, exist_ok=True)
TY = {1: ("b", 1), 2: ("c", 1), 3: (">i2", 2), 4: (">i4", 4), 5: (">f4", 4), 6: (">f8", 8)}
with open(F, "rb") as f: hdr = f.read(200000)
pos = [0]
def rd(n): b = hdr[pos[0]:pos[0] + n]; pos[0] += n; return b
def i4(): return struct.unpack(">i", rd(4))[0]
def name(): n = i4(); s = rd(n).decode(); rd((4 - n % 4) % 4); return s
def atts():
    tag, n = i4(), i4(); out = {}
    for _ in range(n):
        k = name(); t, m = i4(), i4(); sz = TY[t][1] * m; v = rd(sz); rd((4 - sz % 4) % 4)
        out[k] = v.decode(errors="ignore") if t == 2 else np.frombuffer(v, TY[t][0]).tolist()
    return out
magic = rd(4); assert magic[:3] == b"CDF", magic; off8 = magic[3] == 2
nrec = i4(); tag, nd = i4(), i4(); dims = [(name(), i4()) for _ in range(nd)]
gat = atts(); tag, nv = i4(), i4(); V = {}
for _ in range(nv):
    nm = name(); nid = i4(); did = [i4() for _ in range(nid)]; va = atts(); t = i4(); vs = i4()
    beg = struct.unpack(">q", rd(8))[0] if off8 else i4()
    V[nm] = dict(dims=did, t=t, vsize=vs, begin=beg, att=va)
recdim = [i for i, (n, l) in enumerate(dims) if l == 0]
recsize = sum(v["vsize"] for v in V.values() if v["dims"] and v["dims"][0] in recdim)
print("dims", dims, "nrec", nrec); [print(k, [dims[d][0] for d in v["dims"]], v["att"].get("units", "")) for k, v in V.items()]
mm = np.memmap(F, dtype="u1", mode="r")
def get(k):
    v = V[k]; dt = np.dtype(TY[v["t"]][0]); shp = [dims[d][1] for d in v["dims"]]
    if v["dims"] and v["dims"][0] in recdim:
        n1 = int(np.prod(shp[1:])) * dt.itemsize
        def rec(r): return np.frombuffer(mm[v["begin"] + r * recsize: v["begin"] + r * recsize + n1], dt).reshape(shp[1:])
        return rec
    n = int(np.prod(shp)) * dt.itemsize
    return np.frombuffer(mm[v["begin"]: v["begin"] + n], dt).reshape(shp)
latk = [k for k in V if k.upper().startswith("LAT")][0]; lat = get(latk).astype(float); js = np.where(lat <= LATMAX)[0]
latdim = V[latk]["dims"][0]
meta = {}
for k, v in V.items():
    g = get(k)
    if callable(g):
        ax = [d for d in v["dims"][1:]].index(latdim) if latdim in v["dims"][1:] else None
        for y0 in range(0, nrec, 12):
            a = np.array([g(r) for r in range(y0, min(y0 + 12, nrec))])
            if ax is not None: a = np.take(a, js, axis=ax + 1)
            np.save(os.path.join(OUT, f"{k}_m{y0:03d}.npy"), a.astype(np.float32))
        meta[k] = "record"
    else:
        a = g
        if latdim in v["dims"]: a = np.take(a, js, axis=v["dims"].index(latdim))
        np.save(os.path.join(OUT, f"{k}.npy"), a.astype(np.float32) if a.dtype.kind == "f" else a)
        meta[k] = [dims[d][0] for d in v["dims"]]
with open(os.path.join(OUT, "meta.txt"), "w") as f:
    f.write(f"source {F}\nlat <= {LATMAX}\nnrec {nrec}\n"); [f.write(f"{k} {m} {V[k]['att']}\n") for k, m in meta.items()]; f.write(f"global {gat}\n")
print("done", OUT)
