"""Probe (19 Oct 2026): matched-envelope half-power latitude averaged over longitude windows around Drake Passage,
from outputs/a2/tracks_baseline.npz (baseline 200 cm2/s2, 25 %, 100 km; keep = matched in >= 90 % of months).
Compares noise and trend with the single LM Gould section (dk2) and regresses annual means on wind stress."""
import sys, numpy as np, pandas as pd, xarray as xr
REV = "/mnt/user-data/uploads/repository/revision"; sys.path.insert(0, REV)
import rstats
TB = np.load(f"{REV}/outputs/a2/tracks_baseline.npz"); lon = xr.open_dataset(f"{REV}/outputs/clim_allsat.nc").longitude.values
jl = lon[TB["j"]]; hp = TB["lat_hp"].astype(float); keep = TB["keep"]; acc = TB["acc"]
t = 2001 + (np.arange(300) + 0.5) / 12; mon = np.tile(np.arange(1, 13), 25)
S = pd.read_csv("drake/out/dk1/dk1_series.csv")
rows = []
for name, (l0, l1), sel_all in [("passage 70-55W, matched, both fronts", (-70, -55), False), ("passage 70-55W, all ACC envelopes >=50% months", (-70, -55), True),
                                ("upstream 85-70W, matched", (-85, -70), False), ("downstream 55-45W, matched", (-55, -45), False)]:
    m = (jl >= l0) & (jl < l1) & (acc if sel_all else keep) & (TB["lat_hp_clim"] < -52)   # exclude Malvinas-Current envelopes
    if sel_all: m &= np.isfinite(hp).mean(0) >= 0.5
    if m.sum() == 0: rows.append(dict(region=name, n=0)); continue
    a = hp[:, m]; y = np.nanmean(a - np.nanmean(a, 0, keepdims=True), 1) * 111.195
    T = rstats.trend_bcz(t, y, month=mon)
    yr = pd.Series(y).groupby(np.repeat(np.arange(2001, 2026), 12)).mean()
    out = dict(region=name, n=int(m.sum()), cover=float(np.isfinite(a).mean()), std_month_km=float(np.nanstd(y)), std_annual_km=float(yr.std()),
               trend_km_dec=T["slope"] * 10, lo=T["ci_lo"] * 10, hi=T["ci_hi"] * 10, p=T["p"])
    for w in ("wind_up_TAUX", "wind_box_TAUX"):
        wa = S.groupby("year")[w].mean(); wp = 100 * (wa / wa.mean() - 1)
        b = np.polyfit(wp.values, yr.values, 1)[0]; r = np.corrcoef(wp.values, yr.values)[0, 1]
        out[f"r_{w}"] = r; out[f"km_per_pct_{w}"] = b
    rows.append(out)
R = pd.DataFrame(rows); pd.set_option("display.width", 250); print(R.round(3).to_string())
R.to_csv("drake/out/drake_hp_regional_probe.csv", index=False)
