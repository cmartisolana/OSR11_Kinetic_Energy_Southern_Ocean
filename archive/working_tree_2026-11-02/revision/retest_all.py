"""Re-run every reported trend with the bias-corrected, calibrated test (rstats.trend_bcz).

Reads the saved series of each analysis (no field is recomputed) and writes one row per test to
outputs/retest/retest_<family>.csv, with the previously reported p for comparison. Rows already
present are skipped, so the script can be restarted.
Usage: python retest_all.py <family> [<family> ...]   (families: A2 A5 A6 A7 B1 LIU CONTRAST)
       env RETEST_N_OUTER (default 300) = number of calibration nulls (p floor 1/(n_outer+1)).
"""
import os as _os
from pathlib import Path as _P
_REV = str(_P(__file__).resolve().parents[0])  # revision/
_ROOT_OUT = _REV + '/../outputs'; _REV_OUT = _REV + '/outputs'; _DK = _REV + '/plan1019/drake'; _DK_OUT = _DK + '/out'
_CACHE = _os.environ.get('OSR11_CACHE', _REV + '/cache'); _GV = _os.environ.get('GV23_DIR', _REV + '/data/TrendsDrakePassageTransp')
import json
import os
import sys
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rstats as R  # noqa: E402

O = HERE / "outputs"
RT = O / "retest"; RT.mkdir(exist_ok=True)
N_OUTER = int(os.environ.get("RETEST_N_OUTER", 300))
TRENDS_DIR = Path((_ROOT_OUT + "/trends"))


class Runner:
    def __init__(self, family):
        self.f = RT / f"retest_{family}.csv"
        self.done = set(pd.read_csv(self.f).key) if self.f.exists() else set()
        self.family = family

    def __call__(self, key, x, y, month=None, monthly=True, old_p=np.nan, old_test="", scale=1.0, **meta):
        if key in self.done:
            return
        x = np.asarray(x, float); y = np.asarray(y, float)
        o = R.trend_bcz(x, y, month=month, monthly=monthly, n_outer=N_OUTER, seed=0)
        m = float(np.nanmean(y))
        row = dict(key=key, family=self.family, slope=o["slope"] * scale, ci_lo=o["ci_lo"] * scale, ci_hi=o["ci_hi"] * scale,
                   p_new=o["p"], p_aic=o.get("p_aic", np.nan), p_ar12=o.get("p_ar12", np.nan), z=o["z"], z_crit=o["z_crit"],
                   old_p=old_p, old_test=old_test,
                   pct_per_decade=1000 * o["slope"] / m if m else np.nan, n=o["n"], ar_order=o["ar_order"],
                   phi_sum=o["phi_sum"], phi_sum_bc=o["phi_sum_bc"], meta=json.dumps(meta, default=str))
        pd.DataFrame([row]).to_csv(self.f, mode="a", header=not self.f.exists(), index=False)
        self.done.add(key)
        print(self.family, key, round(row["p_new"], 4), "old", old_p if not np.isfinite(old_p) else round(old_p, 4), flush=True)


def fam_A5():
    run = Runner("A5")
    for p in ("allsat", "twosat"):
        S = pd.read_csv(O / "a5" / f"a5_series_{p}.csv")
        T = pd.read_csv(O / "a5" / f"a5_trends_{p}.csv")
        old = {(r.series, r.period): (r.p_cal if np.isfinite(r.p_cal) else r.p_boot, "p_cal" if np.isfinite(r.p_cal) else "sieve")
               for r in T.itertuples() if isinstance(r.series, str)}
        cols = [c for c in S.columns if "_" in c and c.split("_")[-1] in ("KE", "EKE", "LF", "SUB", "PROJ", "SPD", "EKE3")]
        # headline first
        cols = sorted(cols, key=lambda c: (not c.startswith(("ACC_", "Pacific_", "Atlantic_", "Indian_")), c))
        for c in cols:
            for per in ("2001-2025", "2001-2020"):
                if per == "2001-2020" and not c.startswith("box35_65"):
                    continue
                sel = np.ones(len(S), bool) if per == "2001-2025" else S.year.values <= 2020
                op, ot = old.get((c, per), (np.nan, ""))
                run(f"{p}|{c}|{per}", S.x.values[sel], S[c].values[sel], month=S.month.values[sel], old_p=op, old_test=ot,
                    product=p, series=c, period=per)
        B = pd.read_csv(O / "a5" / f"a5_blockKE_{p}.csv")
        bb = B[B.block == "annual"]
        for c in [c for c in bb.columns if c.endswith("KEblock")]:
            r = T[T.series == f"{c}_annual"]
            run(f"{p}|{c}_annual|2001-2025", bb.x.values, bb[c].values, monthly=False,
                old_p=float(r.p_boot.iloc[0]) if len(r) else np.nan, old_test="sieve", product=p, series=f"{c}_annual", period="2001-2025")


def fam_A2():
    run = Runner("A2")
    os.environ.setdefault("OSR11_CACHE", _CACHE)
    import a2_run as A
    lon = np.load(A.CACHE / "lon.npy")
    t = pd.DataFrame({"year": np.repeat(np.arange(2001, 2026), 12), "month": np.tile(np.arange(1, 13), 25)})
    t["x"] = t.year + (t.month - 0.5) / 12
    S = pd.read_csv(O / "a2" / "a2_summary.csv")
    for name in ["baseline"] + sorted(p.stem.replace("tracks_", "") for p in (O / "a2").glob("tracks_*.npz") if "baseline" not in p.stem):
        d = np.load(O / "a2" / f"tracks_{name}.npz", allow_pickle=True)
        keep = d["keep"]; sec = A.sector(lon[d["j"]])
        groups = {"ACC": np.ones(keep.size, bool), "Atlantic": sec == "Atlantic", "Indian": sec == "Indian", "Pacific": sec == "Pacific"}
        fields = ["lat_hp", "peak_lat", "peak_val", "seg_km", "fwhm_km", "w25_km", "sig_km"] if name == "baseline" else ["lat_hp", "w25_km"]
        for g, gm in groups.items():
            if name != "baseline" and g != "ACC":
                continue
            sel = keep & gm
            for f in fields:
                a = d[f][:, sel].astype(float)
                y = np.nanmean(a - np.nanmean(a, 0, keepdims=True), axis=1) * (1e4 if f == "peak_val" else 1.0)
                r = S[(S.set == name) & (S.region == g)]
                op = float(r[f"m_{f}_p"].iloc[0]) if len(r) else np.nan
                run(f"{name}|{g}|{f}", t.x.values, y, month=t.month.values, old_p=op, old_test="sieve", set=name, region=g, field=f)
    # A2f longitude bins (baseline)
    d = np.load(O / "a2" / "tracks_baseline.npz", allow_pickle=True); keep = d["keep"]; lonq = lon[d["j"]]
    LT = pd.read_csv(O / "a2" / "a2f_lon_trends.csv").set_index("lon0")
    for b0 in np.arange(-180, 180, 5):
        sel = keep & (lonq >= b0) & (lonq < b0 + 5)
        if sel.sum() < 5:
            continue
        a = d["lat_hp"][:, sel].astype(float)
        y = np.nanmean(a - np.nanmean(a, 0, keepdims=True), axis=1)
        run(f"lonbin|{b0}", t.x.values, y, month=t.month.values, old_p=float(LT.loc[b0, "p"]) if b0 in LT.index else np.nan,
            old_test="sieve", lon0=b0)


def fam_A6():
    run = Runner("A6")
    A6 = pd.read_csv(O / "a6" / "a6_trends.csv").set_index("series")
    Tc = pd.read_csv(O / "a6" / "a6c_transport_fixed_series.csv")
    oc = pd.read_csv(O / "a6" / "a6c_trends.csv"); oc = {(r.series, r.period): r.p_cal for r in oc.itertuples()}
    for c in [c for c in Tc.columns if c.startswith("T")]:
        for per in ("2001-2025", "2005-2019"):
            sel = np.ones(len(Tc), bool) if per == "2001-2025" else ((Tc.year >= 2005) & (Tc.year <= 2019)).values
            run(f"{c}|{per}", Tc.x.values[sel], Tc[c].values[sel], month=Tc.month.values[sel], old_p=oc.get((c, per), np.nan),
                old_test="p_cal", series=c, period=per)
    T = pd.read_csv(O / "a6" / "a6_transport_series.csv")
    for per in ("2001-2025", "2005-2019"):
        sel = np.ones(len(T), bool) if per == "2001-2025" else ((T.year >= 2005) & (T.year <= 2019)).values
        run(f"dADT|{per}", T.x.values[sel], T.dADT.values[sel], month=T.month.values[sel],
            old_p=float(A6.loc[f"dADT_{per}", "p_cal"]), old_test="p_cal", series="dADT", period=per)
    obs = pd.read_csv(TRENDS_DIR / "drake_passage_total_transport_observations.csv", parse_dates=["date"]).dropna(subset=["total_transport_sv"])
    ann = obs.groupby(obs.date.dt.year).total_transport_sv.mean()
    run("insitu_annual|2005-2019", ann.index.values.astype(float) + 0.5, ann.values, monthly=False,
        old_p=float(A6.loc["insitu_annual_2005-2019", "p_cal"]), old_test="p_cal", series="insitu_annual", period="2005-2019")
    BX = pd.read_csv(O / "a6" / "a6b_drake_box_series.csv"); BX["x"] = BX.year + (BX.month - 0.5) / 12
    ob = pd.read_csv(O / "a6" / "a6b_drake_box_trends.csv"); ob = {(r.series, r.period): r.p_cal for r in ob.itertuples()}
    for mk in ("ice_any", "nomask"):
        for v in ("EKE", "KE"):
            for per in ("2001-2025", "2001-2024", "2005-2019"):
                sel = {"2001-2025": BX.year <= 2025, "2001-2024": BX.year <= 2024, "2005-2019": (BX.year >= 2005) & (BX.year <= 2019)}[per].values
                run(f"Drake_{v}_{mk}|{per}", BX.x.values[sel], BX[f"{mk}_{v}"].values[sel], month=BX.month.values[sel],
                    old_p=ob.get((f"Drake_{v}_{mk}", per), np.nan), old_test="p_cal", series=f"Drake_{v}_{mk}", period=per)
    W = pd.read_csv(O / "a6" / "a6_wind_stress_series.csv")
    W = W[(W.year >= 2001) & (W.year <= 2025)]
    for c in ("Circumpolar_45_60S", "Pacific_45_60S", "Drake_box"):
        for per in ("2001-2025", "2008-2025"):
            sel = (W.year >= (2001 if per == "2001-2025" else 2008)).values
            run(f"tau_{c}|{per}", W.x.values[sel], W[c].values[sel], month=W.month.values[sel],
                old_p=float(A6.loc[f"tau_{c}_{per}", "p_cal"]), old_test="p_cal", series=f"tau_{c}", period=per)


def fam_A7():
    run = Runner("A7")
    FEAT = ["Kerguelen", "SEIR_Macquarie", "PAR", "Drake"]
    VARS = ["EKE", "KE", "SPD", "SUB", "LF"]
    for p in ("allsat", "twosat"):
        S = pd.read_csv(O / "a7" / f"a7_series_{p}.csv")
        old = pd.read_csv(O / "a7" / f"a7_tests_{p}.csv")
        tests = [("primary", "EKE", "WAKE", "REST")] + [("secondary", v, "WAKE", "REST") for v in VARS[1:]]
        for f_ in FEAT:
            for v in ("EKE", "KE"):
                tests += [("secondary", v, f"down_{f_}", f"up_{f_}"), ("secondary", v, f"down_{f_}", "REST")]
        for per in ("2001-2025", "2001-2017"):
            sel = np.ones(len(S), bool) if per == "2001-2025" else S.year.values <= 2017
            for kind, v, a, b in tests:
                if per == "2001-2017" and a != "WAKE":
                    continue
                ya, yb = S[f"{a}_{v}"].values[sel], S[f"{b}_{v}"].values[sel]
                d = ya / ya.mean() - yb / yb.mean()
                r = old[(old.period == per) & (old["var"] == v) & (old.region_a == a) & (old.region_b == b)]
                run(f"{p}|{per}|{v}|{a}-{b}", S.x.values[sel], d, month=S.month.values[sel], scale=1000,
                    old_p=float(r.p_cal.iloc[0]) if len(r) else np.nan, old_test="p_cal",
                    product=p, period=per, kind=kind, var=v, region_a=a, region_b=b)


def fam_CONTRAST():
    run = Runner("CONTRAST")
    old = pd.read_csv(O / "stats" / "contrast_tests.csv").set_index("test")
    for p in ("allsat", "twosat"):
        S = pd.read_csv(O / "a5" / f"a5_series_{p}.csv")
        for reg in ("ACC", "Pacific", "Atlantic", "Indian"):
            d = S[f"{reg}_EKE"] / S[f"{reg}_EKE"].mean() - S[f"{reg}_KE"] / S[f"{reg}_KE"].mean()
            k = f"ACC EKE-KE {p}"
            run(f"{reg} EKE-KE {p}", S.x.values, d.values, month=S.month.values, scale=1000,
                old_p=float(old.loc[k, "p_cal"]) if (reg == "ACC" and k in old.index) else np.nan, old_test="p_cal")
        # sector ranking: Pacific minus each other sector (relative trends), EKE and KE
        for v in ("EKE", "KE"):
            for other in ("Atlantic", "Indian"):
                d = S[f"Pacific_{v}"] / S[f"Pacific_{v}"].mean() - S[f"{other}_{v}"] / S[f"{other}_{v}"].mean()
                run(f"Pacific-{other} {v} {p}", S.x.values, d.values, month=S.month.values, scale=1000)
    BX = pd.read_csv(O / "a6" / "a6b_drake_box_series.csv"); BX["x"] = BX.year + (BX.month - 0.5) / 12
    for mk, k in (("ice_any", "Drake box EKE-KE ice mask"), ("nomask", "Drake box EKE-KE no mask")):
        d = BX[f"{mk}_EKE"] / BX[f"{mk}_EKE"].mean() - BX[f"{mk}_KE"] / BX[f"{mk}_KE"].mean()
        run(f"Drake EKE-KE {mk} allsat", BX.x.values, d.values, month=BX.month.values, scale=1000,
            old_p=float(old.loc[k, "p_cal"]), old_test="p_cal")


def fam_B1():
    run = Runner("B1")
    P = pd.read_csv(O / "b1" / "b1_pilot_series.csv")
    old = pd.read_csv(O / "b1" / "b1_pilot_trends.csv").set_index("series")
    cols = [c for c in P.columns if c not in ("year", "month", "x", "i") and np.isfinite(P[c]).mean() >= 0.5]
    cols = sorted(cols, key=lambda c: (c not in ("up_core_speed", "mea_core_speed", "down_eke", "core_s100_int_abs_kappa"), c))
    for c in cols:
        run(f"pilot|{c}", P.x.values, P[c].values, month=P.month.values,
            old_p=float(old.loc[c, "p_cal"]) if c in old.index else np.nan, old_test="p_cal", series=c)
    import xarray as xr
    SECS = {"PAR": {"upstream": (-150, -135.13), "meander": (-135.13, -109.38), "downstream": (-109.38, -80)},
            "SEIR150": {"upstream": (125, 138), "meander": (140, 150), "downstream": (150, 158)},
            "KERG": {"upstream": (70, 77), "meander": (76, 84), "downstream": (84, 94)}}
    old = pd.read_csv(O / "b1" / "b1_paths_stats.csv")
    yrs = np.arange(2001, 2026)
    for n, secs in SECS.items():
        Pp = xr.open_dataset(O / "b1" / f"paths_{n}.nc")
        ann = pd.DataFrame(Pp.path_lat.values).groupby(np.repeat(yrs, 12)).mean().values
        for s, (a, b) in secs.items():
            m = (Pp.lon.values >= a) & (Pp.lon.values <= b)
            for var, Y in (("path_lat_deg", np.nanmean(ann[:, m], 1)), ("within_year_sd_deg", np.nanmean(Pp.annual_sd.values[:, m], 1))):
                r = old[(old.meander == n) & (old.section == s) & (old["var"] == var)]
                run(f"paths|{n}|{s}|{var}", yrs + 0.5, Y, monthly=False, scale=10, old_p=float(r.p.iloc[0]), old_test="sieve",
                    meander=n, section=s, var=var)


def fam_LIU():
    run = Runner("LIU")
    for p in ("allsat", "twosat"):
        M = pd.read_csv(O / "liu" / f"liu_monthly_metrics_{p}.csv")
        old = pd.read_csv(O / "liu" / f"liu_monthly_trends_{p}.csv").set_index("series")
        x = (M.year + (M.month - 0.5) / 12).values
        for c in [c for c in M.columns if c not in ("year", "month")]:
            run(f"monthly|{p}|{c}|2001-2025", x, M[c].values, month=M.month.values, scale=10,
                old_p=float(old.loc[c, "p_cal"]) if c in old.index else np.nan, old_test="p_cal", product=p, series=c)
    M = pd.read_csv(O / "liu" / "liu_metrics_g025.csv")
    old = pd.read_csv(O / "liu" / "liu_trends_g025.csv")
    x = (M.year + (M.month - 0.5) / 12).values
    head = [c for c in M.columns if c.endswith(("_position", "_width", "_speed", "_halfwidth")) or c.startswith("amplitude")]
    for c in head:
        for per in ("1993-2023(Liu)", "1993-2025"):
            sel = ((M.year < 2023) | ((M.year == 2023) & (M.month == 1))).values if per.startswith("1993-2023") else np.ones(len(M), bool)
            r = old[(old.series == c) & (old.period == per)]
            run(f"daily_g025|{c}|{per}", x[sel], M[c].values[sel], month=M.month.values[sel], scale=10,
                old_p=float(r.p_cal.iloc[0]) if len(r) else np.nan, old_test="p_cal", series=c, period=per)


if __name__ == "__main__":
    for fam in sys.argv[1:]:
        globals()[f"fam_{fam}"]()
    print("done", sys.argv[1:])
