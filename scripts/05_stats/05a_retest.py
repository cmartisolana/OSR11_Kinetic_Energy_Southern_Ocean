"""Re-run every reported trend with the bias-corrected, calibrated test (rstats.trend_bcz).

Reads the saved series of each analysis (no field is recomputed) and writes one row per test to
outputs/retest/retest_<family>.csv, with the previously reported p for comparison. Rows already
present are skipped, so the script can be restarted.
       env RETEST_N_OUTER (default 300) = number of calibration nulls (p floor 1/(n_outer+1)).

Purpose : one calibrated test (rstats.trend_bcz: bias-corrected AR surrogates, 300 null series, larger p of the AIC-order and
          AR(12) noise models) for every reported trend, from the saved series. Families kept in the final pipeline:
            A2        half-power latitude (baseline, sensitivity sets, 5-deg bins)      needs outputs/a2/tracks_<set>.npz + $OSR11_CACHE/lon.npy
            A5        KE/EKE and components, ACC band, sectors, both products           outputs/a5/*
            A7        relative EKE trend downstream of topography vs rest (p = 0.95)    outputs/a7/*
            CONTRAST  EKE-KE relative-trend contrasts (ACC +3.0 points per decade, p = 0.060)  outputs/a5/*, outputs/a6/a6b_*
          The families A6 (withdrawn transport proxies), B1 (meander paths) and LIU (Liu et al. reproduction) are archived;
          their retest_*.csv files are in archive/ and no longer enter retest_summary.csv (05b).
          The 'old_p' column of CONTRAST rows is read from the superseded outputs/stats/contrast_tests.csv when it exists; that
          file is not distributed (it held an older test, p = 0.0099), so on a fresh run old_p is NaN. old_p is a reference
          column only and is never quoted.
Inputs  : as listed per family
Outputs : outputs/retest/retest_<family>.csv (rows already present are skipped, so the shipped files make this a no-op)
Run     : python scripts/05_stats/05a_retest.py A2 A5 A7 CONTRAST          (env RETEST_N_OUTER, default 300)
Seeds   : trend_bcz(seed=0).
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "src"))  # repository src/ (package osr11)
import json
import os
import sys
import numpy as np, pandas as pd

from osr11 import rstats as R, paths as PA  # noqa: E402

O = PA.OUT
RT = PA.out("retest")
N_OUTER = int(os.environ.get("RETEST_N_OUTER", 300))


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
    from osr11 import halfpower as A
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
    f_old = O / "stats" / "contrast_tests.csv"   # superseded older test; reference column only (see docstring)
    old = pd.read_csv(f_old).set_index("test") if f_old.exists() else pd.DataFrame(columns=["p_cal"])
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
            old_p=float(old.loc[k, "p_cal"]) if k in old.index else np.nan, old_test="p_cal")


if __name__ == "__main__":
    for fam in sys.argv[1:]:
        globals()[f"fam_{fam}"]()
    print("done", sys.argv[1:])
