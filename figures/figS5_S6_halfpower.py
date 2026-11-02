"""SI Figures S5 and S6: sensitivity of the half-power (energetic-pathway) analysis. (Former si_figures.py modes S2 and S3;
the obsolete modes S1 mean-KE map, S4 tracked eddies and S5 path-length covariation are not carried over.)
S5  a circumpolar half-power latitude trend 2001-2025 for the baseline and each alternative detection parameter (km per decade,
      positive northward; retest_summary <set>|ACC|lat_hp); b trends in 5 deg longitude bins (retest lonbin|*);
    c percentage of synthetic realisations in which an injected linear shift is detected at p < 0.05 (a3_summary, S1_shift*km;
      S0_null = 0 km); d median recovered shift (dashed 1:1).
S6  half-power latitude (a) and KE at the envelope peak (b) trends 2001-2025 in Products 1 and 4, circumpolar and by sector
    (a2h_twosat_vs_allsat.csv).

Publication repository: ported from figures_final/scripts/figS5_S6_halfpower.py (2 Nov 2026); only paths changed (osr11.paths).
Inputs  : outputs/retest/retest_summary.csv, outputs/a3/a3_summary.csv, outputs/a2/a2h_twosat_vs_allsat.csv
Outputs : rendered/{png,pdf}/FigureS5, FigureS6
Run     : python figures/figS5_S6_halfpower.py [figure dir, default rendered/]
"""
import sys as _sys; from pathlib import Path as _Path  # noqa: E401,E702
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "src"))  # repository src/ (package osr11)
import sys
from pathlib import Path
import numpy as np, pandas as pd
from osr11 import style as S, paths as PA  # noqa: E402
FIGD = Path(sys.argv[1]) if len(sys.argv) > 1 else PA.FIG; OUT = PA.OUT
plt = S.plt; C = S.C; KM = 111.195
A_, B_ = "Product 1 (all-satellite)", "Product 4 (two-satellite)"
COLS = {A_: C["ink"], B_: C["two"]}; MK = {A_: "o", B_: "s"}
R = pd.read_csv(OUT / "retest" / "retest_summary.csv").set_index("key")


def kmrow(key):
    r = R.loc[key]; return dict(pct=r.slope * KM * 10, pct_lo=r.ci_lo * KM * 10, pct_hi=r.ci_hi * KM * 10, p=r.p_new)


# ---------------------------------------------------------------- S5
H = 100.0
fig = S.figure(S.W2, H)
sets = [("baseline", "Baseline (200 cm$^2$ s$^{-2}$, 25 %, 100 km)"), ("thr150", "Threshold 150 cm$^2$ s$^{-2}$"),
        ("thr250", "Threshold 250 cm$^2$ s$^{-2}$"), ("thr300", "Threshold 300 cm$^2$ s$^{-2}$"),
        ("thr400", "Threshold 400 cm$^2$ s$^{-2}$"), ("relthr", "Relative threshold"),
        ("frac20", "Peak fraction 20 %"), ("frac33", "Peak fraction 33 %"), ("frac50", "Peak fraction 50 %"),
        ("minw50", "Minimum width 50 km"), ("minw200", "Minimum width 200 km"), ("minw300", "Minimum width 300 km")]
P = "circumpolar position"
ax = S.axmm(fig, 46, 12, 34, 80)
S.forest(ax, [(lab, {P: kmrow(f"{s}|ACC|lat_hp")}) for s, lab in sets], {P: C["pos"]}, xlim=(-10, 12),
         xlabel="Trend (km per decade,\npositive northward)", gap=0.6)
ax.axhline(-0.5, color=C["light"], lw=0.5)
S.letter(fig, 2, H - 5, "a", "Detection parameters")
ax = S.axmm(fig, 98, 60, 68, 32)
L = sorted([(int(k.split("|")[1]) + 2.5, kmrow(k)) for k in R.index if k.startswith("lonbin|")])
lo_ = np.array([l for l, _ in L]); d = pd.DataFrame([r for _, r in L])
S.zero(ax); ax.vlines(lo_, d.pct_lo, d.pct_hi, color=C["pos"], lw=0.6)
S.profile(ax, lo_, d.pct, d.p, C["pos"], ms=2.6, lw=0)
S.lon_axis(ax, xlabel="Longitude (bin centre)"); ax.set_ylabel("km per decade"); S.yticks_minus(ax)
S.letter(fig, 86, H - 5, "b", "Trends in 5° longitude bins")
A3 = pd.read_csv(OUT / "a3" / "a3_summary.csv").set_index("scenario")
sh = [("S0_null", 0), ("S1_shift5km", 5), ("S1_shift10km", 10), ("S1_shift25km", 25), ("S1_shift50km", 50)]
xs = np.array([s for _, s in sh]); pw = np.array([A3.loc[k, "lat_hp_slope_pow"] for k, _ in sh]) * 100
rec = np.array([A3.loc[k, "lat_hp_slope_med"] for k, _ in sh]) * KM * 25
ax = S.axmm(fig, 98, 12, 28, 30)
ax.plot(xs, pw, "o-", color=C["pos"], ms=2.8, lw=0.8); ax.axhline(5, color=C["grey"], lw=0.5, ls="--")
ax.set_xlabel("Injected shift (km)"); ax.set_ylabel("Detected (%)"); ax.set_ylim(0, 105)
S.letter(fig, 86, 47, "c", "Detection")
ax = S.axmm(fig, 138, 12, 28, 30)
ax.plot(xs, rec, "o-", color=C["pos"], ms=2.8, lw=0.8); ax.plot([0, 50], [0, 50], color=C["grey"], lw=0.5, ls="--")
ax.set_xlabel("Injected shift (km)"); ax.set_ylabel("Recovered (km)")
S.letter(fig, 129, 47, "d", "Recovery")
S.save(fig, FIGD, "FigureS5"); plt.close(fig)
print("S5 lonbins", len(d), "sig", int((d.p.round(3) < 0.05).sum()), "power %", dict(zip(xs.tolist(), pw.round(1).tolist())),
      "recovered km", dict(zip(xs.tolist(), rec.round(1).tolist())), "baseline", {k: round(v, 2) for k, v in kmrow("baseline|ACC|lat_hp").items()})

# ---------------------------------------------------------------- S6
Hh = pd.read_csv(OUT / "a2" / "a2h_twosat_vs_allsat.csv")
H = 50.0
fig = S.figure(S.W2, H)
for x0, (fld, lim, xl, t, l) in zip((22, 108), (("pos", (-12, 16), "Trend (km per decade, positive northward)", "Half-power latitude", "a"),
                                               ("peak", (-4, 12), "Trend (% per decade)", "KE at the envelope peak", "b"))):
    rows = []
    for reg, lab in (("ACC", "Circumpolar"), ("Pacific", "Pacific"), ("Atlantic", "Atlantic"), ("Indian", "Indian")):
        dd = {}
        for prod, nm in (("all-sat", A_), ("two-sat", B_)):
            r = Hh[(Hh["product"] == prod) & (Hh.region == reg)].iloc[0]
            dd[nm] = (dict(pct=r.pos_km_dec, pct_lo=r.pos_lo, pct_hi=r.pos_hi, p=r.pos_p) if fld == "pos" else
                      dict(pct=r.peak_pct_dec, pct_lo=r.peak_lo, pct_hi=r.peak_hi, p=r.peak_p))
        rows.append((lab, dd))
    ax = S.axmm(fig, x0, 17, 58, 26)
    S.forest(ax, rows, COLS, MK, xlim=lim, xlabel=xl)
    S.letter(fig, 2 if l == "a" else 88, H - 4, l, t)
fig.legend(handles=S.sig_handles(COLS, MK), loc="center", ncol=3, bbox_to_anchor=(0.5, 3.0 / H))
S.save(fig, FIGD, "FigureS6"); plt.close(fig)
