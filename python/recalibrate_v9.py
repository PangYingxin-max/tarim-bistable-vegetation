"""One-command recalibration from a fresh Earth Engine export (v6 channel).

    python python/recalibrate_v9.py --annual A.csv --corridor C.csv --region R.csv

Without arguments the files in data/ are used.  The script

  1. fits the dry / wet edges of the NDVI-LST feature space on the regional
     point table (pooled slopes, per-year intercepts, as in Section 2.1),
  2. evaluates the TVDI of the corridor points,
  3. converts corridor NDVI to fractional cover with the pooled 5th / 98th
     percentiles as endmembers,
  4. estimates the linear trends with their standard errors and p-values,
  5. re-derives the model parameters with the rule documented in Table S2,
  6. verifies that the new parameter set still folds, and
  7. writes outputs/calibration_v9.json with every number quoted in the text.
"""
import argparse
import glob
import io
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model_v8 as M
import paths

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# annual ecological water transfer, 10^8 m^3 (Tarim River Basin Authority)
W = {2000: 1.02, 2001: 2.53, 2002: 3.21, 2003: 3.45, 2004: 3.82, 2005: 4.15,
     2006: 4.23, 2007: 4.36, 2008: 4.52, 2009: 4.78, 2010: 5.03, 2011: 5.21,
     2012: 5.47, 2013: 5.86, 2014: 6.12, 2015: 6.35, 2016: 6.58, 2017: 6.82,
     2018: 7.05, 2019: 7.28, 2020: 8.52, 2021: 7.35, 2022: 7.18, 2023: 7.42}
R0, KAPPA = 28.0, 5.0


def newest(pattern):
    hits = sorted(glob.glob(os.path.join(str(paths.DATA), pattern)),
                  key=os.path.getmtime)
    return hits[-1] if hits else None


def sig2(x):
    """Round to two significant figures (0.0498 -> 0.05, 0.197 -> 0.2)."""
    return float("%.2g" % x)


def binned_edge(g, q, binw=0.02, min_n=60, ndvi_max=0.50):
    """Quantile envelope of LST in NDVI bins, fitted by OLS (tvdi_v6 method)."""
    nd = g.NDVI.to_numpy()
    ls = g.LST.to_numpy()
    keep = nd <= ndvi_max
    nd, ls = nd[keep], ls[keep]
    lo = np.floor(nd.min() / binw) * binw
    lab = np.digitize(nd, np.arange(lo, nd.max() + binw, binw))
    xs, ys = [], []
    for k in np.unique(lab):
        m = lab == k
        if m.sum() < min_n:
            continue
        xs.append(nd[m].mean())
        ys.append(np.quantile(ls[m], q))
    xs, ys = np.array(xs), np.array(ys)
    s, i = np.polyfit(xs, ys, 1)
    return float(s), float(i), len(xs)


def edges(reg, q=0.05):
    """Pooled dry / wet edge slopes, as in Section 2.1 and Appendix C.2."""
    sw, iw, nw = binned_edge(reg, q)
    sd, idr, nd = binned_edge(reg, 1 - q)
    return sd, sw, min(nw, nd)


def tvdi(corr, reg, dry_s, wet_s):
    out = []
    for y in sorted(corr.year.unique()):
        r = reg[reg.year == y]
        c = corr[corr.year == y]
        rw = r.LST.to_numpy() - wet_s * r.NDVI.to_numpy()
        rd = r.LST.to_numpy() - dry_s * r.NDVI.to_numpy()
        wi = float(np.quantile(rw, 0.05))
        di = float(np.quantile(rd, 0.95))
        nd, ls = c.NDVI.to_numpy(), c.LST.to_numpy()
        den = (dry_s - wet_s) * nd + (di - wi)
        t = np.clip((ls - (wet_s * nd + wi)) / den, 0.0, 1.0)
        out.append(dict(year=int(y), TVDI=float(t.mean()),
                        TVDI_sd=float(t.std()),
                        frac_at_one=float((t >= 0.999).mean())))
    return pd.DataFrame(out)


def trend(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    s, i = np.polyfit(x, y, 1)
    res = y - (s * x + i)
    n = len(x)
    sxx = float(np.sum((x - x.mean()) ** 2))
    se = float(np.sqrt(np.sum(res ** 2) / (n - 2) / sxx))
    r2 = 1.0 - float(np.sum(res ** 2)) / float(np.sum((y - y.mean()) ** 2))
    t = s / se
    try:
        from scipy import stats
        p = float(2 * stats.t.sf(abs(t), n - 2))
    except Exception:                                    # noqa: BLE001
        p = float("nan")
    return dict(slope=s, se=se, t=t, r2=r2, p=p)


def calibrate(df, tvdi_mean, pmax):
    """Parameter rule of Table S2, re-solved on the new record.

    (i)   sigma0 = mu * TVDI_mean      -> modelled degraded-state stress
    (ii)  K      = pmax + 0.05         -> saturation cover
    (iii) A1 and A2 keep their assumed values
    (iv)  A0 solves R_rec = R at the inflection year of a logistic fit to P
    (v)   the degradation threshold follows; bistability is then checked
    """
    p = dict(M.BASE)
    p["s0"] = p["mu"] * tvdi_mean
    p["K"] = round(float(pmax) + 0.05, 2)          # saturation cover, 2 d.p.

    y = df.year.to_numpy(float)
    P = df.P.to_numpy(float)
    # logistic fit by grid search (no scipy dependency for the fit itself)
    best = None
    for r in np.linspace(0.2, 3.0, 57):
        for t0 in np.linspace(2002, 2021, 77):
            f = 1.0 / (1.0 + np.exp(-r * (y - t0)))
            lo = f.min()
            hi = f.max()
            g = P.min() + (P.max() - P.min()) * (f - lo) / (hi - lo)
            sse = float(np.sum((P - g) ** 2))
            if best is None or sse < best[0]:
                best = (sse, r, t0)
    _, rate, t_infl = best
    year_infl = int(round(t_infl))
    R_infl = float(df.loc[df.year == year_infl, "R"].iloc[0])

    Tb = p["s0"] / p["mu"]
    Th = M.T_of_P(p["K"], p)[0]
    R_target = float(round(R_infl))                # recovery target, 1 mm/yr
    p["A0"] = p["A2"] * R_target - p["A1"] * Tb
    Rg, Rr, Tb, Th = M.thresholds(p)
    return p, dict(logistic_rate=rate, inflection_year=year_infl,
                   R_at_inflection=R_infl, R_target=R_target,
                   R_deg=Rg, R_rec=Rr, width=Rr - Rg,
                   T_bare=Tb, T_healthy=Th,
                   fold_gap=p["K"] - p["A1"] * (Th - Tb))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annual")
    ap.add_argument("--corridor")
    ap.add_argument("--region")
    ap.add_argument("--out", default=str(paths.OUT / "calibration_v9.json"))
    a = ap.parse_args()

    f_ann = a.annual or newest("tarim_annual_ndvi_lst_2000_2023_series*.csv")
    f_cor = a.corridor or newest("tarim_points_corridor_2000_2023*.csv")
    f_reg = a.region or newest("tarim_points_region_2000_2023*.csv")
    print("annual  :", f_ann)
    print("corridor:", f_cor)
    print("region  :", f_reg)

    ann = pd.read_csv(f_ann)
    ann["year"] = ann.year.astype(int)
    if not f_cor or not f_reg:
        # pre-computed mode: the annual table already carries NDVI and TVDI
        print("\n(pre-computed mode: NDVI and TVDI read from the annual table)")
        ren = {"NDVI": "NDVI_mean", "LST": "LST_C"}
        df = ann.rename(columns=ren)
        lo = float(df.P.min()) if "P" in df else float("nan")
        hi = float(df.P.max()) if "P" in df else float("nan")
        if "P" not in df:
            lo_n, hi_n = 0.050, 0.200
            df["P"] = np.clip((df.NDVI_mean - lo_n) / (hi_n - lo_n), 0, 1)
        res_edges = None
        df = df[["year", "NDVI_mean", "NDVI_p50", "NDVI_p90", "LST_C", "TVDI", "P"]]
    else:
        cor = pd.read_csv(f_cor).dropna(subset=["NDVI", "LST"])
        reg = pd.read_csv(f_reg).dropna(subset=["NDVI", "LST"])
        for d in (cor, reg):
            d["year"] = d.year.astype(int)
        dry_s, wet_s, nb = edges(reg)
        print("\nedge slopes: dry %.3f  wet %.3f  (%d NDVI bins)"
              % (dry_s, wet_s, nb))
        tv = tvdi(cor, reg, dry_s, wet_s)
        df = ann[["year", "NDVI_mean", "NDVI_p50", "NDVI_p90",
                  "LST_C"]].merge(tv, on="year")
        lo_raw, hi_raw = (float(v) for v in np.percentile(cor.NDVI, [5, 98]))
        lo, hi = sig2(lo_raw), sig2(hi_raw)      # endmembers as quoted, 2 s.f.
        print("NDVI endmembers: raw %.4f / %.4f  ->  %.3f / %.3f"
              % (lo_raw, hi_raw, lo, hi))
        df["P"] = np.clip((df.NDVI_mean - lo) / (hi - lo), 0, 1)
        res_edges = dict(dry_slope=dry_s, wet_slope=wet_s, bins=nb)
    df["W"] = df.year.map(W)
    df["R"] = R0 + KAPPA * df.W

    print(df.to_string(index=False, float_format=lambda v: "%.3f" % v))

    res = dict(files=dict(annual=f_ann, corridor=f_cor, region=f_reg),
               edges=res_edges,
               ndvi_min=(float(lo) if res_edges else None),
               ndvi_max=(float(hi) if res_edges else None),
               years=[int(v) for v in df.year],
               NDVI=[float(v) for v in df.NDVI_mean],
               LST=[float(v) for v in df.LST_C],
               P=[float(v) for v in df.P],
               TVDI=[float(v) for v in df.TVDI],
               R=[float(v) for v in df.R],
               W_total=float(df.W.sum()))

    for nm, col in (("NDVI", "NDVI_mean"), ("LST", "LST_C"),
                    ("TVDI", "TVDI"), ("P", "P")):
        res["trend_" + nm] = trend(df.year, df[col])

    res["TVDI_mean"] = float(df.TVDI.mean())
    res["TVDI_sd"] = float(df.TVDI.std())
    res["P_max"] = float(df.P.max())
    res["P_max_year"] = int(df.loc[df.P.idxmax(), "year"])

    p, cal = calibrate(df, res["TVDI_mean"], res["P_max"])
    res["parameters"] = {k: float(v) for k, v in p.items()}
    res["calibration"] = cal

    # independent check: does continuation see two transitions?
    found = M.scan_thresholds(p, lo=0.0, hi=max(120.0, 1.6 * df.R.max()))
    res["continuation"] = [list(f) for f in found]

    print("\ntrends")
    for nm in ("NDVI", "LST", "TVDI", "P"):
        t = res["trend_" + nm]
        print("  %-5s %+.5f /yr  se %.5f  t %+.2f  R2 %.2f  p %.3f"
              % (nm, t["slope"], t["se"], t["t"], t["r2"], t["p"]))
    print("\nTVDI mean %.3f +- %.3f   P max %.3f (%d)"
          % (res["TVDI_mean"], res["TVDI_sd"], res["P_max"], res["P_max_year"]))
    print("logistic inflection year %d (R = %.1f)"
          % (cal["inflection_year"], cal["R_at_inflection"]))
    print("thresholds: R_deg %.2f  R_rec %.2f  width %.2f  fold gap %.4f"
          % (cal["R_deg"], cal["R_rec"], cal["width"], cal["fold_gap"]))
    print("continuation:", res["continuation"])

    with io.open(a.out, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)
    print("\nwritten", a.out)

    # annual table in the layout used by make_figs_v8.py
    tbl = df.rename(columns={"NDVI_mean": "NDVI", "LST_C": "LST"})[
        ["year", "W", "R", "NDVI", "NDVI_p50", "NDVI_p90", "P", "LST", "TVDI"]]
    tbl_path = paths.OUT / "tarim_annual_final_2000_2023.csv"
    tbl.to_csv(tbl_path, index=False, float_format="%.4f")
    print("written", tbl_path)


if __name__ == "__main__":
    main()
