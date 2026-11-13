"""Assemble the final annual table (2000-2023) and the calibration report."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
import model_v8 as M

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
DATA = ROOT / "data"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
W = {2000: 1.02, 2001: 2.53, 2002: 3.21, 2003: 3.45, 2004: 3.82, 2005: 4.15,
     2006: 4.23, 2007: 4.36, 2008: 4.52, 2009: 4.78, 2010: 5.03, 2011: 5.21,
     2012: 5.47, 2013: 5.86, 2014: 6.12, 2015: 6.35, 2016: 6.58, 2017: 6.82,
     2018: 7.05, 2019: 7.28, 2020: 8.52, 2021: 7.35, 2022: 7.18, 2023: 7.42}
NDVI_MIN, NDVI_MAX = 0.050, 0.200

ann = pd.read_csv(DATA / "tarim_annual_ndvi_lst_2000_2023_series.csv")
tv = pd.read_csv(OUT / "tarim_annual_ndvi_tvdi_v7.csv")
ann["year"] = ann.year.astype(int)
df = ann[["year", "NDVI_mean", "NDVI_p50", "NDVI_p90", "LST_C"]].merge(
    tv[["year", "TVDI"]], on="year")
df["W"] = df.year.map(W)
df["R"] = 28.0 + 5.0 * df.W
df["P"] = np.clip((df.NDVI_mean - NDVI_MIN) / (NDVI_MAX - NDVI_MIN), 0, 1)
df = df.rename(columns={"NDVI_mean": "NDVI", "LST_C": "LST"})
df = df[["year", "W", "R", "NDVI", "NDVI_p50", "NDVI_p90", "P", "LST", "TVDI"]]

p = M.BASE
Rg, Rr, Tb, Th = M.thresholds(p)

print(df.to_string(index=False, float_format=lambda v: "%.3f" % v))
print("\nW total %.2f (1e8 m3)" % df.W.sum())
print("R %.2f - %.2f mm/yr" % (df.R.min(), df.R.max()))
print("P %.3f - %.3f" % (df.P.min(), df.P.max()))
print("TVDI %.3f +- %.3f" % (df.TVDI.mean(), df.TVDI.std()))
print("\nthresholds R_deg %.2f  R_rec %.2f  window %.2f" % (Rg, Rr, Rr - Rg))
print("model stress: bare %.4f  healthy %.4f" % (Tb, Th))

print("\nmodel stress at the observed cover (healthy branch):")
for y in (2000, 2010, 2020, 2023):
    r = df[df.year == y].iloc[0]
    T = M.T_of_P(float(r.P), p)
    print("  %d  P=%.3f  T_model=%s  TVDI=%.3f" % (y, r.P, ["%.4f" % t for t in T], r.TVDI))

print("\nsensitivity of the thresholds to +-10% (mm/yr):")
names = ["rho", "A0", "A1", "A2", "K", "s0", "a", "g", "eta", "mu"]
print("%-6s %-14s %-14s" % ("par", "-20%", "+20%"))
rows = []
for nm in names:
    out = []
    for f in (0.9, 1.1):
        q = dict(p)
        q[nm] = p[nm] * f
        try:
            g_, r_ = M.thresholds(q)[:2]
        except IndexError:
            g_, r_ = float("nan"), float("nan")
        out.append((g_, r_))
    rows.append((nm, out[0], out[1]))
    print("%-6s %6.1f / %6.1f  %6.1f / %6.1f" % (nm, out[0][0], out[0][1], out[1][0], out[1][1]))

print("\nequilibrium states:")
for R in (20.0, 28.0, 30.0, 33.1, 41.0, 52.0, 60.0, 65.1, 70.0):
    eq = M.equilibria(R, p)
    print("  R=%5.1f  " % R + " | ".join(
        "P=%.3f T=%.3f %s" % (e["P"], e["T"], "S" if e["stable"] else "U") for e in eq))

out_csv = OUT / "tarim_annual_final_2000_2023.csv"
df.to_csv(out_csv, index=False, float_format="%.4f")
print("\nwritten:", out_csv)
