"""TVDI with pooled edge slopes and year-specific edge intercepts.

The pooled fit (tvdi_v6.py) shows that year-to-year differences in LST are
larger than the NDVI-driven part of the feature space, so a single pair of
edges makes TVDI saturate in the hot years.  Here the SLOPES are kept from the
pooled fit (they are the identifiable part of the geometry) and the
INTERCEPTS are re-estimated for every year as the 5th and 95th percentile of
the residuals of the regional sample.  TVDI then measures where the corridor
sits between the wet and the dry edge of its own year.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "outputs" / "tarim_annual_ndvi_tvdi_v7.csv"
OUT.parent.mkdir(exist_ok=True)

# v7 was run on the pre-final GEE exports (" ... (2).csv"); the archived copies
# in data/ carry the same columns, so a re-run may differ in the last decimal.
corr = pd.read_csv(DATA / "tarim_points_corridor_2000_2023.csv").dropna(subset=["NDVI", "LST"])
reg = pd.read_csv(DATA / "tarim_points_region_2000_2023.csv").dropna(subset=["NDVI", "LST"])

W_S, D_S = 4.034, -1.160          # pooled quantile-5% slopes, tvdi_v6.py

rows = []
for y in sorted(corr.year.unique()):
    r = reg[reg.year == y]
    c = corr[corr.year == y]
    rw = r.LST.to_numpy() - W_S * r.NDVI.to_numpy()
    rd = r.LST.to_numpy() - D_S * r.NDVI.to_numpy()
    wi = np.quantile(rw, 0.05)
    di = np.quantile(rd, 0.95)
    nd, ls = c.NDVI.to_numpy(), c.LST.to_numpy()
    den = (D_S - W_S) * nd + (di - wi)
    t = np.clip((ls - (W_S * nd + wi)) / den, 0, 1)
    rows.append(dict(year=int(y), NDVI=nd.mean(), LST_C=ls.mean(), TVDI=t.mean(),
                     TVDI_p10=np.percentile(t, 10), TVDI_p90=np.percentile(t, 90),
                     frac_at_1=float((t >= 0.999).mean()), wet_i=wi, dry_i=di))

s = pd.DataFrame(rows)
s.to_csv(OUT, index=False, float_format="%.4f")
print(s.to_string(index=False, float_format=lambda v: "%.3f" % v))
print("\nwritten:", OUT)
