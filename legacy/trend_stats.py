"""Trend statistics of the v7 annual series (NDVI, TVDI, LST)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ROOT = Path(__file__).resolve().parents[1]
ann = pd.read_csv(ROOT / "data" / "tarim_annual_ndvi_lst_2000_2023_series.csv")
tv = pd.read_csv(ROOT / "outputs" / "tarim_annual_ndvi_tvdi_v7.csv")
m = ann.merge(tv[["year", "TVDI"]], on="year")


def trend(x, y):
    s, i = np.polyfit(x, y, 1)
    pred = s * x + i
    n = len(x)
    res = y - pred
    sse = float(np.sum(res ** 2))
    sxx = float(np.sum((x - x.mean()) ** 2))
    se = np.sqrt(sse / (n - 2) / sxx)
    r2 = 1 - sse / float(np.sum((y - y.mean()) ** 2))
    return s, se, s / se, r2


for name in ("NDVI_mean", "LST_C", "TVDI"):
    s, se, t, r2 = trend(m.year.to_numpy(), m[name].to_numpy())
    print("%-10s slope %+.5f /yr  se %.5f  t %+.2f  R2 %.2f" % (name, s, se, t, r2))

print()
print("NDVI 2000 %.4f  2010 %.4f  2020 %.4f  2023 %.4f" % tuple(
    m.set_index("year").loc[[2000, 2010, 2020, 2023], "NDVI_mean"]))
print("TVDI 2000 %.3f  2010 %.3f  2020 %.3f  2023 %.3f" % tuple(
    m.set_index("year").loc[[2000, 2010, 2020, 2023], "TVDI"]))
print("LST  2000 %.2f  2010 %.2f  2020 %.2f  2023 %.2f" % tuple(
    m.set_index("year").loc[[2000, 2010, 2020, 2023], "LST_C"]))
print("TVDI mean %.3f sd %.3f min %.3f max %.3f" % (m.TVDI.mean(), m.TVDI.std(), m.TVDI.min(), m.TVDI.max()))
nd = (m.NDVI_mean.to_numpy() - 0.05) / 0.15
print("P    2000 %.3f  2010 %.3f  2020 %.3f  2023 %.3f  max %.3f" %
      (nd[0], nd[10], nd[20], nd[23], nd.max()))
