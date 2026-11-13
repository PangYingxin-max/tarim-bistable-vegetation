"""Repository-relative paths, so every script runs from a clean clone.

Import this module instead of hard-coding absolute paths::

    from paths import DATA, FIG, OUT
    df = pd.read_csv(DATA / "tarim_annual_final_2000_2023.csv")

OUT is created on import; everything the scripts write goes there, so the
tracked inputs in data/ and figures/ are never overwritten.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "figures"
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

ANNUAL = DATA / "tarim_annual_final_2000_2023.csv"
SERIES = DATA / "tarim_annual_ndvi_lst_2000_2023_series.csv"
POINTS_CORRIDOR = DATA / "tarim_points_corridor_2000_2023.csv"
POINTS_REGION = DATA / "tarim_points_region_2000_2023.csv"
CALIBRATION = DATA / "calibration_v10.json"
