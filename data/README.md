# Data files

| File | Rows | Produced by | Description |
|---|---|---|---|
| tarim_annual_ndvi_lst_2000_2023_series.csv | 24 | GEE export 1 | Growing-season corridor statistics from MOD13Q1 and MOD11A2, 2000-2023 |
| tarim_points_corridor_2000_2023.csv | 72,000 | GEE export 2 | 3,000 random points per year inside the 1 km corridor buffer, seed 42 |
| tarim_points_region_2000_2023.csv | 72,000 | GEE export 3 | 3,000 random points per year in the surrounding region, used to fit the TVDI edges |
| tarim_annual_final_2000_2023.csv | 24 | recalibrate_v9.py | Assembled annual table: W, R, NDVI, P, LST, TVDI |
| calibration_v10.json | - | anchor_v10.py | Published parameters and every fitted statistic quoted in the paper |
| table_observable_thresholds.csv | 9 | observable_thresholds.py | Table S4: thresholds in delivered volume and critical NDVI |
| bimodality_test.json | - | bimodality_test.py | Appendix D.5 mixture screen |

The three exports come from gee/tarim_gee_v6_run.js, run at
<https://code.earthengine.google.com>. The corridor is the 1 km buffer on each
side of the mapped Tarim channel between the Daxihaizi Reservoir
(87.3795 E, 40.6474 N) and Taitema Lake (88.3685 E, 39.4709 N).

The point tables are the original exports and are not modified by any script;
everything derived from them is written to outputs/.
