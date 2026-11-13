# Bistable vegetation-water model of the lower Tarim River

Code and data for

> Pang, Y. & Shen, J. (2026). *Modelling hysteresis and resilience in arid
> forest ecosystems: a bistable vegetation-water model of the Populus
> euphratica dieback-recovery cycle.* (manuscript)

The model is a local two-state ODE for a corridor-scaled cover index P and a
water-stress index T:

```
dP/dt = rho P (P - A) (K - P)
dT/dt = s0 - a P/(1 + g T) + eta P - mu T
A(T,R) = A0 + A1 T - A2 R
```

With the baseline parameters the equilibrium curve is S-shaped, with two
saddle-node bifurcations at R_deg = 36.0 and R_rec = 55.0 mm/yr.

## Data sources

No MODIS data are redistributed here: the repository holds the Earth Engine
script that extracts them and the tables it exported. Everything below is
public.

| Input | Product | Access |
|---|---|---|
| NDVI, 250 m, 16-day | MOD13Q1 v061 | <https://doi.org/10.5067/MODIS/MOD13Q1.061> · [GEE catalogue](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13Q1) |
| Land surface temperature, 1 km, 8-day | MOD11A2 v061 | <https://doi.org/10.5067/MODIS/MOD11A2.061> · [GEE catalogue](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD11A2) |
| River centreline | OpenStreetMap, waterway=river | <https://www.openstreetmap.org> (Overpass API; the 43 vertices used are printed in section 1 of the script) |
| Annual ecological water transfer | Tarim River Basin Authority | no machine-readable release; the series used is Table S1 of the paper and is carried in the W dictionary of recalibrate_v9.py |

The MODIS products are distributed by the NASA LP DAAC and are in the public
domain. The OpenStreetMap extract is © OpenStreetMap contributors, available
under the ODbL.

## Quick start

```bash
python -m pip install -r requirements.txt
python python/model_v8.py              # thresholds and equilibria
python python/make_figs_v8.py          # Figures 2-4  -> outputs/
python python/identifiability_v8.py    # Appendix D     -> outputs/
```

Regenerating the satellite numbers from scratch needs a free Google Earth
Engine account:

1. Open <https://code.earthengine.google.com> and paste gee/tarim_gee_v6_run.js.
2. Run the three exports one at a time (the free tier fails if all three are
   queued at once) and download the CSVs from Google Drive.
3. Then, with the three files in data/:

```bash
python python/recalibrate_v9.py   # edges, TVDI, trends -> outputs/calibration_v9.json
python python/anchor_v10.py       # published parameters, checked against data/calibration_v10.json
```

## Which script produces which number

| Result in the paper | Script | Output |
|---|---|---|
| Edge slopes, endmembers, TVDI, trends, Table 1 | recalibrate_v9.py | outputs/calibration_v9.json |
| Published parameters, Table 1 and Table S2 | anchor_v10.py | outputs/parameters_anchored.json |
| Thresholds R_deg, R_rec, width | model_v8.py | stdout |
| Appendix D, Monte-Carlo robustness, Table S5 | identifiability_v8.py | outputs/identifiability_v8.json |
| Table S4, critical NDVI | observable_thresholds.py | outputs/table_observable_thresholds.csv |
| Figures 2, 3, 4 | make_figs_v8.py | outputs/figure*.pdf |
| Figure 5, Appendix A.3 | pattern_csh_v9.py | outputs/figure5_patterns.* |
| Appendix D.5 screen | bimodality_test.py | outputs/bimodality_test.json |

## Repository layout

```
gee/     Earth Engine script (single source of the MODIS extraction)
python/  analysis and figure scripts
data/    exported and derived tables, plus the published calibration
figures/ the five figures as they appear in the manuscript
legacy/  superseded single-purpose scripts, kept for audit only
outputs/ written at run time (not tracked)
```

## Notes

- The pattern simulation in Figure 5 is a phenomenological complex
  Swift-Hohenberg model. It is not a spatial extension of the ODE and it is not
  fitted to the corridor.
- recalibrate_v9.py re-derives the parameters with the first published rule
  (K = P_max + 0.05, R_rec at the logistic inflection). It is retained as an
  audit trail; the parameters used in the paper are produced by anchor_v10.py,
  which anchors the window at 36 and 55 mm/yr. Both are reproducible from the
  same inputs, and anchor_v10.py checks itself against
  data/calibration_v10.json.
- The water-transfer series has no machine-readable official release. It was
  transcribed from the Authority's annual conveyance announcements and is
  fixed in Table S1 of the paper; recalibrate_v9.py carries the same values in
  the W dictionary.

## Licence

Code: MIT (LICENSE). Data tables and figures: CC BY 4.0 (LICENSE-DATA.md).

## What goes in the repository — and what does not

Everything here is what a reader needs to regenerate the numbers, and nothing
that only exists on the author's account:

| Included | Why |
|---|---|
| `gee/tarim_gee_v6_run.js` | the extraction script itself; it can be re-run by anyone with a free Earth Engine account |
| `data/*.csv` | the exports it produced, so the analysis can be checked without re-running GEE |
| `python/*.py` | every table, threshold and figure, each mapped to a result in the table above |
| `figures/*.pdf` | the artwork as submitted |
| `CITATION.cff`, `.zenodo.json` | machine-readable citation and archive metadata |

Not included: the Earth Engine *share link* of a script or asset. Such links
are tied to the author's account and can be withdrawn, so a reviewer may not be
able to open them. The deposited artefact is the script text plus the CSV
exports, which is what the README documents and what the Data Availability
statement should point to. The DOI of the archived release, not a Drive or
Earth Engine URL, is the citable locator.

Suggested Data Availability statement once the repository is public and
archived (replace the placeholders):

> The model code, the Earth Engine extraction script and all derived tables are
> available at https://github.com/PangYingxin-max/tarim-bistable-vegetation and
> archived at https://doi.org/10.5281/zenodo.24191944 (version 某个版本.0:
> https://doi.org/10.5281/zenodo.24191945). MODIS NDVI (MOD13Q1
> v061) and land-surface temperature (MOD11A2 v061) are available from the NASA
> LP DAAC (<https://doi.org/10.5067/MODIS/MOD13Q1.061>,
> <https://doi.org/10.5067/MODIS/MOD11A2.061>). The ecological water-transfer
> series is given in Table S1.

## Citation

See CITATION.cff. Please cite the paper, and the archived release of this
repository (DOI) when referring to the code or data.
