# Legacy scripts (audit only)

These three single-purpose scripts predate the consolidated pipeline and are
kept only so that the history of the numbers is traceable. **They are not part
of the reproduction path.** Paths were rewritten to repository-relative form
when the archive was assembled, so they run from a clean clone, but they read
the final exports in data/ rather than the intermediate files of the v7/v8
runs; treat their printed numbers as historical, not as the published values.

| Script | Superseded by |
|---|---|
| tvdi_v7.py | recalibrate_v9.py (edge fitting and TVDI) |
| trend_stats.py | recalibrate_v9.py (trend estimation) |
| build_tables_v8.py | recalibrate_v9.py (annual table assembly) |

Run order, if you do want to replay them:

```bash
python legacy/tvdi_v7.py       # outputs/tarim_annual_ndvi_tvdi_v7.csv
python legacy/trend_stats.py   # trend table printed to stdout
python legacy/build_tables_v8.py
```
