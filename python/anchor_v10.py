"""Anchor the parameters at the observed window (Table 1 of the paper).

recalibrate_v9.py re-derives the record from the Earth Engine exports using
the *first* parameter rule (K = P_max + 0.05, R_rec = logistic inflection).
That rule is kept for audit only.  The published parameters use the rule
documented in Section 2.4:

    K    = P_max + 0.01                     saturation just above the record
    s0   = mu * TVDI_mean                   degraded-state stress = observed mean
    T*   = s0 / mu                          bare-branch stress
    T**  = stress balance at P = K          healthy-branch stress
    A2   = [K - A1 (T** - T*)] / dR         dR = 19 mm/yr, the observed width
    A0   = A2 R_rec - A1 T*                 R_rec = 55 mm/yr (2012 crossing)

With dR and R_rec fixed, R_deg = R_rec - dR = 36 mm/yr.  This script applies
the rule to outputs/calibration_v9.json and checks the result against
data/calibration_v10.json.

Usage:  python python/anchor_v10.py
"""
import json
import sys

import paths
from model_v8 import T_of_P, load_parameters, thresholds

R_REC = 55.0      # mm/yr, the input first reached in 2012
WIDTH = 19.0      # mm/yr, R_rec - R_deg


def anchor(rec):
    p = load_parameters()
    p["K"] = round(float(rec["P_max"]) + 0.01, 2)
    p["s0"] = p["mu"] * float(rec["TVDI_mean"])
    t_bare = p["s0"] / p["mu"]
    t_healthy = T_of_P(p["K"], p)[0]
    p["A2"] = (p["K"] - p["A1"] * (t_healthy - t_bare)) / WIDTH
    p["A0"] = p["A2"] * R_REC - p["A1"] * t_bare
    return p, t_bare, t_healthy


def main():
    with open(paths.OUT / "calibration_v9.json", encoding="utf-8") as fh:
        rec = json.load(fh)
    p, t_bare, t_healthy = anchor(rec)
    r_deg, r_rec, _, _ = thresholds(p)

    print("anchored parameters")
    for k, v in p.items():
        print("  %-4s %.8g" % (k, v))
    print("T* = %.4f   T** = %.4f" % (t_bare, t_healthy))
    print("R_deg = %.2f   R_rec = %.2f   width = %.2f"
          % (r_deg, r_rec, r_rec - r_deg))

    with open(paths.CALIBRATION, encoding="utf-8") as fh:
        pub = {k: float(v) for k, v in json.load(fh)["parameters"].items()}
    worst = max(abs(p[k] - pub[k]) / max(abs(pub[k]), 1e-12) for k in pub)
    print("largest relative difference from data/calibration_v10.json: %.2e"
          % worst)

    dest = paths.OUT / "parameters_anchored.json"
    with open(dest, "w", encoding="utf-8") as fh:
        json.dump(p, fh, indent=2)
    print("written", dest)
    if worst > 1e-6:
        sys.exit("anchored parameters do not reproduce the published set")


if __name__ == "__main__":
    main()
