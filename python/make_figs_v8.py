"""Figures 2-4 of the manuscript, from the 2000-2023 MODIS record.

Reads data/calibration_v10.json and data/tarim_annual_final_2000_2023.csv,
writes outputs/figure{2,3,4}_*.pdf.
"""
import json
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import model_v8 as M
import paths

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
OUT = paths.OUT

with open(paths.CALIBRATION, encoding="utf-8") as fh:
    p = {k: float(v) for k, v in json.load(fh)["parameters"].items()}
obs = pd.read_csv(paths.ANNUAL)
print("figures use", paths.CALIBRATION)
Rg, Rr, Tb, Th = M.thresholds(p)
print("thresholds: R_deg %.2f  R_rec %.2f" % (Rg, Rr))

# ---------------- Figure 2: bifurcation diagram ----------------
Rs = np.linspace(0, 80, 801)
low = [(R, 0.0) for R in Rs]
up = [(R, p["K"]) for R in Rs]
mid = [e for R in Rs for e in M.equilibria(R, p) if 1e-6 < e["P"] < p["K"] - 1e-6]

fig, ax = plt.subplots(figsize=(6.6, 4.8))
ax.axvspan(Rg, Rr, color="0.90", zorder=0)
ax.text(0.5 * (Rg + Rr), 0.455, "bistable window", ha="center", fontsize=9, color="0.35")
ax.plot([v[0] for v in low], [v[1] for v in low], color="tab:blue", lw=2.2)
ax.plot([v[0] for v in up], [v[1] for v in up], color="tab:blue", lw=2.2)
# middle branch: R from R_deg to R_rec
Rmid, Pmid = [], []
for R in np.linspace(Rg + 1e-3, Rr - 1e-3, 400):
    for e in M.equilibria(R, p):
        if 1e-6 < e["P"] < p["K"] - 1e-6:
            Rmid.append(R); Pmid.append(e["P"])
ax.plot(Rmid, Pmid, "--", color="tab:red", lw=2, label="unstable (saddle)")

ax.plot(obs.R, obs.P, "o-", color="tab:green", ms=3.6, lw=1.1,
        mfc="white", label="observed 2000-2023")
for y in (2000, 2013, 2023):
    r = obs[obs.year == y].iloc[0]
    ax.annotate(str(y), (r.R, r.P), textcoords="offset points", xytext=(4, 5), fontsize=7.5,
                color="tab:green")
ax.axvline(Rg, color="0.35", ls=":", lw=1)
ax.axvline(Rr, color="0.35", ls=":", lw=1)
ax.text(Rg, -0.052, r"$R_{\rm deg}$=%.0f" % Rg, ha="center", fontsize=8.5)
ax.text(Rr, -0.052, r"$R_{\rm rec}$=%.0f" % Rr, ha="center", fontsize=8.5)
ax.annotate("", xy=(Rr, 0.30), xytext=(Rg, 0.30),
            arrowprops=dict(arrowstyle="->", color="0.45", lw=1.1))
ax.text(0.5 * (Rg + Rr), 0.312, "recovery", ha="center", fontsize=8, color="0.45")
ax.annotate("", xy=(Rg, 0.045), xytext=(Rr, 0.045),
            arrowprops=dict(arrowstyle="->", color="0.45", lw=1.1))
ax.text(0.5 * (Rg + Rr), 0.057, "collapse", ha="center", fontsize=8, color="0.45")
ax.set_xlabel(r"Effective water input  $R$  (mm yr$^{-1}$)")
ax.set_ylabel(r"Vegetation cover $P^{*}$ (equilibrium) and observed $P$")
ax.set_xlim(0, 80); ax.set_ylim(-0.075, 0.50)
ax.legend(loc="upper left", frameon=False, fontsize=8.5)
fig.tight_layout()
fig.savefig(OUT / "figure2_bistability.pdf")
fig.savefig(OUT / "figure2_bistability.png", dpi=150)
plt.close(fig)

# ---------------- Figure 3: phase portrait ----------------
Rc = float(round(0.5 * (Rg + Rr)))      # middle of the bistable window
print("phase portrait at R = %.1f (quote this value in Section 3.2)" % Rc)
pcol = np.linspace(-0.01, 0.45, 24)
tcol = np.linspace(0.84, 0.95, 20)
PP, TT = np.meshgrid(pcol, tcol)
UP = M.F(PP, TT, Rc, p)
UT = np.array([[M.G(x, y, p) for x in pcol] for y in tcol])
N = np.hypot(UP, UT); N[N == 0] = 1
fig, ax = plt.subplots(figsize=(6.2, 4.9))
ax.quiver(PP, TT, UP / N, UT / N, color="0.75", pivot="mid", scale=26, width=0.0035)
tt = np.linspace(0.84, 0.95, 400)
ax.plot([M.A_of(t, Rc, p) for t in tt], tt, color="tab:green", lw=1.9,
        label=r"$P=A(T,R)$")
ax.axvline(0.0, color="tab:green", lw=1.1)
ax.axvline(p["K"], color="tab:green", lw=1.1)
den = p["a"] / (1 + p["g"] * tt) - p["eta"]
Pn = (p["s0"] - p["mu"] * tt) / den
m = (Pn > -0.02) & (Pn < 0.46)
ax.plot(Pn[m], tt[m], color="tab:orange", lw=1.9, label=r"$T$-nullcline")
for e in M.equilibria(Rc, p):
    ax.plot(e["P"], e["T"], "o", ms=9.5,
            mfc=("tab:blue" if e["stable"] else "white"),
            mec=("tab:blue" if e["stable"] else "tab:red"), mew=2, zorder=5)
ax.set_xlabel(r"Corridor-scaled cover index  $P$")
ax.set_ylabel(r"Water-stress index  $T$")
ax.set_title(r"Phase portrait at $R=%.0f$ mm yr$^{-1}$" % Rc, fontsize=10)
ax.set_xlim(-0.02, 0.46); ax.set_ylim(0.84, 0.95)
ax.legend(loc="lower right", fontsize=8.5, frameon=False)
fig.tight_layout()
fig.savefig(OUT / "figure3_phase.pdf")
fig.savefig(OUT / "figure3_phase.png", dpi=150)
plt.close(fig)

# ---------------- Figure 4: forcing and observed state ----------------
fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.0))
ax = axs[0]
ax.axhspan(Rg, Rr, color="0.90", zorder=0)
ax.plot(obs.year, obs.R, "o-", color="tab:blue", ms=4, lw=1.4)
ax.axhline(Rg, color="0.35", ls=":", lw=1)
ax.axhline(Rr, color="0.35", ls=":", lw=1)
ax.text(2000.3, Rg + 0.9, r"$R_{\rm deg}=%.0f$" % Rg, fontsize=8, color="0.35")
ax.text(2000.3, Rr + 0.9, r"$R_{\rm rec}=%.0f$" % Rr, fontsize=8, color="0.35")
ax.text(2014, 40.5, "bistable", fontsize=8, color="0.45")
ax.set_xlabel("Year"); ax.set_ylabel(r"Effective water input  $R$  (mm yr$^{-1}$)")
ax.set_title("(a) Forcing", fontsize=10)
ax = axs[1]
ax.plot(obs.year, obs.P, "o-", color="tab:green", ms=4, lw=1.4)
ax.set_xlabel("Year"); ax.set_ylabel(r"Corridor-scaled cover index  $P$", color="tab:green")
ax.tick_params(axis="y", labelcolor="tab:green")
ax.set_ylim(0, 0.45)
ax2 = ax.twinx()
ax2.plot(obs.year, obs.TVDI, "s--", color="tab:brown", ms=3.4, lw=1.1)
ax2.text(2013.5, 0.808, "model stress range %.2f-%.2f" % (Tb, Th),
         fontsize=7.5, color="tab:brown")
ax2.set_ylabel(r"Water-stress index  $T$ (TVDI)", color="tab:brown")
ax2.tick_params(axis="y", labelcolor="tab:brown")
ax2.set_ylim(0.80, 1.00)
ax.set_title("(b) Observed state", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "figure4_forcing_state.pdf")
fig.savefig(OUT / "figure4_forcing_state.png", dpi=150)
plt.close(fig)
print("figures written to", OUT)
