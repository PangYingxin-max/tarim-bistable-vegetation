"""Figure 5: phenomenological spatial pattern simulation (v9).

The equation solved here is the complex Swift-Hohenberg (CSH) equation

    dA/dt = mu A - (1 + i c1) (lap + k0^2)^2 A - (1 + i c3) |A|^2 A ,

which is the finite-wavenumber form of the Complex Ginzburg-Landau equation.
The plain CGLE that was used in the earlier version of the figure,

    dA/dt = mu A + (1 + i c1) lap A - (1 + i c3) |A|^2 A ,

has its maximum linear growth rate at k = 0, so it selects no wavelength: what
it produces is domain coarsening of the uniform state, and the apparent
"pattern" changes with the grid and with the seed.  The (lap + k0^2)^2 term
restores a preferred wavenumber k0, so the simulation selects a well-defined
wavelength that can be measured from the spectrum.  Near onset the envelope of
the CSH solution obeys the CGLE, which is why the phenomenological reading of
the figure is unchanged.

Integration is a Strang split: the nonlinear term is solved exactly, since the
ODE dA/dt = -(1 + i c3)|A|^2 A has the closed-form solution
A(t) = A(0) exp(-(1 + i c3)|A(0)|^2 t), and the linear term is applied exactly
in Fourier space.  The scheme is therefore unconditionally stable in the linear
part and does not need the small time steps of the explicit ETD1 scheme used
before.

Outputs (written to outputs/):
    figure5_patterns.png / .pdf   the four-panel figure
    figure5_patterns.npz          the fields, for archiving
    pattern_v9.json               measured wavelength, amplitude, convergence
"""
import io
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import paths

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

OUT = paths.OUT

# ---------------------------------------------------------------- parameters
MU = 0.50            # distance above onset
C1 = 0.16            # real part of the linear coefficient (as in the text)
C3 = -0.67           # imaginary part of the nonlinear coefficient
LAMBDA0 = 5.0        # preferred wavelength, dimensionless
K0 = 2.0 * np.pi / LAMBDA0
N = 256              # grid points per direction (2D) and in 1D
L = 100.0            # domain size, 20 wavelengths across
DT = 0.05
T2D = 200.0          # 2D integration time
T1D = 100.0          # 1D integration time
SEED = 20261018


def linear_operator(k2):
    """mu - (1 + i c1)(k0^2 - k^2)^2, the CSH linear symbol."""
    return MU - (1.0 + 1j * C1) * (K0 ** 2 - k2) ** 2


def nonlinear_step(A, dt):
    """Exact solution of dA/dt = -(1 + i c3)|A|^2 A over dt."""
    return A * np.exp(-(1.0 + 1j * C3) * np.abs(A) ** 2 * dt)


def run_2d():
    rng = np.random.default_rng(SEED)
    x = np.linspace(-L / 2, L / 2, N, endpoint=False)
    kx = 2.0 * np.pi * np.fft.fftfreq(N, d=x[1] - x[0])
    KX, KY = np.meshgrid(kx, kx)
    Lop = linear_operator(KX ** 2 + KY ** 2)
    expL = np.exp(Lop * DT)

    A = 0.05 * (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N)))
    nsteps = int(round(T2D / DT))
    snapshots = {}
    for n in range(nsteps + 1):
        if n % 1000 == 0:
            snapshots[int(round(n * DT))] = A.copy()
        A = nonlinear_step(A, 0.5 * DT)
        A = np.fft.ifft2(np.fft.fft2(A) * expL)
        A = nonlinear_step(A, 0.5 * DT)

    absA = np.abs(A)
    kmag = np.sqrt(KX ** 2 + KY ** 2)
    bins = np.linspace(0.0, kmag.max(), 400)
    idx = np.digitize(kmag.ravel(), bins)

    def radial(field):
        P = np.abs(np.fft.fft2(field)) ** 2
        p = np.bincount(idx, weights=P.ravel(), minlength=len(bins) + 1)
        c = np.bincount(idx, minlength=len(bins) + 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            p = p[:len(bins)] / np.maximum(c[:len(bins)], 1)
        return p[:-1]

    kc = 0.5 * (bins[:-1] + bins[1:])
    power = radial(A)
    power_100 = radial(snapshots[int(T2D - 100)])
    kpeak = float(kc[np.argmax(power)])
    kpeak_100 = float(kc[np.argmax(power_100)])
    # spectral change between t-100 and t: a robust stationarity measure that
    # is not polluted by the slow drift of individual defects
    rel_spec = float(np.sum(np.abs(power - power_100)) / np.sum(power_100))
    return dict(x=x, A=A, absA=absA, kc=kc, power=power, power_100=power_100,
                kpeak=kpeak, kpeak_100=kpeak_100, rel_spec=rel_spec,
                snapshots=snapshots, Lop=Lop, expL=expL)


def run_1d():
    rng = np.random.default_rng(SEED + 1)
    x = np.linspace(-L / 2, L / 2, N, endpoint=False)
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=x[1] - x[0])
    Lop = linear_operator(k ** 2)
    expL = np.exp(Lop * DT)

    A = 0.05 * (rng.standard_normal(N) + 1j * rng.standard_normal(N))
    nsteps = int(round(T1D / DT))
    every = 20
    field = []
    times = []
    for n in range(nsteps + 1):
        if n % every == 0:
            field.append(np.real(A).copy())
            times.append(n * DT)
        A = nonlinear_step(A, 0.5 * DT)
        A = np.fft.ifft(np.fft.fft(A) * expL)
        A = nonlinear_step(A, 0.5 * DT)
    power1d = np.abs(np.fft.fft(A)) ** 2
    return dict(x=x, times=np.array(times), field=np.array(field),
                final=np.real(A), power1d=power1d, k=k)


def main():
    print("2D CSH run, grid %d x %d, domain %.0f, %.0f time units" % (N, N, L, T2D))
    sys.stdout.flush()
    r2 = run_2d()
    print("  peak wavenumber k = %.4f  ->  wavelength %.3f" % (r2["kpeak"], 2 * np.pi / r2["kpeak"]))
    print("  peak at t-100: k = %.4f, spectral change %.4f" % (r2["kpeak_100"], r2["rel_spec"]))
    sys.stdout.flush()
    print("1D CSH run")
    sys.stdout.flush()
    r1 = run_1d()

    lam_meas = 2.0 * np.pi / r2["kpeak"]
    amp = float(np.mean(r2["absA"]))

    fig, ax = plt.subplots(2, 2, figsize=(10.5, 8.6))

    im = ax[0, 0].imshow(np.real(r2["A"]), cmap="RdBu_r", extent=[-L / 2, L / 2, -L / 2, L / 2],
                         origin="lower", vmin=-0.8, vmax=0.8)
    ax[0, 0].set_title(r"(a) real part $\mathrm{Re}\,A$, 2D", fontsize=11)
    ax[0, 0].set_xlabel("space X")
    ax[0, 0].set_ylabel("space Y")
    plt.colorbar(im, ax=ax[0, 0], fraction=0.046)

    ax[0, 1].plot(r2["kc"], r2["power"], lw=1.4, color="#1f4e79", label=r"$t=%.0f$" % T2D)
    ax[0, 1].plot(r2["kc"], r2["power_100"], lw=1.0, ls="--", color="#8faadc",
                  label=r"$t=%.0f$" % (T2D - 100))
    ax[0, 1].axvline(K0, ls="--", lw=1.0, color="crimson",
                     label=r"imposed $k_0=2\pi/\lambda_0$")
    ax[0, 1].axvline(r2["kpeak"], ls=":", lw=1.2, color="black",
                     label=r"measured peak, $\lambda=%.1f$" % lam_meas)
    ax[0, 1].set_xlim(0, 3 * K0)
    ax[0, 1].set_title(r"(b) radially averaged power spectrum", fontsize=11)
    ax[0, 1].set_xlabel("wavenumber $k$")
    ax[0, 1].set_ylabel("power")
    ax[0, 1].legend(fontsize=8, frameon=False, loc="upper right")

    im2 = ax[1, 0].imshow(r1["field"].T, aspect="auto", origin="lower", cmap="magma",
                          extent=[r1["times"][0], r1["times"][-1], -L / 2, L / 2])
    ax[1, 0].set_title(r"(c) 1D space--time diagram, $\mathrm{Re}\,A(X,t)$", fontsize=11)
    ax[1, 0].set_xlabel("time $t$")
    ax[1, 0].set_ylabel("space X")
    plt.colorbar(im2, ax=ax[1, 0], fraction=0.046)

    ax[1, 1].plot(r1["x"], r1["final"], lw=1.0, color="#1f4e79")
    ax[1, 1].set_xlim(-2.5 * LAMBDA0, 2.5 * LAMBDA0)
    ax[1, 1].set_title(r"(d) final 1D profile ($\lambda\approx%.1f$)" % lam_meas, fontsize=11)
    ax[1, 1].set_xlabel("space X")
    ax[1, 1].set_ylabel(r"$\mathrm{Re}\,A$")

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / ("figure5_patterns.%s" % ext), dpi=300, bbox_inches="tight")
    plt.close(fig)

    np.savez_compressed(OUT / "figure5_patterns.npz",
                        abs2d=r2["absA"], final1d=r1["final"],
                        x=r1["x"], times=r1["times"], field1d=r1["field"],
                        k=r2["kc"], power=r2["power"])

    diag = dict(
        equation="dA/dt = mu A - (1 + i c1)(lap + k0^2)^2 A - (1 + i c3)|A|^2 A",
        mu=MU, c1=C1, c3=C3, lambda0=LAMBDA0, k0=K0,
        grid=N, domain=L, dt=DT, t_2d=T2D, t_1d=T1D, seed=SEED,
        measured_wavelength=lam_meas,
        measured_peak_wavenumber=r2["kpeak"],
        wavelength_error_percent=100.0 * abs(lam_meas - LAMBDA0) / LAMBDA0,
        mean_amplitude=amp,
        std_amplitude=float(np.std(r2["absA"])),
        max_amplitude=float(r2["absA"].max()),
        kpeak_at_t_minus_100=r2["kpeak_100"],
        radial_spectrum_change_last_100=float(r2["rel_spec"]),
        benjamin_feir_1_plus_c1c3=1.0 + C1 * C3,
        note=("The plain CGLE of the earlier version has its growth maximum at "
              "k=0 and selects no wavelength; the (lap + k0^2)^2 term is what "
              "makes the pattern wavelength well defined and measurable."),
    )
    with open(OUT / "pattern_v9.json", "w", encoding="utf-8") as fh:
        json.dump(diag, fh, indent=2)
    print(json.dumps({k: v for k, v in diag.items() if k != "note"}, indent=2))
    print("wrote figure5_patterns.png/.pdf, figure5_patterns.npz, pattern_v9.json to", OUT)


if __name__ == "__main__":
    main()
