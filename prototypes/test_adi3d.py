"""Tests for the 3D ADI conduction solver (prototype).

The decisive test: a separable sine mode sin(pi x/Lx) sin(pi y/Ly) sin(pi z/Lz) sampled on a
uniform grid with Dirichlet-zero walls is an EXACT eigenvector of each 1D second-difference
operator. So the LOD/ADI step must multiply it by exactly 1/prod_d (1 + dt K mu_d), where
mu_d = (2/h_d^2)(1 - cos(pi/(n_d-1))) is the discrete eigenvalue. That gives a closed-form
prediction for the numerical field after N steps, checkable to machine precision -- it isolates
solver-implementation correctness from spatial discretization error.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from adi3d import ADI3DSolver  # noqa: E402


def _sine_mode(nx, ny, nz):
    ix, iy, iz = np.arange(nx), np.arange(ny), np.arange(nz)
    sx = np.sin(np.pi * ix / (nx - 1))
    sy = np.sin(np.pi * iy / (ny - 1))
    sz = np.sin(np.pi * iz / (nz - 1))
    return sx[:, None, None] * sy[None, :, None] * sz[None, None, :]


def _mu(h, n):
    return (2.0 / h**2) * (1.0 - np.cos(np.pi / (n - 1)))


def test_adi_matches_analytic_discrete_eigenmode():
    nx, ny, nz = 21, 17, 25
    Lx, Ly, Lz = 1.0, 0.8, 0.5
    hx, hy, hz = Lx / (nx - 1), Ly / (ny - 1), Lz / (nz - 1)
    K, dt, N = 7.7e-10, 25.0, 50

    T0 = _sine_mode(nx, ny, nz)
    solver = ADI3DSolver((nx, ny, nz), (hx, hy, hz), K, dt)

    f = 1.0 / ((1 + dt * K * _mu(hx, nx))
               * (1 + dt * K * _mu(hy, ny))
               * (1 + dt * K * _mu(hz, nz)))

    T = T0.copy()
    for _ in range(N):
        T = solver.step(T)

    predicted = f**N * T0
    err = np.max(np.abs(T - predicted)) / np.max(np.abs(T0))
    assert err < 1e-10, f"eigenmode decay mismatch: rel err {err:.2e}"


def test_unconditionally_stable_at_huge_dt():
    nx = ny = nz = 21
    T0 = _sine_mode(nx, ny, nz)
    solver = ADI3DSolver((nx, ny, nz), (0.05, 0.05, 0.05), K=7.7e-10, dt=1e9)
    T = T0.copy()
    for _ in range(100):
        T = solver.step(T)
    assert np.all(np.isfinite(T)), "solution went non-finite at huge dt"
    # a diffusion solver must never amplify a decaying mode
    assert np.max(np.abs(T)) <= np.max(np.abs(T0)) + 1e-12


def test_dirichlet_walls_stay_zero():
    nx, ny, nz = 12, 12, 12
    T0 = _sine_mode(nx, ny, nz)
    solver = ADI3DSolver((nx, ny, nz), (0.05, 0.05, 0.05), K=7.7e-10, dt=25.0)
    T = T0.copy()
    for _ in range(20):
        T = solver.step(T)
    for face in (T[0], T[-1], T[:, 0], T[:, -1], T[:, :, 0], T[:, :, -1]):
        assert np.max(np.abs(face)) < 1e-14


if __name__ == "__main__":
    # Runnable without pytest (not installed in the thermospec env).
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failures = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"FAIL {t.__name__}: {e}")
        except Exception as e:  # e.g. ImportError while solver missing
            failures += 1
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failures}/{len(tests)} passed")
    sys.exit(1 if failures else 0)
