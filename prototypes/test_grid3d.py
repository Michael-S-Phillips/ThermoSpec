"""Tests for VolumeGrid: the integrated 3D conduction step built on ThermoSpec's own grid.

Correctness anchors:
  1. With lateral conduction OFF, the 3D LOD step must reproduce the existing 1D banded solve
     column-for-column (the z-sweep reuses LayerGrid.diag).
  2. A laterally-uniform field has zero lateral gradient, so with lateral conduction ON the step
     must STILL match the 1D solve per column -- lateral coupling can't corrupt the 1D case.
  3. The lateral operator conserves heat (Neumann walls) and relaxes a contrast toward uniform.
"""
import os
import sys

import numpy as np
from scipy.linalg import solve_banded

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import SimulationConfig      # noqa: E402
from grid import LayerGrid               # noqa: E402
from grid3d import VolumeGrid, _lateral_banded_neumann, build_vertical_diag  # noqa: E402


def _base():
    cfg = SimulationConfig(
        use_RTE=False, single_layer=True,
        depth_dependent_properties=False, temperature_dependent_properties=False,
        auto_dt=False, tsteps_day=1000, ndays=1,
        dust_thickness=0.05, Et=1000.0, geometric_spacing=False, bottom_bc='dirichlet',
        k_dust=5.5e-4, rho_dust=1100.0, cp_dust=825.0,
    )
    return cfg, LayerGrid(cfg)


def test_lateral_r_equals_physical_diffusivity():
    """Pin the tau conversion quantitatively: the lateral coefficient the code builds must equal
    the PHYSICAL r = dt * (k/(rho*cp)) / dx_m^2. A wrong Et power (the single most error-prone
    line) would be off by Et^2 = 1e6 and fail this, though the uniform-field reduce-to-1D tests
    would not notice (their lateral operator is identity for any r)."""
    cfg, g = _base()
    dx_m = 0.02
    vg = VolumeGrid(g, nx=5, ny=5, dx_m=dx_m, dy_m=dx_m, lateral_k=None)
    kmid = g.x_num // 2                                   # an interior (real) depth node
    r_used = (vg._abx[kmid][1, 2] - 1.0) / 2.0            # interior main diag = 1 + 2r
    r_phys = g.dt * (cfg.k_dust / (cfg.rho_dust * cfg.cp_dust)) / dx_m**2
    assert abs(r_used - r_phys) / r_phys < 1e-10, f"lateral r={r_used:.3e} vs physical {r_phys:.3e}"


def test_build_vertical_diag_matches_layergrid():
    """Pin build_vertical_diag against the 1D LayerGrid operator so the temp-dependent per-column
    stencil cannot silently drift from grid.py's formula."""
    cfg, g = _base()
    ab = build_vertical_diag(g.l_thick, g.dens, g.heat, g.cond, g.dt)
    assert np.max(np.abs(ab - g.diag)) < 1e-12


def test_z_sweep_matches_1d_when_lateral_off():
    cfg, g = _base()
    nz = g.x_num
    vg = VolumeGrid(g, nx=3, ny=4, dx_m=0.01, dy_m=0.01, lateral_k=0.0)
    rng = np.random.default_rng(0)
    col, src = rng.random(nz), rng.random(nz)
    T = np.broadcast_to(col, (3, 4, nz)).copy()
    S = np.broadcast_to(src, (3, 4, nz)).copy()
    out = vg.step(T, S)
    ref = solve_banded((1, 1), g.diag, col + g.dt * src)
    assert np.max(np.abs(out - ref[None, None, :])) < 1e-12


def test_lateral_uniform_field_matches_1d_with_lateral_on():
    cfg, g = _base()
    nz = g.x_num
    vg = VolumeGrid(g, nx=3, ny=4, dx_m=0.01, dy_m=0.01, lateral_k=None)  # isotropic lateral
    rng = np.random.default_rng(1)
    col = rng.random(nz)
    T = np.broadcast_to(col, (3, 4, nz)).copy()
    out = vg.step(T)
    ref = solve_banded((1, 1), g.diag, col)
    assert np.max(np.abs(out - ref[None, None, :])) < 1e-12


def test_lateral_operator_conserves_and_relaxes():
    n, r = 8, 0.3
    ab = _lateral_banded_neumann(n, r)
    T = np.zeros(n)
    T[:4] = 1.0                        # lateral step contrast, mean 0.5
    contrast0 = T.max() - T.min()
    for _ in range(200):
        T = solve_banded((1, 1), ab, T)
    assert abs(T.sum() - 4.0) < 1e-10, f"not conservative: sum={T.sum()}"
    assert (T.max() - T.min()) < 0.02 * contrast0, f"not relaxed: contrast={T.max()-T.min()}"
    assert np.allclose(T, 0.5, atol=0.02), f"not uniform: {T}"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failures = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except AssertionError as e:
            failures += 1
            print(f"FAIL {t.__name__}: {e}")
        except Exception as e:
            failures += 1
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failures}/{len(tests)} passed")
    sys.exit(1 if failures else 0)
