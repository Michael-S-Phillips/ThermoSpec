"""Integration tests for Simulator3D (non-RTE diurnal conduction).

Decisive anchor: with laterally-uniform insolation the 3D run must reproduce the real 1D
`modelmain.Simulator` column-for-column to ~machine precision -- this simultaneously validates
the LOD stepping, the reused vertical operator, AND the vectorized nonlinear surface energy
balance (the Newton BC that the linear prototype could not exercise).

Physics anchor: with a permanently shadowed half, lateral conduction must warm the shadowed
columns relative to running them in isolation (lateral off).
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import SimulationConfig      # noqa: E402
from modelmain import Simulator          # noqa: E402
from sim3d import Simulator3D            # noqa: E402


def _cfg():
    return SimulationConfig(
        use_RTE=False, single_layer=True, diurnal=True, sun=True,
        depth_dependent_properties=False, temperature_dependent_properties=False,
        auto_dt=False, tsteps_day=500, ndays=2,
        dust_thickness=0.05, Et=1000.0, geometric_spacing=False,
        bottom_bc='dirichlet', T_bottom=250.0,
        latitude=0.0, dec=0.0, P=88775.0, S=1361.0, R=1.0,
        em=0.95, albedo=0.1, k_dust=0.01, rho_dust=1500.0, cp_dust=800.0,
    )


def test_laterally_uniform_matches_1d_simulator():
    cfg = _cfg()
    sim1d = Simulator(cfg)
    sim1d.run()

    sim3d = Simulator3D(_cfg(), nx=3, ny=2, dx_m=0.02, dy_m=0.02, lateral_k=None)
    sim3d.run()

    dT = np.max(np.abs(sim3d.T - sim1d.T[None, None, :]))
    dTs = np.max(np.abs(sim3d.T_surf - sim1d.T_surf))
    assert dT < 1e-6, f"3D field diverged from 1D by {dT:.2e} K"
    assert dTs < 1e-6, f"3D surface T diverged from 1D by {dTs:.2e} K"


def _run_match(cfg_kwargs):
    base = dict(
        use_RTE=False, single_layer=True, diurnal=True, sun=True,
        temperature_dependent_properties=False,
        auto_dt=False, tsteps_day=500, ndays=2, dust_thickness=0.05, Et=1000.0,
        bottom_bc='dirichlet', T_bottom=250.0, latitude=0.0, dec=0.0, P=88775.0,
        S=1361.0, R=1.0, em=0.95, albedo=0.1, k_dust=0.01, rho_dust=1500.0, cp_dust=800.0,
    )
    base.update(cfg_kwargs)
    sim1d = Simulator(SimulationConfig(**base))
    sim1d.run()
    sim3d = Simulator3D(SimulationConfig(**base), nx=3, ny=2, dx_m=0.02, dy_m=0.02)
    sim3d.run()
    return np.max(np.abs(sim3d.T - sim1d.T[None, None, :])), \
        np.max(np.abs(sim3d.T_surf - sim1d.T_surf))


def test_matches_1d_with_geometric_nonuniform_grid():
    dT, dTs = _run_match(dict(geometric_spacing=True, depth_dependent_properties=False))
    assert dT < 1e-6 and dTs < 1e-6, f"geometric grid: dT={dT:.2e}, dTs={dTs:.2e}"


def test_matches_1d_with_depth_dependent_properties():
    dT, dTs = _run_match(dict(
        geometric_spacing=True, depth_dependent_properties=True,
        rho_surface=1100.0, rho_deep=1800.0, rho_particle=3000.0,
        k_surface=7.4e-4, k_deep=3.4e-3, density_scale_height=0.06))
    assert dT < 1e-6 and dTs < 1e-6, f"depth-dependent: dT={dT:.2e}, dTs={dTs:.2e}"


def test_matches_1d_with_temperature_dependent_cp():
    # cp(T) drifts as the surface heats/cools, so the operator is rebuilt mid-run; the 3D
    # per-column operators must track the 1D single operator (identical for a uniform field).
    dT, dTs = _run_match(dict(
        geometric_spacing=True,
        temperature_dependent_properties=True, temp_dependent_cp=True,
        temp_change_threshold=1.0,
        cp_coeffs=[-3.6125, 2.7431, 2.3616e-3, -1.2340e-5, 8.9093e-9]))
    assert dT < 1e-5 and dTs < 1e-5, f"temp-dependent cp: dT={dT:.2e}, dTs={dTs:.2e}"


def test_lateral_conduction_warms_shadowed_columns():
    # nx=6 columns; the x<3 half is permanently shadowed (F_gate=0).
    def build(lateral_k):
        s = Simulator3D(_cfg(), nx=6, ny=1, dx_m=0.01, dy_m=0.01, lateral_k=lateral_k)
        s.F_gate[:3, :] = 0.0          # left half never sees the sun
        return s

    lit = build(lateral_k=None); lit.run()      # lateral conduction ON
    iso = build(lateral_k=0.0); iso.run()       # columns thermally isolated

    shadow_on = lit.T_surf[:3].mean()
    shadow_off = iso.T_surf[:3].mean()
    lit_on = lit.T_surf[3:].mean()
    lit_off = iso.T_surf[3:].mean()

    assert shadow_on > shadow_off + 0.5, (
        f"lateral conduction did not warm shadowed columns: {shadow_on:.2f} vs {shadow_off:.2f}")
    assert lit_on < lit_off, "lit columns should give up heat to the shadowed side"


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
            import traceback
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
            traceback.print_exc()
    print(f"\n{len(tests) - failures}/{len(tests)} passed")
    sys.exit(1 if failures else 0)
