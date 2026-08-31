"""Regression for the equilibrium initial condition (HANDOFF 2026-08-29, CS steps 2-3).

`cfg.equilibrium_ic=True` initializes each column with the steady conductive profile that carries the
basal geothermal flux upward (dT/dz = geothermal_flux/k), instead of a uniform T_bottom. This is what
lets a permanently-shadowed cold-trap column START already draining F_geo, so the expensive spin-up does
not have to grow the whole geothermal gradient from flat. The two science-acceptance gates that policed
the old uniform-IC artifact (tools/check_science_gates.py G2/G3) must pass on the IC itself:

  G2 IC drainage   -- conductive flux through the dust cap approaches F_geo (was 20x under the bad IC)
  G3 monotone      -- the column rises monotonically with depth, peak at the base

Acceptance: with equilibrium_ic on, the node-gradient conductive flux equals F_geo everywhere interior
and the temperature peak is at the deepest node; with it off, the column is flat at T_bottom.
"""
import os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)                                       # optical-property relative paths resolve from repo root
from config import SimulationConfig                  # noqa: E402
from modelmain import Simulator                       # noqa: E402
F_GEO = 0.018


def _cfg():
    c = SimulationConfig()
    c.use_RTE = False
    c.crater = False
    c.temperature_dependent_properties = False
    c.bottom_bc = "geothermal"
    c.geothermal_flux = F_GEO
    c.T_bottom = 110.0
    c.dust_thickness = 0.10
    c.k_dust, c.rho_dust, c.cp_dust = 5.5e-4, 1100.0, 825.0
    c.k_rock, c.rho_rock, c.cp_rock = 2.0, 1500.0, 800.0
    c.rock_thickness = 1.0
    c.Et = 1000.0
    return c


def _profile(sim):
    Et = np.asarray(sim.cfg.Et, dtype=float)
    z = sim.grid.x / Et
    k = sim.grid.cond / (Et ** 2)
    return sim.T.copy(), z, k


def test_uniform_ic_is_flat():
    c = _cfg(); c.equilibrium_ic = False
    T, _, _ = _profile(Simulator(c))
    assert np.ptp(T) < 1e-6, f"uniform IC should be flat at T_bottom, got span {np.ptp(T):.3g} K"


def test_equilibrium_ic_drains_geothermal_and_is_monotone():
    c = _cfg(); c.equilibrium_ic = True
    T, z, k = _profile(Simulator(c))
    # G3: monotone rising, peak at the deepest node
    peak_frac = int(np.argmax(T)) / (len(T) - 1)
    assert peak_frac > 0.95, f"G3: profile peak at frac {peak_frac:.2f} (want >0.95)"
    assert np.all(np.diff(T) >= -1e-9), "G3: profile must be non-decreasing with depth"
    # G2: conductive flux == F_geo at every interior node (dT/dz = F_geo/k by construction)
    good = z > 1e-6
    q = np.abs(k * np.gradient(T, z))[good]
    assert np.allclose(q, F_GEO, rtol=0.02), \
        f"G2: interior conductive flux {q.min():.4f}..{q.max():.4f} W/m2 != F_geo={F_GEO}"
    # the dust cap (thin, low-k) carries the whole gradient: F_geo*L/k_dust
    expected_cap_dT = F_GEO * c.dust_thickness / c.k_dust
    surf_to_cap = T[z <= c.dust_thickness + 1e-9]
    assert abs(np.ptp(surf_to_cap) - expected_cap_dT) < 0.2 * expected_cap_dT, \
        f"cap dT {np.ptp(surf_to_cap):.2f} K != expected {expected_cap_dT:.2f} K"


def test_offset_zero_without_geothermal_bc():
    c = _cfg(); c.equilibrium_ic = True; c.bottom_bc = "neumann"
    T, _, _ = _profile(Simulator(c))
    assert np.ptp(T) < 1e-6, "equilibrium IC must be inert unless bottom_bc=='geothermal'"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failures = 0
    for t in tests:
        try:
            t(); print(f"PASS {t.__name__}")
        except AssertionError as e:
            failures += 1; print(f"FAIL {t.__name__}: {e}")
        except Exception as e:
            failures += 1; import traceback; print(f"ERROR {t.__name__}: {type(e).__name__}: {e}"); traceback.print_exc()
    print(f"\n{len(tests) - failures}/{len(tests)} passed")
    sys.exit(1 if failures else 0)
