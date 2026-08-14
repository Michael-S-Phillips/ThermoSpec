"""The crater/terrain fail-fast stability guard (CS ASK 2).

When a facet diverges at too-large dt its temperature leaves DISORT's valid range; the guard must
raise a clear one-line error naming the facet and step, instead of DISORT flooding millions of
temper-range warnings and the conduction solve then NaN-crashing.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from config import SimulationConfig       # noqa: E402
from modelmain import Simulator           # noqa: E402

MESH = os.path.join(ROOT, "Roughness_files", "new_crater2.txt")
VF = os.path.join(ROOT, "Roughness_files", "new_crater2_selfheating_list.txt")


def _sim():
    cfg = SimulationConfig(
        crater=True, use_RTE=False, diurnal=True, sun=True, single_layer=True,
        crater_mesh=MESH, crater_selfheating=VF,
        auto_dt=False, tsteps_day=500, ndays=1, dust_thickness=0.05, Et=1000.0,
        geometric_spacing=True, bottom_bc='dirichlet', T_bottom=250.0,
        latitude=0.0, dec=0.0, P=88775.0, S=1361.0,
        em=0.95, albedo=0.1, k_dust=0.01, rho_dust=1500.0, cp_dust=800.0)
    return Simulator(cfg)


def test_guard_passes_on_finite_in_range_field():
    sim = _sim()                              # T_crater initialised to T_bottom (finite, in range)
    sim._assert_crater_finite(1)              # must not raise


def test_guard_fails_fast_on_negative_temperature():
    sim = _sim()
    sim.T_crater[3, 7] = -5.0                 # a facet undershoots below 0 (DISORT temper error)
    try:
        sim._assert_crater_finite(99)
    except RuntimeError as e:
        msg = str(e)
        assert "facet 7" in msg and "step 99" in msg and "reduce dt" in msg, msg
    else:
        assert False, "guard did not fail fast on a negative facet temperature"


def test_guard_fails_fast_on_nan_temperature():
    sim = _sim()
    sim.T_crater[10, 42] = np.nan
    try:
        sim._assert_crater_finite(7)
    except RuntimeError as e:
        assert "facet 42" in str(e), str(e)
    else:
        assert False, "guard did not fail fast on a NaN facet temperature"


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
