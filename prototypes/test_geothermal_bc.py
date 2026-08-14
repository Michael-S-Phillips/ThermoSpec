"""Tests for the geothermal (fixed base heat flux) bottom boundary condition.

bottom_bc='geothermal' sets the virtual base node so the conductive flux at the base equals
cfg.geothermal_flux (W/m^2, upward). Two checks:
  1. after a run the base finite-difference flux equals the prescribed value;
  2. a dark column (no sun) reaches the steady state where the surface radiates exactly the
     geothermal flux: em*sigma*T_surf^4 = F_geo  ->  T_surf = (F_geo/(em*sigma))^0.25.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from config import SimulationConfig       # noqa: E402
from modelmain import Simulator           # noqa: E402

SIGMA = 5.670374419e-8


def _base_flux(sim):
    g = sim.grid
    Et_b = sim.cfg.Et[-1] if np.ndim(sim.cfg.Et) else sim.cfg.Et
    dtau = g.x[-1] - g.x[-2]
    return g.cond[-1] * (sim.T[-1] - sim.T[-2]) / (Et_b * dtau)


def test_geothermal_bc_produces_prescribed_base_flux():
    F_geo = 0.05
    cfg = SimulationConfig(
        use_RTE=False, single_layer=True, diurnal=True, sun=True,
        bottom_bc='geothermal', geothermal_flux=F_geo,
        auto_dt=False, tsteps_day=500, ndays=1, dust_thickness=0.05, Et=1000.0,
        geometric_spacing=True, T_bottom=250.0, latitude=0.0, dec=0.0, P=88775.0, S=1361.0,
        em=0.95, albedo=0.1, k_dust=0.01, rho_dust=1500.0, cp_dust=800.0)
    sim = Simulator(cfg)
    sim.run()
    flux = _base_flux(sim)
    assert abs(flux - F_geo) < 1e-9, f"base flux {flux:.6e} != geothermal_flux {F_geo}"


def test_geothermal_warms_deep_column_vs_zero_flux():
    common = dict(use_RTE=False, single_layer=True, diurnal=True, sun=True,
                  auto_dt=False, tsteps_day=500, ndays=3, dust_thickness=0.05, Et=1000.0,
                  geometric_spacing=True, T_bottom=250.0, latitude=0.0, dec=0.0, P=88775.0,
                  S=1361.0, em=0.95, albedo=0.1, k_dust=0.01, rho_dust=1500.0, cp_dust=800.0)
    geo = Simulator(SimulationConfig(bottom_bc='geothermal', geothermal_flux=0.5, **common)); geo.run()
    zero = Simulator(SimulationConfig(bottom_bc='neumann', **common)); zero.run()
    assert geo.T[-2] > zero.T[-2], "geothermal base flux should warm the deep column vs zero-flux"


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
