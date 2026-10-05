"""Regression test for the two-layer grid's bottom ghost node (bug fixed 2026-10-05).
The ghost layer thickness was computed after the ghost node was appended, giving 2*(L - x_ghost) < 0; the
banded operator then coupled the last real node to the ghost with the WRONG sign (anti-diffusive), so a
geothermal/Neumann base acted as a heat sink of ~0.1 W/m2 and every two-layer column drained from the
bottom. Checks, for a dust-over-dust (dry control) and a dust-over-ice column:
  1. every layer thickness is positive and the ghost node mirrors the last real node about the edge;
  2. the operator's last-real-node coupling to the ghost is negative (diffusive), like every other coupling;
  3. pure conduction with an insulated top and a geothermal base GAINS energy at F_geo (within 20%) over
     a year and keeps the base warmer than the surface.
"""
import os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); os.chdir(ROOT)
from config import SimulationConfig
from modelmain import Simulator

def _cfg(ice):
    c = SimulationConfig()
    c.P = 2551443.0; c.use_RTE = True; c.RTE_solver = 'disort'; c.thermal_evolution_mode = 'two_wave'
    c.mie_file = 'Optical_props/enst_300K_mie_combined.txt'; c.wn_bounds = 'Optical_props/enst_300K_wn_bounds.txt'
    c.mie_file_out = c.mie_file; c.wn_bounds_out = c.wn_bounds
    c.Et = 1000.0; c.single_layer = False; c.dust_thickness = 0.05 if ice else 0.10
    c.k_dust, c.rho_dust, c.cp_dust = 5.5e-4, 1100.0, 825.0
    if ice:
        c.k_rock, c.rho_rock, c.cp_rock = 2.0, 1500.0, 800.0
    else:
        c.k_rock, c.rho_rock, c.cp_rock = 5.5e-4, 1100.0, 825.0
    c.rock_thickness = 1.0; c.temperature_dependent_properties = False
    c.bottom_bc = 'geothermal'; c.geothermal_flux = 0.018; c.T_bottom = 45.0; c.equilibrium_ic = True
    c.auto_dt = False; c.tsteps_day = 40000; c.ndays = 1; c.crater = False; c.latitude = np.radians(-89.3)
    return c

def _check(ice):
    sim = Simulator(_cfg(ice)); g = sim.grid
    L = (sim.cfg.dust_thickness + sim.cfg.rock_thickness) * 1000.0
    assert np.all(g.l_thick > 0), f"negative layer thickness: {g.l_thick[g.l_thick <= 0]}"
    assert abs((g.x[-1] - L) - (L - g.x[-2])) < 1e-6, "ghost node must mirror the last real node about the edge"
    assert abs(g.l_thick[-1] - 2 * (L - g.x[-2])) < 1e-6
    assert g.diag[0, -1] < 0, f"last-real-node -> ghost coupling must be diffusive (negative), got {g.diag[0,-1]:+.3e}"
    T = sim.T.copy(); T, _ = sim._bc(T, 45.0); z = np.zeros(g.x_num)
    cap = g.dens[1:-1] * g.heat[1:-1] * g.l_thick[1:-1] / 1000.0
    E0 = np.sum(T[1:-1] * cap); n = 40000 * 12
    for j in range(n):
        T = sim._fd1d_heat_implicit_diag(T, z); T, _ = sim._bc(T, 45.0)
    rate = (np.sum(T[1:-1] * cap) - E0) / (n * g.dt)
    print(f"  {'ice' if ice else 'dry'} column: ghost l_thick {g.l_thick[-1]:+.1f} tau, coupling {g.diag[0,-1]:+.2e}, "
          f"dE/dt {rate*1e3:+.1f} mW/m2 (F_geo 18), base-surf {T[-2]-T[1]:+.2f} K after 1 yr")
    assert 0.8 * 0.018 < rate < 1.2 * 0.018, f"column must gain energy at F_geo, got {rate:+.4f} W/m2"
    assert T[-2] > T[1], "base must stay warmer than the surface under an upward geothermal flux"

def test_dry_two_layer_column_bottom_is_a_source_not_a_sink():
    _check(ice=False)

def test_ice_two_layer_column_bottom_is_a_source_not_a_sink():
    _check(ice=True)

if __name__ == "__main__":
    fails = 0
    for name, fn in [(n, f) for n, f in globals().items() if n.startswith("test_") and callable(f)]:
        try:
            print(f"[RUN ] {name}"); fn(); print(f"[PASS] {name}")
        except Exception as e:
            fails += 1; print(f"[FAIL] {name}: {e}")
    sys.exit(1 if fails else 0)
