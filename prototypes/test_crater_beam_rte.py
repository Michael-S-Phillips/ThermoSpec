"""Regression: the RTE (DISORT two_wave) crater direct beam must heat the column even when the
inter-facet SCATTERED term is zero.

Second beam-dead blocker (HANDOFF 2026-08-21): in the two_wave/hybrid crater path the visible DISORT
solver carries BOTH the direct beam (via mu_solar_facets + illuminated) and the scattered light
(Q=Q_scat), and `_bc` discards Q_dir for RTE ("already accounted for in the RTE solver"). The solver
call was gated on `np.any(Q_scat > 1e-2)`, so on a near-flat / coplanar crater (mutual view factors ~0
=> Q_scat = 0) the direct beam was dropped and sunlit facets froze cold despite Q_dir being hundreds
of W/m^2. The gate must open on illumination, not only on scattered light.

Setup here forces the exact trigger: a uniformly TILTED plane (all facets coplanar -> zero mutual view
factors -> Q_scat == 0), lit by an overhead injected sun. Pre-fix: frozen near T_bottom. Post-fix: the
sunlit facets warm to a few hundred K.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import SimulationConfig                        # noqa: E402
from modelmain import Simulator                            # noqa: E402
from topography import DEMMesh                              # noqa: E402
from view_factors import compute_view_factors, ViewFactorList  # noqa: E402
from crater import ShadowTester, CraterRadiativeTransfer   # noqa: E402


def _tilted_plane(n=8, dx=1.0, tilt_deg=20.0):
    # E = col * dx * tan(tilt): a single flat plane tilted about the row axis. All facets coplanar.
    j = np.arange(n)
    E = np.tile(j * dx * np.tan(np.radians(tilt_deg)), (n, 1))
    return E


def test_rte_crater_direct_beam_heats_when_scatter_is_zero():
    dx = 1.0
    mesh = DEMMesh(_tilted_plane(8, dx, 20.0), dx=dx, dy=dx, origin="centroid")
    nfac = len(mesh.normals)
    F = compute_view_factors(mesh, occlusion=True)
    assert F.max() < 1e-12, "precondition: coplanar tilted plane must have ~zero mutual view factors"

    # overhead sun -> lit, but coplanar so Q_scat stays 0 (the bug trigger)
    sun = np.array([0.0, 0.0, 1.0])
    illum = ShadowTester(mesh).illuminated_facets(sun)
    rt = CraterRadiativeTransfer(mesh, ViewFactorList(F))
    Q_dir, Q_scat, _, cosines = rt.compute_fluxes(sun, illum, np.zeros(nfac), 0.0, 1.0, 1366.0, 1)
    assert illum.sum() > 0 and Q_dir.max() > 100.0, "precondition: facets lit with a real direct beam"
    assert (Q_scat / np.pi).max() <= 1e-2, "precondition: scattered term must be below the old gate"

    cfg = SimulationConfig()
    cfg.use_RTE = True
    cfg.RTE_solver = "disort"
    cfg.thermal_evolution_mode = "two_wave"
    cfg.output_radiance_mode = "two_wave"
    cfg.mie_file = os.path.join(ROOT, "Optical_props", "enst_300K_mie_combined.txt")
    cfg.wn_bounds = os.path.join(ROOT, "Optical_props", "enst_300K_wn_bounds.txt")
    cfg.mie_file_out = cfg.mie_file
    cfg.wn_bounds_out = cfg.wn_bounds
    cfg.Et = 1000.0
    cfg.single_layer = True
    cfg.dust_thickness = 0.05
    cfg.k_dust, cfg.rho_dust, cfg.cp_dust = 5.5e-4, 1100.0, 825.0
    cfg.crater = True
    cfg.illum_freq = 10
    cfg.bottom_bc = "dirichlet"
    cfg.T_bottom = 100.0
    cfg.latitude = np.radians(0.0)
    cfg.dec = 0.0
    cfg.P = 88775.0
    cfg.auto_dt = False
    cfg.tsteps_day = 1500          # dt ~ 59 s (stable side of the surface radiative BC)
    cfg.ndays = 1
    cfg.freq_out = 24

    sv = np.tile(sun, (cfg.tsteps_day * cfg.ndays, 1))
    sim = Simulator(cfg, crater_mesh=mesh, crater_selfheating=ViewFactorList(F), sun_vectors=sv)
    sim.run()
    Tmax = np.nanmax(sim.T_surf_crater_out)
    assert Tmax > 250.0, \
        f"RTE crater direct beam did not heat the column (Tmax {Tmax:.1f} K, T_bottom 100) -- beam dropped"


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
