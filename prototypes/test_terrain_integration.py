"""Integration / reduction tests: DEM mesh + generated view factors through the crater engine.

The headline validation (spec section 5): generated view factors reproduce the file-based crater
thermal result, and a DEM-derived mesh runs the whole facet-coupled radiative + conduction engine.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import SimulationConfig                        # noqa: E402
from modelmain import Simulator                            # noqa: E402
from crater import CraterMesh, SelfHeatingList             # noqa: E402
from view_factors import compute_view_factors, ViewFactorList  # noqa: E402
from topography import DEMMesh                              # noqa: E402

MESH = os.path.join(ROOT, "Roughness_files", "new_crater2.txt")
VF = os.path.join(ROOT, "Roughness_files", "new_crater2_selfheating_list.txt")


def _cfg(**over):
    base = dict(
        crater=True, use_RTE=False, diurnal=True, sun=True, single_layer=True,
        crater_mesh=MESH, crater_selfheating=VF,
        auto_dt=False, tsteps_day=200, ndays=1, dust_thickness=0.05, Et=1000.0,
        geometric_spacing=True, bottom_bc='dirichlet', T_bottom=250.0,
        latitude=np.radians(0.0), dec=0.0, P=88775.0, S=1361.0,
        em=0.95, albedo=0.1, k_dust=0.01, rho_dust=1500.0, cp_dust=800.0)
    base.update(over)
    return SimulationConfig(**base)


def test_injection_hook_is_transparent():
    # Injecting the very objects the config would load must reproduce the file-path route exactly.
    a = Simulator(_cfg()); a.run()
    b = Simulator(_cfg(), crater_mesh=CraterMesh(MESH),
                  crater_selfheating=SelfHeatingList(VF)); b.run()
    assert np.max(np.abs(a.T_crater_out - b.T_crater_out)) < 1e-12


def test_generated_view_factors_agree_with_file_within_method_tolerance():
    # Reduction/consistency: generated VFs run through the SAME crater engine as the supplied file
    # VFs (so the engine's behavior cancels), isolating the view-factor method difference. Our VFs
    # are exactly reciprocal; the reference file is only self-reciprocal to 0.064, and uses a finer
    # near-field integration, so a few-K difference on this small deep crater is expected -- the
    # rigorous VF-correctness proof is the physics tests in test_view_factors.py (reciprocity,
    # closed-enclosure closure). dt-robust (~5 K at 89 s and 22 s alike; not a numerical artifact).
    m = CraterMesh(MESH)
    F = compute_view_factors(m, occlusion=True, refine=True)
    gen = Simulator(_cfg(), crater_mesh=m, crater_selfheating=ViewFactorList(F)); gen.run()
    ref = Simulator(_cfg()); ref.run()
    dT = np.max(np.abs(gen.T_surf_crater_out - ref.T_surf_crater_out))
    assert dT < 6.0, f"generated-VF vs file-VF crater temps differ by {dT:.3f} K"


def _bowl(n=15, R=7.0, depth=5.0):
    ax = np.arange(n) - (n - 1) / 2.0
    X, Y = np.meshgrid(ax, ax)
    d = np.sqrt(X**2 + Y**2)
    return np.where(d < R, -depth * (1.0 - (d / R)**2), 0.0)


def test_dem_mesh_drives_crater_engine():
    # End-to-end: a DEM-derived mesh + generated view factors drive the full crater conduction +
    # radiative engine and produce a finite temperature field for every facet.
    bowl = DEMMesh(_bowl(), dx=1.0, dy=1.0)
    F = compute_view_factors(bowl, occlusion=True)
    sim = Simulator(_cfg(tsteps_day=1000), crater_mesh=bowl,
                    crater_selfheating=ViewFactorList(F))
    sim.run()
    assert sim.T_surf_crater_out.shape[0] == len(bowl.normals)
    assert np.all(np.isfinite(sim.T_surf_crater_out))


def test_flat_dem_crater_reduces_to_smooth_model():
    # THE flux-fix validation: a flat DEM has zero mutual view factors, so after the
    # multiple-scattering fix (Q_scat purely scattered) the crater engine must reduce EXACTLY to
    # the smooth (non-crater) flat-surface model -- same physics, same surface temperature. Before
    # the fix the direct-beam double-count made the crater flat 456 K vs the smooth 381 K.
    common = dict(use_RTE=False, diurnal=True, sun=True, single_layer=True, auto_dt=False,
                  tsteps_day=2000, ndays=1, dust_thickness=0.05, Et=1000.0, geometric_spacing=True,
                  bottom_bc='dirichlet', T_bottom=250.0, latitude=np.radians(0.0), dec=0.0,
                  P=88775.0, S=1361.0, em=0.95, albedo=0.1, k_dust=0.01, rho_dust=1500.0, cp_dust=800.0)
    smooth = Simulator(SimulationConfig(crater=False, **common)); smooth.run()
    flat = DEMMesh(np.zeros((10, 10)), dx=1.0, dy=1.0)
    F = compute_view_factors(flat, occlusion=True)
    assert F.max() == 0.0, "coplanar flat facets must have zero mutual view factors"
    crat = Simulator(SimulationConfig(crater=True, crater_mesh=MESH, crater_selfheating=VF, **common),
                     crater_mesh=flat, crater_selfheating=ViewFactorList(F)); crat.run()
    T = crat.T_surf_crater_out
    assert np.all(np.ptp(T, axis=0) < 1e-6), "identical flat facets diverged"
    assert abs(T.max() - smooth.T_surf_out.max()) < 0.1, \
        f"crater flat {T.max():.2f} K != smooth flat {smooth.T_surf_out.max():.2f} K"


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
