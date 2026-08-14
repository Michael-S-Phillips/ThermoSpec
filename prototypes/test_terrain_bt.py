"""Tests for the observer-geometry per-facet brightness-temperature helper (terrain_bt)."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import SimulationConfig                 # noqa: E402
from grid import LayerGrid                           # noqa: E402
from topography import DEMMesh                        # noqa: E402
from radiance3d import compute_spectra               # noqa: E402
from terrain_bt import TerrainObserver, terrain_bt_cube  # noqa: E402

OPT = os.path.join(ROOT, "Optical_props")


def _cfg():
    enst = os.path.join(OPT, "enst_300K_mie_combined.txt")
    wnb = os.path.join(OPT, "enst_300K_wn_bounds.txt")
    sub = os.path.join(OPT, "sabel_enstatite.txt")
    return SimulationConfig(
        use_RTE=True, RTE_solver='disort', thermal_evolution_mode='two_wave',
        output_radiance_mode='hybrid', single_layer=True, diurnal=True, sun=True,
        depth_dependent_properties=False, temperature_dependent_properties=False,
        auto_dt=False, tsteps_day=1000, ndays=1, dust_thickness=0.5, Et=1000.0,
        geometric_spacing=True, bottom_bc='dirichlet', T_bottom=250.0,
        P=88775.0, S=1361.0, nstr=4, nmom=4, nstr_out=4, nmom_out=4,
        ssalb_therm=0.1, ssalb_vis=0.5, eta=1.0, g_therm=0.0, g_vis=0.0, R_base=0.0,
        fill_frac=0.37, radius=14.0e-6, use_spec=False,
        k_dust=7.4e-4, rho_dust=1100.0, cp_dust=825.0,
        mie_file=enst, mie_file_out=enst, wn_bounds=wnb, wn_bounds_out=wnb,
        substrate_spectrum=sub, substrate_spectrum_out=sub)


def _bowl(n=9, R=5.0, depth=4.0):
    ax = np.arange(n) - (n - 1) / 2.0
    X, Y = np.meshgrid(ax, ax)
    d = np.hypot(X, Y)
    return np.where(d < R, -depth * (1.0 - (d / R)**2), 0.0)


def test_flat_nadir_matches_radiance3d():
    # Flat DEM at nadir: every facet's emission cosine is 1, so a facet's per-band radiance must
    # equal radiance3d.compute_spectra on that column at observer_mu=1 (the validated path).
    cfg = _cfg()
    g = LayerGrid(cfg)
    nz = g.x_num
    flat = DEMMesh(np.zeros((4, 4)), dx=1.0, dy=1.0)
    prof = np.linspace(360.0, 250.0, nz)
    T_facets = np.repeat(prof[:, None], len(flat.normals), axis=1)     # [nz, n_facets]

    obs = TerrainObserver(cfg, g, flat, mu_grid=np.linspace(0.1,1.0,4))                                # nadir default
    assert np.allclose(obs.observer_vec, [0, 0, 1])
    assert np.all(np.isclose(obs.mu_obs, 1.0)) and np.all(obs.visible)
    rad, BT = obs.brightness_temperature(T_facets)

    # compute_spectra needs >=2 columns (the DISORT n_cols=1 path is unsupported); use 2 identical.
    _, rad3d, _ = compute_spectra(cfg, g, np.repeat(prof[None, None, :], 2, axis=0), observer_mu=1.0)
    ref = np.asarray(rad3d)[0, 0]                                       # column 0 [nwave]
    err = np.max(np.abs(rad[0] - ref)) / np.max(np.abs(ref))
    assert err < 1e-6, f"facet radiance at mu=1 vs radiance3d: {err:.2e}"
    assert np.all(BT[obs.visible] > 200.0)                             # warm daytime column


def test_cube_shape_and_bands():
    cfg = _cfg()
    g = LayerGrid(cfg)
    nz = g.x_num
    flat = DEMMesh(np.zeros((4, 4)), dx=1.0, dy=1.0)
    nfac = len(flat.normals)
    T = np.repeat(np.linspace(350, 250, nz)[:, None], nfac, axis=1)[:, :, None]
    T = np.repeat(T, 3, axis=2)                                        # [nz, nfac, 3 times]
    out = terrain_bt_cube(cfg, g, flat, T, mu_grid=np.linspace(0.1,1.0,4))
    nb = len(out['wavenumbers'])
    assert out['BT'].shape == (nfac, nb, 3)
    assert out['wavelengths_um'].shape == (nb,)
    assert np.all(np.isfinite(out['BT']))


def test_observer_visibility_masks_hidden_facets():
    # A bowl viewed from a low oblique angle: some facets face away or are occluded -> BT is NaN.
    cfg = _cfg()
    g = LayerGrid(cfg)
    nz = g.x_num
    bowl = DEMMesh(_bowl(), dx=1.0, dy=1.0)
    T = np.repeat(np.linspace(300, 250, nz)[:, None], len(bowl.normals), axis=1)
    low_obs = np.array([np.cos(np.radians(20.0)), 0.0, np.sin(np.radians(20.0))])
    obs = TerrainObserver(cfg, g, bowl, observer_vec=low_obs, mu_grid=np.linspace(0.1,1.0,4))
    assert not np.all(obs.visible), "some bowl facets should be hidden from a low observer"
    _, BT = obs.brightness_temperature(T)
    assert np.all(np.isnan(BT[~obs.visible])), "hidden facets must be NaN"
    assert np.all(np.isfinite(BT[obs.visible]))


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
