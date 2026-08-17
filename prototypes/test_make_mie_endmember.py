"""Tests for the Mie endmember generator (Preprocessing/make_mie_endmember.py) and the labradorite
plagioclase endmember it produced. Physics sanity + reproducibility + drop-in DISORT compatibility."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "Preprocessing"))

import make_mie_endmember as mk          # noqa: E402

OPT = os.path.join(ROOT, "Optical_props")
NK = os.path.join(OPT, "incoming_labradorite", "labradorite_nk_oriented.txt")
BOUNDS = os.path.join(OPT, "enst_300K_wn_bounds.txt")
PLAG = os.path.join(OPT, "plag_labradorite_300K_mie_combined.txt")


def _regen():
    d = np.loadtxt(NK)
    n_iso, k_iso = mk.orientation_average(d[:, 2:5], d[:, 5:8], "eps")
    edges = np.loadtxt(BOUNDS)
    centers = 0.5 * (edges[:-1] + edges[1:])
    n_g = mk.bin_average(d[:, 0], n_iso, edges)
    k_g = mk.bin_average(d[:, 0], k_iso, edges)
    lam = 1.0e4 / centers
    g, Cext, Csca, ssalb = mk.mie_table(lam, n_g, k_g, 14.0)
    tab = np.column_stack([lam, g, Cext, Csca, ssalb])
    return tab[tab[:, 0].argsort()]


def test_mie_absorption_sign_and_limits():
    # non-absorbing sphere -> ssalb=1; strongly absorbing -> ssalb<1 (correct miepython sign)
    lam = np.array([12.0, 12.0])
    g, Cext, Csca, ssalb = mk.mie_table(lam, np.array([1.5, 1.5]), np.array([0.0, 0.5]), 14.0)
    assert abs(ssalb[0] - 1.0) < 1e-6, f"k=0 ssalb {ssalb[0]}"
    assert ssalb[1] < 0.7, f"k=0.5 ssalb should be well below 1: {ssalb[1]}"
    assert np.all(Cext >= Csca - 1e-9)


def test_orientation_average_isotropic_is_identity():
    # identical principal components -> average equals the component
    n = np.full((5, 3), 1.6); k = np.full((5, 3), 0.2)
    ni, ki = mk.orientation_average(n, k, "eps")
    assert np.allclose(ni, 1.6) and np.allclose(ki, 0.2)


def test_committed_table_matches_regeneration():
    # the committed endmember must be reproducible from the staged n,k (guards accidental edits)
    committed = np.loadtxt(PLAG)
    regen = _regen()
    assert committed.shape == regen.shape == (916, 5)
    assert np.max(np.abs(committed - regen)) < 1e-5, "committed table != regenerated"


def test_table_physics_and_grid():
    t = np.loadtxt(PLAG)
    lam, g, Cext, Csca, ssalb = t.T
    assert t.shape == (916, 5)
    assert 6.6 < lam.min() and lam.max() < 25.1                 # enst thermal grid
    assert np.all((ssalb >= 0) & (ssalb <= 1.0000001))
    assert np.all(Cext >= Csca - 1e-9)
    assert np.all((g > -1) & (g < 1))
    # wn_bounds companion exists and matches enst grid length
    b = np.loadtxt(os.path.join(OPT, "plag_labradorite_300K_wn_bounds.txt"))
    assert len(b) == 917


def test_drop_in_disort_solves_and_differs_from_enstatite():
    from config import SimulationConfig
    from grid import LayerGrid
    from radiance3d import compute_spectra

    def cfg(mie, wnb):
        return SimulationConfig(
            use_RTE=True, RTE_solver='disort', thermal_evolution_mode='hybrid',
            output_radiance_mode='hybrid', single_layer=True, diurnal=True, sun=True,
            auto_dt=False, tsteps_day=1000, ndays=1, dust_thickness=0.5, Et=1000.0,
            radius=14e-6, fill_frac=0.63, geometric_spacing=True, bottom_bc='dirichlet',
            T_bottom=250.0, P=88775.0, S=1361.0, nstr=4, nmom=4, nstr_out=8, nmom_out=8,
            ssalb_therm=0.1, ssalb_vis=0.5, eta=1.0, g_therm=0.0, g_vis=0.0, R_base=0.0,
            use_spec=False, k_dust=7.4e-4, rho_dust=1100.0, cp_dust=825.0,
            mie_file=f'{OPT}/{mie}', mie_file_out=f'{OPT}/{mie}',
            wn_bounds=f'{OPT}/{wnb}', wn_bounds_out=f'{OPT}/{wnb}',
            substrate_spectrum=f'{OPT}/sabel_enstatite.txt',
            substrate_spectrum_out=f'{OPT}/sabel_enstatite.txt')

    g = LayerGrid(cfg('enst_300K_mie_combined.txt', 'enst_300K_wn_bounds.txt'))
    T = np.broadcast_to(np.linspace(340, 250, g.x_num), (2, 3, g.x_num)).copy()
    _, _, BTe = compute_spectra(cfg('enst_300K_mie_combined.txt', 'enst_300K_wn_bounds.txt'), g, T, observer_mu=1.0)
    _, _, BTl = compute_spectra(cfg('plag_labradorite_300K_mie_combined.txt', 'plag_labradorite_300K_wn_bounds.txt'), g, T, observer_mu=1.0)
    BTe = np.asarray(BTe)[0, 0]; BTl = np.asarray(BTl)[0, 0]
    assert np.all(np.isfinite(BTl)), "labradorite produced non-finite BT"
    assert np.nanmax(np.abs(BTl - BTe)) > 1.0, "labradorite spectrum should differ from enstatite"


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
