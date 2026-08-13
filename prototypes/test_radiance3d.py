"""Test 3D emissivity/brightness-temperature spectra against the 1D radiance processor.

The emergent thermal radiance is computed per column (radiation is vertical), so the 3D spectra
must equal ThermoSpec's own radiance_processor output column-for-column for the same temperature
profile. Uses the enstatite optics (hybrid thermal output; no solar file needed).
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import SimulationConfig                 # noqa: E402
from grid import LayerGrid                          # noqa: E402
from radiance3d import compute_spectra, band_brightness_temperature  # noqa: E402
from radiance_processor import calculate_radiances_from_results      # noqa: E402

OPT = os.path.join(ROOT, "Optical_props")


def _cfg():
    enst_mie = os.path.join(OPT, "enst_300K_mie_combined.txt")
    enst_wnb = os.path.join(OPT, "enst_300K_wn_bounds.txt")
    sub = os.path.join(OPT, "sabel_enstatite.txt")
    return SimulationConfig(
        use_RTE=True, RTE_solver='disort', thermal_evolution_mode='two_wave',
        output_radiance_mode='hybrid', single_layer=True, diurnal=True, sun=True,
        depth_dependent_properties=False, temperature_dependent_properties=False,
        auto_dt=False, tsteps_day=1000, ndays=1, dust_thickness=0.5, Et=1000.0,
        geometric_spacing=True, bottom_bc='dirichlet', T_bottom=250.0,
        P=88775.0, S=1361.0, nstr=4, nmom=4, nstr_out=16, nmom_out=16,
        ssalb_therm=0.1, ssalb_vis=0.5, eta=1.0, g_therm=0.0, g_vis=0.0, R_base=0.0,
        fill_frac=0.37, radius=14.0e-6, use_spec=False,
        k_dust=7.4e-4, rho_dust=1100.0, cp_dust=825.0,
        mie_file=enst_mie, mie_file_out=enst_mie, wn_bounds=enst_wnb, wn_bounds_out=enst_wnb,
        substrate_spectrum=sub, substrate_spectrum_out=sub,
    )


def _profile(nz):
    return np.linspace(360.0, 250.0, nz)               # a plausible daytime column


def test_3d_spectra_match_1d_radiance_processor():
    cfg = _cfg()
    g = LayerGrid(cfg)
    nz = g.x_num
    prof = _profile(nz)
    T_field = np.broadcast_to(prof, (2, 3, nz)).copy()

    wn, rad3d, BT3d = compute_spectra(cfg, g, T_field, observer_mu=1.0)

    # 1D reference: radiance_processor on the same single profile
    res = calculate_radiances_from_results(
        (prof[:, None], np.array([prof[1]]), np.array([0.0])), config=cfg,
        observer_mu=1.0, spectral_mode='hybrid', time_indices=[0], grid=g)
    rad1d = np.asarray(res['radiance_thermal'])[0, :, 0]      # [nwave]

    assert rad3d.shape == (2, 3, len(wn))
    err = np.max(np.abs(rad3d - rad1d[None, None, :]))
    scale = np.max(np.abs(rad1d))
    assert err / scale < 1e-9, f"3D spectra vs 1D radiance processor: rel err {err/scale:.2e}"


def test_phase_spectra_end_to_end():
    # Run a short DISORT 3D diurnal sim, record noon, and produce per-column noon spectra.
    from sim3d import Simulator3D
    sim = Simulator3D(_cfg(), nx=2, ny=1, dx_m=0.02, dy_m=0.02)
    sim.run(record_phases=True)
    wn, rad, BT = sim.phase_spectra('noon', observer_mu=1.0)
    assert rad.shape[:2] == (2, 1) and rad.shape[2] == len(wn)
    assert np.all(np.isfinite(rad)) and np.all(rad >= 0)
    finite_bt = BT[np.isfinite(BT)]
    assert finite_bt.size > 0 and finite_bt.max() > 200.0    # daytime noon is warm


def test_band_brightness_temperature_inverts_planck():
    # a blackbody at 300 K over a narrow band must invert back to 300 K
    wn_lo, wn_hi = np.array([900.0]), np.array([901.0])
    h, c, k = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    wn_m = 900.5 * 100.0
    B = (2 * h * c**2 * wn_m**3) / (np.exp(h * c * wn_m / (k * 300.0)) - 1.0)  # spectral
    L = B * (wn_hi - wn_lo) * 100.0                                            # band-integrated
    T = band_brightness_temperature(wn_lo, wn_hi, L)[0]
    assert abs(T - 300.0) < 0.05, f"BT inversion gave {T:.3f}"


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
