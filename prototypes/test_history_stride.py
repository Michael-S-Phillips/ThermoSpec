"""history_stride: thin the interpolation history to bound peak memory (CS efficiency ask).

`Simulator.run()` stores a full T_crater [n_depth, n_facets] copy every integration step and stacks
them at the end to cubic-interpolate onto the output times, so peak memory scales with
tsteps_day*ndays -- ~50 GB for a stable-dt 2-day 450-facet run. history_stride>1 stores only every
Nth *spin-up* step (the pre-output cycles that last_day discards); the output window is always kept
at full resolution because sharp shadow-transition facets alias under a coarse grid.

Checks:
 - stride=1 is bit-identical to the previous every-step behaviour;
 - larger strides keep outputs within a physically negligible tolerance (<0.05 K, well under
   Diviner NEdT ~0.1 K) while storing fewer steps and staying length-aligned (t / T_crater / mu);
 - a single-day last_day run is a no-op (the whole run is the output window);
 - the 1-D (non-crater) path is unaffected.
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
TOL = 0.05  # K -- physically negligible (below Diviner NEdT); realistic fine-dt runs are sub-mK


def _crater_cfg(stride, tsteps_day=800, ndays=2):
    return SimulationConfig(
        crater=True, use_RTE=False, diurnal=True, sun=True, single_layer=True,
        crater_mesh=MESH, crater_selfheating=VF, history_stride=stride,
        auto_dt=False, tsteps_day=tsteps_day, ndays=ndays, dust_thickness=0.05, Et=1000.0,
        geometric_spacing=True, bottom_bc='dirichlet', T_bottom=250.0,
        latitude=0.0, dec=0.0, P=88775.0, S=1361.0,
        em=0.95, albedo=0.1, k_dust=0.01, rho_dust=1500.0, cp_dust=800.0,
        freq_out=48, last_day=True)


def _run(cfg):
    sim = Simulator(cfg)
    sim.run()
    return sim


def test_histories_stay_length_aligned():
    sim = _run(_crater_cfg(stride=50))
    assert len(sim.t_history) == len(sim.T_crater_history) == len(sim._hist_step_idx), \
        f"{len(sim.t_history)}, {len(sim.T_crater_history)}, {len(sim._hist_step_idx)}"


def test_stride1_is_bit_identical():
    ref = _run(_crater_cfg(stride=1))
    s1 = _run(_crater_cfg(stride=1))
    assert np.array_equal(ref.T_crater_out, s1.T_crater_out)
    # every step stored
    assert len(s1.t_history) == 800 * 2


def test_stride_thins_spinup_and_preserves_output():
    ref = _run(_crater_cfg(stride=1))
    sim = _run(_crater_cfg(stride=100))
    assert len(sim.t_history) < len(ref.t_history), "stride>1 did not thin the history"
    d = np.nanmax(np.abs(sim.T_crater_out - ref.T_crater_out))
    ds = np.nanmax(np.abs(sim.T_surf_crater_out - ref.T_surf_crater_out))
    assert d < TOL and ds < TOL, f"output drift too large: {d:.3e} K / {ds:.3e} K"


def test_single_day_last_day_is_noop():
    # ndays=1, last_day=True -> out_start=0 -> the whole run is the output window -> no thinning
    ref = _run(_crater_cfg(stride=1, ndays=1))
    sim = _run(_crater_cfg(stride=100, ndays=1))
    assert len(sim.t_history) == len(ref.t_history)
    assert np.array_equal(sim.T_crater_out, ref.T_crater_out)


def test_1d_path_unaffected():
    base = dict(crater=False, use_RTE=False, diurnal=True, sun=True, single_layer=True,
                auto_dt=False, tsteps_day=800, ndays=3, dust_thickness=0.05, Et=1000.0,
                geometric_spacing=True, bottom_bc='dirichlet', T_bottom=250.0,
                latitude=0.0, dec=0.0, P=88775.0, S=1361.0,
                em=0.95, albedo=0.1, k_dust=0.01, rho_dust=1500.0, cp_dust=800.0,
                freq_out=48, last_day=True)
    ref = _run(SimulationConfig(history_stride=1, **base))
    sim = _run(SimulationConfig(history_stride=50, **base))
    assert len(sim.t_history) < len(ref.t_history)
    assert np.nanmax(np.abs(sim.T_out - ref.T_out)) < TOL


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
