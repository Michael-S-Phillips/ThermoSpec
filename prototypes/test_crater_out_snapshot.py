"""cfg.crater_out_snapshot=True must reproduce the interpolated T_crater_out (to the <= dt/2 timing
difference) while storing NO full-depth history, and leave the surface outputs byte-identical."""
import os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); os.chdir(ROOT)
from config import SimulationConfig
from modelmain import Simulator
from topography import DEMMesh
from view_factors import compute_view_factors, ViewFactorList

def _bowl(nx=8, L=800.0, R=300.0, depth=60.0):
    dx = L / nx
    xx, yy = np.meshgrid(np.arange(nx) * dx - L/2, np.arange(nx) * dx - L/2)
    r = np.hypot(xx, yy)
    return np.where(r < R, -depth * (1 - (r/R)**2), 0.0), dx

def _run(snapshot):
    E, dx = _bowl(); m = DEMMesh(E, dx=dx, origin="centroid")
    F = compute_view_factors(m, occlusion=True)
    cfg = SimulationConfig(crater=True, use_RTE=True, RTE_solver='disort', thermal_evolution_mode='two_wave',
                           diurnal=True, last_day=True, auto_dt=False, tsteps_day=4000, ndays=2, freq_out=24,
                           history_stride=10, crater_out_snapshot=snapshot, latitude=np.radians(-80.0), dec=np.radians(1.5),
                           P=255144.3, T_bottom=90.0, illum_freq=5, compute_crater_radiance=False)
    sim = Simulator(cfg, crater_mesh=m, crater_selfheating=ViewFactorList(F)); sim.run()
    return sim

def test_snapshot_matches_interpolation_and_stores_no_depth_history():
    a = _run(False); b = _run(True)
    assert len(b.T_crater_history) == 0, "snapshot mode must not keep the depth history"
    assert len(a.T_crater_history) > 0
    assert np.array_equal(a.T_surf_crater_out, b.T_surf_crater_out), "surface outputs must be identical"
    assert a.T_crater_out.shape == b.T_crater_out.shape
    # the snapshot must be EXACTLY the field the interpolating run stored at the same step
    idx = np.asarray(a._hist_step_idx); H = np.stack(a.T_crater_history, 2)
    for j, ks in b._crater_snap_steps.items():
        h = np.where(idx == j)[0][0]
        for k in ks:
            assert np.array_equal(H[:, :, h], b.T_crater_out[:, :, k]), f"snapshot at step {j} != stored history"
    # and it can differ from the cubic interpolation by at most one step of evolution (the output time sits
    # up to dt/2 off the snap step; at a shadow turn-on the interpolant reports mid-jump, the snapshot post-jump)
    d = np.abs(a.T_crater_out - b.T_crater_out)
    win = H[:, :, idx >= b._hist_out_start / (a.t[1] - a.t[0])]
    one_step = np.abs(np.diff(win, axis=2)).max()
    print(f"  max|snapshot-interp| = {d.max():.3f} K (median {np.median(d):.4f}); max one-step change in window = {one_step:.3f} K")
    assert d.max() <= one_step + 1e-9
    assert np.all(b.T_crater_out > 0), "every output frame must have been snapshotted"

if __name__ == "__main__":
    fails = 0
    for name, fn in [(n, f) for n, f in globals().items() if n.startswith("test_") and callable(f)]:
        try:
            print(f"[RUN ] {name}"); fn(); print(f"[PASS] {name}")
        except Exception as e:
            fails += 1; print(f"[FAIL] {name}: {e}")
    sys.exit(1 if fails else 0)
