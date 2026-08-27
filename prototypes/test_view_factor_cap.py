"""Regression for the far-field view-factor guard (CS audit 2026-08-27): F = cos_i cos_j A_j/(pi r^2)
can give F_ij>1 / row-sum>1 for adjacent facets on steep/high-res meshes -> radiosity divergence.
The guard caps F_ij<=1 reciprocity-preservingly (no-op on well-resolved meshes) and warns on row-sum>1."""
import os, sys, warnings
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from view_factors import compute_view_factors                       # noqa: E402
from topography import DEMMesh                                      # noqa: E402


def _bowl(n=15, R=7.0, depth=4.0):
    ax = np.arange(n) - (n - 1) / 2.0
    X, Y = np.meshgrid(ax, ax); d = np.sqrt(X**2 + Y**2)
    return np.where(d < R, -depth * (1.0 - (d / R)**2), 0.0)


def test_production_like_mesh_unchanged_no_cap():
    # A well-resolved bowl (production regime): no pair should exceed F=1, guard is a no-op.
    m = DEMMesh(_bowl(), dx=1.0, dy=1.0, origin="centroid")
    with warnings.catch_warnings():
        warnings.simplefilter("error")                              # any guard warning -> failure
        F = compute_view_factors(m, occlusion=True)
    assert F.max() <= 1.0, f"F_ij>1 on a well-resolved mesh: {F.max()}"
    assert F.sum(1).max() <= 1.0 + 1e-9, f"row-sum>1 on a well-resolved mesh: {F.sum(1).max()}"


def test_steep_mesh_caps_F_at_one_and_preserves_reciprocity():
    # A deliberately steep, coarse mesh (tiny cells + huge relief) drives the far-field kernel past F=1.
    n = 6
    ax = np.arange(n) - (n - 1) / 2.0
    X, Y = np.meshgrid(ax, ax)
    Z = 40.0 * np.abs(X)                # near-vertical walls, cell size 1 -> r/sqrt(A) << 1
    m = DEMMesh(Z, dx=1.0, dy=1.0, origin="centroid")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        F = compute_view_factors(m, occlusion=True)
    assert F.max() <= 1.0 + 1e-12, f"guard failed to cap F_ij<=1: max {F.max()}"
    # reciprocity must survive the cap: A_i F_ij == A_j F_ji
    A = np.asarray(m.areas)
    resid = np.max(np.abs(A[:, None] * F - (A[:, None] * F).T))
    assert resid < 1e-9, f"cap broke reciprocity: residual {resid:.2e}"
    msgs = " ".join(str(x.message) for x in w)
    assert "capped" in msgs or "row-sum" in msgs, "steep mesh should emit a guard warning"


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
