"""Tests for the facet view-factor generator.

Exact physics checks (reciprocity, row-sum bound) plus calibrated geometry checks (closed-
enclosure closure, and reproduction of the reference new_crater2 view factors).
"""
import os
import sys

import numpy as np
import trimesh

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from crater import CraterMesh, SelfHeatingList          # noqa: E402
from view_factors import (                               # noqa: E402
    compute_view_factors, write_selfheating_list, view_factor_matrix_from_file)

CRATER = os.path.join(ROOT, "Roughness_files", "new_crater2.txt")
CRATER_VF = os.path.join(ROOT, "Roughness_files", "new_crater2_selfheating_list.txt")


class _Mesh:
    """Minimal CraterMesh-compatible mesh (normals, areas, centroids, vertices, faces)."""
    def __init__(self, tm, inward=False):
        self.vertices = np.asarray(tm.vertices, float)
        self.faces = np.asarray(tm.faces)
        self.centroids = np.asarray(tm.triangles_center, float)
        self.areas = np.asarray(tm.area_faces, float)
        n = np.asarray(tm.face_normals, float)
        self.normals = -n if inward else n


def test_reciprocity_and_row_sum_bound_on_crater():
    m = CraterMesh(CRATER)
    F = compute_view_factors(m, occlusion=True)
    A = m.areas
    resid = np.max(np.abs(A[:, None] * F - (A[:, None] * F).T))
    assert resid < 1e-10, f"reciprocity residual {resid:.2e}"
    assert np.all(np.diag(F) == 0.0)
    rs = F.sum(axis=1)
    assert rs.max() <= 1.0 + 1e-9, f"row sum exceeds 1: {rs.max()}"
    assert rs.min() > 0.0, "some facet sees nothing"


def test_closed_convex_enclosure_row_sums_near_one():
    # inward-faced icosphere: a closed convex enclosure -> each facet's radiation lands entirely
    # on other facets, so row sums -> 1 (to the point-approximation discretization error).
    ico = trimesh.creation.icosphere(subdivisions=3, radius=1.0)
    m = _Mesh(ico, inward=True)
    F = compute_view_factors(m, occlusion=False)
    rs = F.sum(axis=1)
    assert abs(rs.mean() - 1.0) < 0.03, f"closed enclosure mean row sum {rs.mean():.4f}"
    # convex -> occlusion removes nothing; occlusion path must agree
    Focc = compute_view_factors(m, occlusion=True)
    assert np.max(np.abs(F - Focc)) < 1e-9


def test_reproduces_new_crater2_row_sums():
    m = CraterMesh(CRATER)
    F = compute_view_factors(m, occlusion=True)
    ref = SelfHeatingList(CRATER_VF).as_view_matrix(len(m.normals))
    # the reference file is only self-reciprocal to 0.064, so compare row sums (sky-view deficit),
    # which are the physically meaningful, method-robust quantity.
    assert np.max(np.abs(F.sum(1) - ref.sum(1))) < 0.05, \
        f"row-sum mismatch {np.max(np.abs(F.sum(1)-ref.sum(1))):.3f}"


def test_numpy_and_trimesh_occlusion_backends_agree():
    # the default numpy Moller-Trumbore backend must match the trimesh ray backend exactly
    m = CraterMesh(CRATER)
    Fn = compute_view_factors(m, occlusion=True, occlusion_backend='numpy')
    Ft = compute_view_factors(m, occlusion=True, occlusion_backend='trimesh')
    assert np.max(np.abs(Fn - Ft)) < 1e-12


def test_numba_backend_matches_numpy_on_dem_bowl():
    # the numba grid-DDA occluder must be bit-identical to the numpy full scan (it's a broadphase +
    # the same Moller-Trumbore test), on a concave DEM bowl where occlusion actually happens.
    import view_factors as vf
    if not vf._HAS_NUMBA:
        print("  (numba not installed -- skipping)")
        return
    from topography import DEMMesh
    n = 20
    ax = np.arange(n) - (n - 1) / 2.0
    X, Y = np.meshgrid(ax, ax)
    d = np.hypot(X, Y)
    dem = np.where(d < n * 0.45, -4.0 * (1.0 - (d / (n * 0.45)) ** 2), 0.0)
    m = DEMMesh(dem, dx=10.0, dy=10.0)
    Fnp = compute_view_factors(m, occlusion=True, occlusion_backend='numpy')
    Fnb = compute_view_factors(m, occlusion=True, occlusion_backend='numba')
    assert np.max(np.abs(Fnp - Fnb)) == 0.0, f"numba != numpy: {np.max(np.abs(Fnp-Fnb)):.2e}"
    Fauto = compute_view_factors(m, occlusion=True, occlusion_backend='auto')
    assert np.array_equal(Fauto, Fnb), "auto backend should pick numba and match"


def test_write_reload_roundtrip(tmp_path=None):
    m = CraterMesh(CRATER)
    F = compute_view_factors(m, occlusion=True)
    out = os.path.join(ROOT, "prototypes", "_vf_roundtrip.tmp")
    try:
        write_selfheating_list(F, out)
        F2 = view_factor_matrix_from_file(out, len(m.normals))
        assert np.max(np.abs(F - F2)) < 1e-7
    finally:
        if os.path.exists(out):
            os.remove(out)


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
