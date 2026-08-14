"""Tests for the DEM -> mesh loader (terrain sub-project 2).

DEMMesh must be a drop-in for CraterMesh: same attributes, so ShadowTester / CraterRadiativeTransfer
and the view-factor generator all consume it unchanged.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from crater import ShadowTester, CraterRadiativeTransfer   # noqa: E402
from view_factors import compute_view_factors               # noqa: E402
from topography import DEMMesh, load_dem                     # noqa: E402


def _bowl(n=21, R=10.0, depth=6.0):
    """Synthetic circular bowl heightfield centered in an n x n grid, dx=dy=1."""
    ax = np.arange(n) - (n - 1) / 2.0
    X, Y = np.meshgrid(ax, ax)
    d = np.sqrt(X**2 + Y**2)
    z = np.where(d < R, -depth * (1.0 - (d / R)**2), 0.0)
    return z


def test_flat_dem_geometry():
    E = np.zeros((5, 4))
    m = DEMMesh(E, dx=10.0, dy=10.0)
    n_expected = 2 * (5 - 1) * (4 - 1)
    assert len(m.faces) == n_expected
    assert np.allclose(m.normals[:, 2], 1.0), "flat DEM facet normals must point +z"
    assert np.allclose(m.areas, 0.5 * 10.0 * 10.0), "each triangle area = dx*dy/2"


def test_dem_mesh_is_crater_dropin_shadowtester():
    m = DEMMesh(np.zeros((6, 6)), dx=1.0, dy=1.0)
    # required CraterMesh attributes present
    for attr in ("normals", "areas", "centroids", "vertices", "faces",
                 "sub_vertices", "sub_faces", "sub_face_index", "sub_centroids",
                 "sub_normals", "tangent1", "tangent2"):
        assert hasattr(m, attr), f"missing {attr}"
    st = ShadowTester(m)
    illum = st.illuminated_facets(np.array([0.0, 0.0, 1.0]))   # overhead sun
    assert np.all(illum > 0.99), "flat DEM under overhead sun should be fully lit"


def test_bowl_view_factors_reciprocal_and_floor_sees_walls():
    m = DEMMesh(_bowl(), dx=1.0, dy=1.0)
    F = compute_view_factors(m, occlusion=True)
    A = m.areas
    assert np.max(np.abs(A[:, None] * F - (A[:, None] * F).T)) < 1e-10
    assert F.sum(1).max() <= 1.0 + 1e-9
    # the deepest (floor) facet should see other facets (walls): nonzero row sum
    floor = int(np.argmin(m.centroids[:, 2]))
    assert F[floor].sum() > 0.05, "bowl floor should have a view factor to the walls"


def test_bowl_casts_shadows_at_low_sun():
    m = DEMMesh(_bowl(), dx=1.0, dy=1.0)
    st = ShadowTester(m)
    low_sun = np.array([np.cos(np.radians(10.0)), 0.0, np.sin(np.radians(10.0))])  # 10 deg elev
    illum = st.illuminated_facets(low_sun)
    assert illum.min() < 0.5, "a bowl at 10 deg solar elevation must shadow some facets"
    assert illum.max() > 0.5, "sunward rim/walls should still be lit"


def test_load_dem_ascii_roundtrip():
    E = np.arange(12, dtype=float).reshape(3, 4)
    path = os.path.join(ROOT, "prototypes", "_dem_rt.asc")
    try:
        with open(path, "w") as fh:
            fh.write("ncols 4\nnrows 3\nxllcorner 0\nyllcorner 0\ncellsize 5\nNODATA_value -9999\n")
            for row in E:
                fh.write(" ".join(str(v) for v in row) + "\n")
        Eout, dx, dy = load_dem(path)
        assert np.allclose(Eout, E) and dx == 5.0 and dy == 5.0
    finally:
        if os.path.exists(path):
            os.remove(path)


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
