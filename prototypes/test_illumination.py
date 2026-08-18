"""Guard the direct-beam illumination path. The production bug: ShadowTester.illuminated_facets went
through trimesh's rtree ray engine, which silently returned no hits when libspatialindex failed to
load -> every facet read as shadowed, Q_direct=0, the whole scene froze at ~46 K. The fix routes it
through a dependency-free numpy first-hit ray caster. These tests assert illumination is non-zero for
a sun above the horizon and (where trimesh works) bit-identical to the old ray engine."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from crater import CraterMesh, ShadowTester       # noqa: E402

MESH = os.path.join(ROOT, "Roughness_files", "new_crater2.txt")


def _unit(v):
    v = np.asarray(v, float); return v / np.linalg.norm(v)


def test_sun_above_horizon_illuminates_facets():
    st = ShadowTester(CraterMesh(MESH))
    overhead = st.illuminated_facets(_unit([0.1, 0.1, 1.0]))
    assert np.count_nonzero(overhead > 0) > 0.5 * len(overhead), \
        "overhead sun must illuminate most facets (the bug left ALL facets shadowed)"
    low = st.illuminated_facets(_unit([1.0, 0.3, 0.09]))          # ~5 deg elevation
    assert np.count_nonzero(low > 0) > 0, "a low sun above the horizon must light some facets"
    assert np.count_nonzero(low > 0) < np.count_nonzero(overhead > 0), "low sun lights fewer than overhead"


def test_no_dependency_on_trimesh_ray():
    # ShadowTester must not build a trimesh ray index at all (that was the rtree failure point)
    st = ShadowTester(CraterMesh(MESH))
    assert not hasattr(st, "sub_mesh") and not hasattr(st, "mesh"), \
        "ShadowTester should no longer hold trimesh meshes for the ray engine"


def test_matches_trimesh_reference_when_available():
    try:
        import trimesh
        trimesh.creation.box().ray.intersects_first(np.array([[0, 0, 5.0]]), np.array([[0, 0, -1.0]]))
    except Exception:
        print("  (trimesh ray unavailable -- skipping bit-exact reference check)")
        return
    m = CraterMesh(MESH)
    st = ShadowTester(m)
    for s in (_unit([0.1, 0.1, 1.0]), _unit([1.0, 0.0, 1.0]), _unit([1.0, 0.3, 0.09])):
        sub = trimesh.Trimesh(vertices=m.sub_vertices, faces=m.sub_faces, process=False)
        n = m.sub_centroids.shape[0]
        idx = sub.ray.intersects_first(m.sub_centroids + 40.0 * s, np.tile(-s, (n, 1)))
        illum = (np.arange(n) == idx) & (np.dot(m.sub_normals, s) > 0)
        ref = np.array([np.sum(illum[m.sub_face_index[fi]]) for fi in range(len(m.centroids))], float)
        ref /= len(idx) / len(m.sub_face_index)
        assert np.max(np.abs(st.illuminated_facets(s) - ref)) < 1e-9, "numpy != trimesh reference"


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
