"""Regression for the km-scale shadow-test leak (HANDOFF 2026-08-29): ShadowTester moved the ray
source a hardcoded 40 m along the sun, which is INSIDE a km-scale DEM crater, so rim occluders beyond
40 m (esp. at low sun elevation) were never tested and shadowed floor facets were falsely lit — the
"summer beam leak" that corrupted the sunlit-epoch PSR runs. Fix: scale the source offset to the mesh
bounding box. Acceptance: a deep-crater floor facet whose rim horizon exceeds the sun elevation must
get illuminated == 0.
"""
import os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from topography import DEMMesh                       # noqa: E402
from crater import ShadowTester, _sun_first_hit_numpy  # noqa: E402


def _deep_bowl(nx=24, dx=125.0, depth=600.0, R_frac=0.42):
    # a km-scale crater: nx*dx wide, `depth` m deep -> steep rim, floor deeply shadowed at low sun
    ax = np.arange(nx) - (nx - 1) / 2.0
    X, Y = np.meshgrid(ax, ax); d = np.sqrt(X**2 + Y**2); R = R_frac * nx
    return np.where(d < R, -depth * (1.0 - (d / R)**2), 0.0), dx


def _sun(elev_deg, az_deg=35.0):
    e, a = np.radians(elev_deg), np.radians(az_deg)
    v = np.array([np.cos(e)*np.cos(a), np.cos(e)*np.sin(a), np.sin(e)])
    return v / np.linalg.norm(v)


def test_low_sun_leaves_kmscale_floor_shadowed():
    # Robust invariants that separate the GROSS 40 m bug (deep floor fully lit) from the legitimate
    # fractional shadow-boundary. Measured: 40 m -> deepest-decile mean illum 0.60, 274 shadowed facets
    # fully lit; fix -> 0.00 and 0 fully lit (only ~17 partial-boundary facets remain).
    Z, dx = _deep_bowl()
    m = DEMMesh(Z, dx=dx, dy=dx, origin="centroid")
    ill = ShadowTester(m).illuminated_facets(_sun(4.0))   # sun far below the ~18 deg rim horizon
    sun_low = _sun(4.0)
    deepest = m.centroids[:, 2] < np.percentile(m.centroids[:, 2], 10)  # unambiguously below the rim
    assert ill[deepest].mean() < 0.05, \
        f"deep-floor facets falsely lit at 4 deg sun (mean illum {ill[deepest].mean():.2f}); ray source not outside the mesh"
    # no facet that is truly shadowed (forward-ray ground truth) may be FULLY lit
    cen, nrm = np.asarray(m.centroids), np.asarray(m.normals)
    truth_lit = (_sun_first_hit_numpy(np.asarray(m.vertices), np.asarray(m.faces),
                                      cen + 1e-3*dx*nrm, sun_low) < 0) & (nrm @ sun_low > 0)
    full_false = int(((ill > 0.9) & ~truth_lit).sum())
    assert full_false == 0, f"{full_false} truly-shadowed facets are FULLY lit (gross beam leak)"


def test_high_sun_lights_the_floor():
    # sanity contrast: a high sun (above the rim horizon) DOES light the floor
    Z, dx = _deep_bowl()
    m = DEMMesh(Z, dx=dx, dy=dx, origin="centroid")
    ill = ShadowTester(m).illuminated_facets(_sun(55.0))
    floor = m.centroids[:, 2] < np.percentile(m.centroids[:, 2], 20)
    assert (ill[floor] > 0).sum() > 0.5 * floor.sum(), "high sun should light most of the floor"


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
