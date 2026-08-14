"""Validate the crater multiple-scattering fix: Q_scat is the purely-scattered increment.

Two-facet analytic case (CS-requested): facet 0 sunlit, facet 1 shadowed, mutual view factor vf,
albedo A. Multiple scattering gives radiosity G0 = A(F_sun c + vf G1), G1 = A(vf G0), so the TOTAL
incident is G/A and the SCATTERED increment (what compute_multiple_scattered_sunlight must now
return) is (G/A - direct):
    scat_0 = F_sun c / (1 - A^2 vf^2) - F_sun c        (light scattered back from facet 1)
    scat_1 = A vf F_sun c / (1 - A^2 vf^2)             (facet 1 has no direct beam)
And a facet with no view factors must receive zero scattered light.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from crater import compute_multiple_scattered_sunlight   # noqa: E402


def test_two_facet_scattering_matches_analytic():
    A, F_sun, vf, c = 0.3, 1361.0, 0.4, 0.8
    Alb = np.full((2, 1), A)
    illum = np.array([[1.0], [0.0]])        # facet 0 lit, facet 1 shadowed
    cos = np.array([[c], [0.5]])            # facet 1 cosine irrelevant (illum=0)
    V = np.array([[0.0, vf], [vf, 0.0]])

    Fscat = compute_multiple_scattered_sunlight(Alb, F_sun, illum, cos, V)
    denom = 1.0 - A**2 * vf**2
    scat_0 = F_sun * c / denom - F_sun * c
    scat_1 = A * vf * F_sun * c / denom
    # tolerance at the iterative solver's convergence level (tol=1e-5), not machine precision
    assert abs(Fscat[0, 0] - scat_0) < 0.05, f"lit facet scattered {Fscat[0,0]} vs {scat_0}"
    assert abs(Fscat[1, 0] - scat_1) < 0.05, f"shadow facet scattered {Fscat[1,0]} vs {scat_1}"


def test_no_view_factors_gives_zero_scattered():
    A, F_sun = 0.1, 1361.0
    Alb = np.full((2, 1), A)
    illum = np.array([[1.0], [0.0]])
    cos = np.array([[1.0], [0.0]])
    V = np.zeros((2, 2))
    Fscat = compute_multiple_scattered_sunlight(Alb, F_sun, illum, cos, V)
    assert np.allclose(Fscat, 0.0), f"no-VF facets got scattered light: {Fscat.ravel()}"


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
