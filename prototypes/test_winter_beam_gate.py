"""Regression for the site-level solar-gate bug (CS finding 2026-09-04), a surviving instance of the
beam-dead bug class. modelmain gated all crater direct-beam on the site's flat-horizon elevation
(`self.F = sun_z > 0.001`); at ~89 S the Sun is never >~2.2 deg above horizontal, so for the ~52% of a
polar year the Sun is below the flat horizon that gate discarded the beam on EVERY sunward-tilted slope
(rim/upper walls) -- starving the shadowed floor of wall IR and pushing the winter (annual-minimum) floor
~20 K too cold. The fix gates PER FACET: a facet is lit when its own cos(incidence)=n.s>0, so the shadow
test runs whenever ANY facet faces the Sun. This guards that a Sun below the flat horizon still lights
sunward facets (which the old scalar gate would have zeroed).
"""
import os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); os.chdir(ROOT)
from topography import DEMMesh                 # noqa: E402
from crater import ShadowTester                # noqa: E402


def _bowl(nx=24, dx=125.0, depth=500.0, R_frac=0.45):
    ax = np.arange(nx) - (nx - 1) / 2.0
    X, Y = np.meshgrid(ax, ax); d = np.sqrt(X**2 + Y**2); R = R_frac * nx
    return np.where(d < R, -depth * (1.0 - (d / R)**2), 0.0), dx


def _sun(elev_deg, az_deg=35.0):
    e, a = np.radians(elev_deg), np.radians(az_deg)
    v = np.array([np.cos(e)*np.cos(a), np.cos(e)*np.sin(a), np.sin(e)])
    return v / np.linalg.norm(v)


def test_sun_below_flat_horizon_still_lights_sunward_facets():
    """A Sun 1 deg BELOW the flat horizon (old gate: F=0, whole mesh dark) must still light the
    sunward-tilted slopes: the per-facet gate (any n.s>0) fires and the shadow test lights some facets."""
    Z, dx = _bowl()
    m = DEMMesh(Z, dx=dx, dy=dx, origin="centroid")
    s = _sun(-1.0)
    old_scalar_gate = bool(s[2] > 0.001)                 # the buggy site-level gate
    new_perfacet_gate = bool(np.any(m.normals @ s > 1e-9))  # the fix
    assert old_scalar_gate is False, "test sun should be below the flat horizon (old gate would zero it)"
    assert new_perfacet_gate is True, "per-facet gate must fire: sunward slopes face a below-horizon Sun"
    lit = (ShadowTester(m).illuminated_facets(s) > 0).sum()
    assert lit > 0, "sunward slopes must be lit by a below-flat-horizon Sun (beam-dead regression)"


def test_deep_night_lights_nothing():
    """Sanity contrast: a Sun far below the horizon (steeper than any slope) lights nothing, and the
    per-facet gate correctly does NOT fire."""
    Z, dx = _bowl()
    m = DEMMesh(Z, dx=dx, dy=dx, origin="centroid")
    s = _sun(-60.0)
    assert not np.any(m.normals @ s > 1e-9), "no facet should face a Sun 60 deg below horizontal"


def test_summer_unchanged():
    """A Sun above the flat horizon lights facets under both gates (no regression on the summer path)."""
    Z, dx = _bowl()
    m = DEMMesh(Z, dx=dx, dy=dx, origin="centroid")
    s = _sun(+2.0)
    assert (s[2] > 0.001) and np.any(m.normals @ s > 1e-9)
    assert (ShadowTester(m).illuminated_facets(s) > 0).sum() > 0


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
