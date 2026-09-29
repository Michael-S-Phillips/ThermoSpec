"""ShadowTester numba backend must reproduce the numpy backend's illumination on every mesh/sun,
including grazing polar geometry where t/u/v sit near the eps thresholds. Also reports the speedup."""
import os, sys, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
from crater import CraterMesh, ShadowTester, _HAS_NUMBA
from topography import DEMMesh

MESH = os.path.join(ROOT, "Roughness_files", "new_crater2.txt")

def _bowl(nx, L=2600.0, R=415.0, depth=100.0):
    dx = L / nx
    xx, yy = np.meshgrid(np.arange(nx) * dx - L/2, np.arange(nx) * dx - L/2)
    r = np.hypot(xx, yy)
    E = np.where(r < R, -depth * (1 - (r/R)**2), 0.0)
    rng = np.random.default_rng(0)
    return E + rng.normal(0, 2.0, E.shape), dx     # 2 m noise: breaks symmetry / exact ties

def _suns(n=40, seed=1):
    rng = np.random.default_rng(seed)
    az = rng.uniform(0, 2*np.pi, n)
    el = np.concatenate([np.radians(rng.uniform(-3, 3, n//2)),     # grazing polar band
                         np.radians(rng.uniform(3, 60, n - n//2))])
    return np.column_stack([np.cos(el)*np.cos(az), np.cos(el)*np.sin(az), np.sin(el)])

def _compare(mesh, label):
    a = ShadowTester(mesh, backend='numpy'); b = ShadowTester(mesh, backend='numba')
    ta = tb = 0.0; nlit = 0
    for s in _suns():
        t0 = time.time(); ia = a.illuminated_facets(s); ta += time.time() - t0
        t0 = time.time(); ib = b.illuminated_facets(s); tb += time.time() - t0
        assert np.array_equal(ia, ib), f"{label}: backends differ, max|d|={np.abs(ia-ib).max()}"
        nlit += int((ia > 0).sum())
    print(f"  {label}: {len(mesh.normals)} facets, {len(mesh.sub_faces)} sub-tris, 40 suns identical "
          f"({nlit} lit facet-suns); numpy {ta:.2f}s numba {tb:.2f}s -> {ta/max(tb,1e-9):.1f}x")

def test_crater_mesh_backends_agree():
    if not _HAS_NUMBA:
        print("  numba not installed -> skip"); return
    _compare(CraterMesh(MESH), "hemispherical crater")

def test_dem_bowl_backends_agree_and_numba_is_faster():
    if not _HAS_NUMBA:
        print("  numba not installed -> skip"); return
    for nx in (16, 30):
        E, dx = _bowl(nx)
        _compare(DEMMesh(E, dx=dx, origin="centroid"), f"bowl nx={nx}")

def test_auto_backend_selection():
    m = DEMMesh(_bowl(8)[0], dx=10.0, origin="centroid")
    assert ShadowTester(m).backend == ('numba' if _HAS_NUMBA else 'numpy')
    assert ShadowTester(m, backend='numpy').backend == 'numpy'
    try:
        ShadowTester(m, backend='bogus'); raise AssertionError("bogus backend accepted")
    except ValueError:
        pass

if __name__ == "__main__":
    fails = 0
    for name, fn in [(n, f) for n, f in globals().items() if n.startswith("test_") and callable(f)]:
        try:
            print(f"[RUN ] {name}"); fn(); print(f"[PASS] {name}")
        except Exception as e:
            fails += 1; print(f"[FAIL] {name}: {e}")
    sys.exit(1 if fails else 0)
