"""Regression guard: the vectorized crater self-heating / multiple-scattering (BLAS matmuls) must
stay bit-identical (to machine precision) to the original per-facet Python loops. Locks in the
nx=24 perf fix so a future refactor can't silently change the radiative coupling."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from crater import (CraterMesh, SelfHeatingList, CraterRadiativeTransfer,  # noqa: E402
                    compute_multiple_scattered_sunlight)

MESH = os.path.join(ROOT, "Roughness_files", "new_crater2.txt")
VF = os.path.join(ROOT, "Roughness_files", "new_crater2_selfheating_list.txt")
SUN = np.array([0.3, 0.4, 0.86603]); SUN = SUN / np.linalg.norm(SUN)


def _ref_mss(Alb, F_sun, illum, cos, vm, max_iter=100, tol=1e-5):
    n = len(illum); G = Alb * F_sun * illum * cos
    for _ in range(max_iter):
        Gn = np.zeros_like(G)
        for i in range(n):
            Gn[i] = Alb[i] * (F_sun * illum[i] * cos[i] + np.dot(vm[i], G))
        if np.allclose(Gn, G, rtol=tol, atol=tol):
            break
        G = Gn
    F = G.copy(); direct = F_sun * illum * cos
    F[Alb > 0.0] /= Alb[Alb > 0.0]; F[Alb > 0.0] -= direct[Alb > 0.0]
    return F


def _ref_selfheat(sh, therm_flux, n_facets, n_waves):
    Q = np.zeros((n_facets, n_waves))
    for i in range(n_facets):
        idxs = sh.indices[i]
        vfs = sh.view_factors[i] if therm_flux.ndim == 1 else sh.view_factors[i][:, None]
        Q[i] = np.sum(therm_flux[list(idxs)] * vfs, axis=0)
    return Q


def _setup():
    mesh = CraterMesh(MESH); sh = SelfHeatingList(VF)
    return mesh, sh, CraterRadiativeTransfer(mesh, sh), len(mesh.normals)


def test_multiple_scatter_matches_loop():
    mesh, sh, crt, N = _setup()
    for n_waves in (1, 3):
        cos = np.clip(mesh.normals @ SUN, 0, None)[:, None] * np.ones((1, n_waves))
        illum = (np.arange(N) % 3 == 0).astype(float)[:, None] * np.ones((1, n_waves))
        Alb = np.full((N, n_waves), 0.12)
        new = compute_multiple_scattered_sunlight(Alb.copy(), 1361.0, illum.copy(), cos.copy(), crt.view_matrix)
        ref = _ref_mss(Alb.copy(), 1361.0, illum.copy(), cos.copy(), crt.view_matrix)
        assert np.max(np.abs(new - ref)) < 1e-9, f"n_waves={n_waves}: {np.max(np.abs(new-ref)):.2e}"


def test_selfheat_matches_loop():
    mesh, sh, crt, N = _setup()
    for n_waves in (1, 3):
        therm = 0.95 * 5.67e-8 * (np.linspace(120, 260, N) ** 4)
        if n_waves > 1:
            therm = np.repeat(therm[:, None], n_waves, axis=1)
        _, _, Q_self, _ = crt.compute_fluxes(SUN, (np.arange(N) % 3 == 0).astype(float),
                                             therm, np.full(n_waves, 0.12), 0.95, 1361.0, n_waves=n_waves)
        ref = _ref_selfheat(sh, therm, N, n_waves)
        ref = ref[:, 0] if n_waves == 1 else ref
        assert np.max(np.abs(Q_self - ref)) < 1e-9, f"n_waves={n_waves}: {np.max(np.abs(Q_self-ref)):.2e}"


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
