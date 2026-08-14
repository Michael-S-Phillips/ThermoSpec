# Design: facet view-factor generator (terrain sub-project 1 of 4)

Part of the real-topography PSR terrain-thermal effort (see the analysis handoff
`psr_topography_implementation_spec.md`). This is the first and riskiest piece: the repo can
only *read* facet view factors (`crater.SelfHeatingList`), nothing generates them. For a real
DEM-derived mesh we must compute them.

## Purpose
Given any triangular mesh with the `CraterMesh` interface (`normals`, `areas`, `centroids`,
`vertices`, `faces`), compute the facet-to-facet view-factor matrix `F[i,j]` (fraction of
radiation leaving facet i that reaches facet j), honoring line-of-sight occlusion by other
facets. Write it in the exact sparse text format `SelfHeatingList` already parses, so the whole
downstream radiative chain (`CraterRadiativeTransfer`, self-heating, multiple scattering) is
untouched.

## Interface (new module `view_factors.py`)
- `compute_view_factors(mesh, occlusion=True, lift=None) -> np.ndarray [N,N]`
  - `F[i,j]` = view factor i→j; diagonal 0.
- `write_selfheating_list(F, fname, threshold=1e-8)` — write one line per facet:
  `n idx1..idxn vf1..vfn` with **1-based** indices (matching `SelfHeatingList.__init__`, which
  subtracts 1), listing only j with `F[i,j] > threshold`.
- `view_factor_matrix_from_file(fname, N)` — thin convenience = `SelfHeatingList(fname).as_view_matrix(N)`.

## Algorithm
Point-to-point (facets as small patches at their centroids):
`F_ij = (cosθ_i · cosθ_j · A_j) / (π r²)`, where `r = |c_j - c_i|`, `cosθ_i = n_i·û_ij`,
`cosθ_j = n_j·(-û_ij)`, `û_ij` the unit vector c_i→c_j. Set 0 if `cosθ_i ≤ 0` or `cosθ_j ≤ 0`
(facets not facing each other).
- **Occlusion:** cast a ray from `c_i` lifted slightly along `n_i` toward `c_j`; F_ij=0 unless the
  first mesh intersection is facet j (trimesh `ray.intersects_first`, the same engine
  `ShadowTester` uses). Batched per source facet (N rays), O(N²) total — a one-time per-site
  precompute, not per timestep.
- `lift` defaults to a small fraction of the median edge length to avoid self-intersection.

## Validation (tests, `prototypes/test_view_factors.py`)
- **Reciprocity (exact, hard):** `A_i F_ij = A_j F_ji` to <1e-10 (true by construction if occlusion
  is symmetric; also checks the occlusion is symmetric).
- **Physical bound (hard):** every row sum `Σ_j F_ij ≤ 1 + 1e-9`.
- **Closed-enclosure closure:** an inward-faced convex mesh (trimesh icosphere) with no occlusion
  gives row sums ≈ 1 to the point-approximation discretization tolerance (empirically calibrated,
  a few %).
- **Reproduce `new_crater2` (tolerance):** on `Roughness_files/new_crater2.txt`, generated row sums
  match the reference `new_crater2_selfheating_list.txt` row sums (~0.49) to a few % — the
  reference is only self-reciprocal to 0.064, so per-element exactness is not expected; row sums and
  overall pattern are the comparison.

## Notes / scope
- Config crater paths point at `/Users/ryan/...`; the local mesh is `Roughness_files/new_crater2.txt`
  (repoint via kwargs/local path, as with the optical files — no `config.py` edit).
- Only the generator here; DEM→mesh loader, SPICE sun, and geothermal BC are later sub-projects.
- `trimesh 5.0.0` now installed. rasterio/spiceypy deferred to sub-projects 2/3.
