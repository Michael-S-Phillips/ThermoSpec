# Finding: crater engine double-counts direct sunlight on illuminated facets

Surfaced by the terrain integration/reduction test (a clean flat-DEM case). **This is in the
pre-existing crater radiative engine, not the new terrain code** (view_factors / topography /
the injection hook are validated correct). Flagged for the analysis session / original author,
because it affects absolute temperatures of *sunlit* terrain.

## Symptom
A **flat** DEM surface (identical facets, zero mutual view factors) run through `crater=True`
peaks at **465 K** at equatorial noon, versus the correct flat-plate subsolar value **388 K**
from the smooth (non-crater) model with identical parameters. `(465/388)^4 ≈ 2.06` — roughly a
factor-of-two excess flux. dt-independent (checked to dt=22 s); facets stay identical to 0 K, so
the mesh is fine.

## Root cause
`crater.compute_multiple_scattered_sunlight` returns, per facet,
`F_SCAT = F_sun·illum·cos + Σ_j VF_ij·G_j` — i.e. **total incident solar = direct + scattered**,
not the scattered component alone. Instrumented on the flat facet: `Q_direct = 1224.9` (correct),
`Q_selfheat = 0` (correct), but `Q_scattered = 1361 = F_sun` even though every view factor is 0.

The surface BC (`modelmain.py:1007`) then sums
`Q_dir + Q_scat·(1-albedo)·π + Q_selfheat·emissivity·π`. Since `Q_scat` already contains the
direct beam, the **direct beam is counted twice** (once in `Q_dir`, once inside `Q_scat`), and
carries an extra `π`. That inflates directly-lit facets.

## Scope / who is affected
- **Directly-illuminated facets (sunlit slopes): inflated** — the terrain model's daytime lit
  temperatures are too high until this is resolved.
- **Shadowed facets (PSR floors — the science target): `Q_dir = 0`,** so `F_SCAT` is pure
  scattered-from-walls. The double-count term vanishes; the floor energy balance is
  scattered+self-heat+(earthshine)+geothermal as intended. So the PSR-floor result is likely
  **unaffected** — but the `π` factors in the BC assembly still warrant an independent
  radiometric check before trusting absolute floor temperatures.

## Suggested fix (for review — do not apply without the author's sign-off)
Either (a) have `compute_multiple_scattered_sunlight` return only the *scattered* increment
(subtract the direct term `F_sun·illum·cos`), so `Q_scat` is purely inter-facet scattered light;
or (b) drop the separate `Q_dir` term and let `Q_scat·(1-albedo)` be the single "absorbed total
incident solar" source — but then reconcile the `π` and the `(1-albedo)` vs `Q_dir`'s own
`(1-albedo)`. The clean flat-DEM case (must give 388 K) and a two-facet analytic scattering case
are good regression targets for whichever convention is chosen.

## RESOLUTION (2026-08-14, fix approved by CS and applied)
`compute_multiple_scattered_sunlight` now subtracts the direct beam (`F_sun*illum*cos`) before
returning, so `Q_scat` is the **purely-scattered** increment. The BC's separate `Q_dir` is the
only direct-beam term.

**Validation (all green):**
- **Flat DEM reduces to the smooth model EXACTLY:** crater flat `T_surf_crater` max = smooth
  (non-crater) flat `T_surf` max = **381.23 K, difference 0.000 K** (both are the diurnal peak;
  381 < 388 because thermal inertia lowers the peak below the inertia-free instantaneous
  equilibrium). Before the fix the crater flat was 456 K. (`test_terrain_integration.py::
  test_flat_dem_crater_reduces_to_smooth_model`)
- **Two-facet analytic scattering** reproduces `A vf F_sun c/(1-A^2 vf^2)` for the scattered term
  (confirms the scattered physics survives the subtraction), and a no-view-factor facet gets
  exactly zero scattered light. (`test_crater_scattering.py`)

**Baseline shift on the hemispherical `new_crater2` crater=True run** (dt=44 s, equatorial diurnal,
non-RTE, albedo 0.1, em 0.95):

| quantity | before fix | after fix |
|---|---|---|
| T_surf min | 185.77 K | 185.67 K |
| T_surf max | 499.43 K | 418.34 K |
| T_surf mean | 291.32 K | 260.82 K |

Per-facet shift: mean **−30.5 K**, max |Δ| **83.7 K** (the most strongly self-illuminated /
sunlit-tilted facets fall most; shadow-side facets barely move). Any prior published `crater=True`
lit-facet results shift by this magnitude.

**Still open (CS-flagged, does not block the fix):** the `π` and `(1-albedo)` factors in the BC
assembly `Q_dir + Q_scat*(1-A)*π + Q_selfheat*em*π`. The flat case cannot test these (they multiply
zero terms). A facet-pair case with a known analytic *equilibrium* (absorbed = emitted) is the way
to confirm the scattered/self-heat radiometry; deferred as a follow-up.

## New terrain code (unchanged, validated independently)
View-factor reciprocity + closed-enclosure closure; flat-DEM facets identical to round-off;
injection hook bit-transparent. The fix above is in the pre-existing engine, not the terrain code.
