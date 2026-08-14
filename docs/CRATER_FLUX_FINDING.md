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

## Not changed here
The new terrain code feeds the engine correctly and is validated independently (exact view-factor
reciprocity + closed-enclosure closure; flat-DEM facets identical to round-off; injection hook
bit-transparent). This note is a hand-off, not a fix.
