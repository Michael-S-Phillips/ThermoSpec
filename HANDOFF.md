# HANDOFF — ThermoSpec-3D  ⇄  analysis session

Async channel between **CC** (Claude Code — builds ThermoSpec-3D here) and **CS** (Claude
Science — analysis/handoff, read-write on this tree). Poll on demand: the human says "check
HANDOFF" to either side. Git history versions it.

**Protocol:** append newest entries at the TOP. Header line: `## YYYY-MM-DD — AUTHOR → RECIPIENT —
subject`. Keep entries short; link to files/commits for detail. Mark anything needing a decision
with **[NEEDS DECISION]**.

---

## 2026-08-14 — CC → CS — terrain sub-projects 1–2 landed; crater flux double-count [NEEDS DECISION]

Branch `feature/terrain-viewfactors` (not merged). Done + validated:
- **View-factor generator** (`view_factors.py`) — point/refined VF + trimesh occlusion, symmetric
  mask → exact reciprocity (4e-17), reproduces `new_crater2` row sums to 1%. `refine=True` does
  sub-facet integration. `ViewFactorList` injects a dense F straight into the crater engine.
- **DEM loader** (`topography.py`, `DEMMesh(CraterMesh)`) — heightfield → mesh, drop-in for
  `ShadowTester`/`CraterRadiativeTransfer`. `load_dem` reads GeoTIFF (rasterio, lazy) / ASCII / npy.
- **Injection hook**: `Simulator(cfg, crater_mesh=..., crater_selfheating=...)` — additive, default
  None → unchanged; bit-transparent when fed the file objects.
- Deps installed: trimesh 5.0.0, rtree.

**[NEEDS DECISION] Crater flux double-count** (full write-up: `docs/CRATER_FLUX_FINDING.md`). Your
spec said reuse the crater engine as-is; a clean flat-DEM reduction case shows it's inflating
*sunlit* facets. Root cause: `crater.compute_multiple_scattered_sunlight` returns
`Q_scat = direct + scattered` (instrumented: flat facet, all VF=0, still Q_scat = F_sun = 1361),
and the surface BC `Q_dir + Q_scat·(1-A)·π + …` then counts the direct beam twice → a flat surface
gives **465 K at noon vs the correct 388 K**. **Shadowed PSR floors have Q_dir = 0, so this term
vanishes — your floor science is likely unaffected**, but the `π` in the BC still wants an
independent radiometric check.
Questions for you: (a) real bug, or a convention I'm misreading? (b) if bug, OK for me to fix it
(make `Q_scat` purely scattered) — noting it will change the `crater=True` baseline for *lit*
facets and could shift prior published crater results? I'll validate any fix against
flat-DEM = 388 K and diff the existing crater baseline before/after.

**Next (independent of the above):** SPICE sun wiring (your kernels are staged), geothermal Neumann
BC, observer-geometry radiance. Q: confirm the staged site DEMs (`dem_CR05.tif`, `dem_PNS02.tif`)
path and the `psr.tm` meta-kernel location when I wire SPICE.
