# Prototype plan: 3D conduction in ThermoSpec (structured grid + ADI)

Status: **plan only, no code written.** Fork clone lives at
`/Users/phillipsm/Documents/Software/ThermoSpec-3D` (main @ 2173190, core files identical to
andrewryan2727 origin). Goal: extend conduction from 1D-per-column to structured 3D, reusing
~80% of the codebase.

## The governing insight (why this is tractable)

1. **RTE stays 1D-per-column.** In regolith the photon mean free path is microns, so radiation
   escapes vertically. `rte_disort.disort_run` / `rte_hapke.compute_source` compute a flux
   divergence along one vertical column (modelmain.py:848 / :841). In 3D you run that per
   surface pixel and inject the source into that column's top cells. **`rte_disort.py`,
   `rte_hapke.py`, all optical-property loading, and radiance/spectra post-processing do not
   change.**
2. **The conduction solve is one banded call.** The whole per-step conduction update is
   `self.T = self._fd1d_heat_implicit_diag(self.T, source)` → a single
   `scipy.linalg.solve_banded((1,1), self.grid.diag, b)` (modelmain.py:461, :873). Going 3D
   means replacing that one call with an ADI triple — three banded solves (z, x, y).
3. **ADI reuses the existing primitive.** The vertical stencil is already an implicit
   non-uniform tridiagonal solve (Kieffer 2013, `stencils.fd1d_heat_implicit_diagonal_
   nonuniform_kieffer`). Douglas–Gunn ADI turns a 3D implicit step into three 1D tridiagonal
   sweeps, each still `solve_banded`. So "3D conduction" becomes "two more tridiagonal sweeps,"
   not a new solver core.

## The two real risks (what the prototype exists to retire)

- **Wall-time.** Lunar cadence is ~100k steps/day × 8 days = 800k steps. A 3D ADI step is 3
  batched tridiagonal solves over Nx·Ny·Nz unknowns. The 1D baseline does ~10⁴ steps/s on 177
  nodes. Will a realistic 3D grid finish a converged lunation in acceptable wall time?
- **τ-vs-meter units.** The vertical coordinate is **optical depth τ**, not meters
  (`grid.x` is τ; depth_m = x/Et; diffusivity carries `Et²`, grid.py:415, :431). Lateral
  spacing is naturally in **meters**. The lateral conduction coefficients must be put in units
  consistent with the τ-scaled vertical ones. This is the main conceptual gotcha of the whole
  effort.

The prototype sidesteps unit-mixing (works in pure meters, no RTE) to isolate the ADI
correctness + wall-time questions first. Units are confronted only at integration (Phase 2).

---

## Phase 1 — standalone ADI prototype (throwaway; touches NO existing code)

New file `prototypes/adi3d_prototype.py`. Does not import ThermoSpec classes. Answers "does ADI
work and is it fast enough" in isolation.

**Domain/grid:** structured block, **all meters**, constant properties. Regolith values
`k=7.4e-4, rho=1100, cp=825` → K = k/(ρ·cp) ≈ 7.7e-10 m²/s. Lateral 1 m × 1 m, depth 0.5 m.
Start uniform in all three axes; Nz≈100 (match real vertical resolution), Nx=Ny≈32. Store T as a
`[Nz, Nx, Ny]` array.

**Scheme — Douglas ADI (3 sweeps/step):** each sweep is an implicit 1D solve along one axis for
all lines in the other two, i.e. a *batched* `solve_banded((1,1), ab, B)` where B stacks all
lines as columns (solve_banded already accepts a 2D RHS). Constant coefficients → build the
three direction matrices once. Reshape/transpose so the swept axis is leading, batch-solve,
reshape back. This is the single piece to get right.

**Validation — exact 3D decay mode.** On `[0,Lx]×[0,Ly]×[0,Lz]` with homogeneous Dirichlet
(T=0) walls,
`T(x,y,z,t) = exp(-K·λ·t)·sin(πx/Lx)·sin(πy/Ly)·sin(πz/Lz)`, `λ = π²(1/Lx²+1/Ly²+1/Lz²)`
is an exact solution. Initialize with the mode, march, compare to the analytic exponential.
Checks to pass:
- L2 error decreases at 2nd order in Δx and expected order in Δt;
- **unconditional stability** — take Δt far past the explicit limit and confirm it stays bounded
  and accurate (the reason for an implicit scheme at all);
- energy/symmetry sanity.

**Wall-time benchmark.** Time 10³–10⁴ steps at Δt≈25 s (P/100000) on Nx=Ny=32, Nz=100
(~10⁵ unknowns) and on 64². Report steps/s and the projected wall time for 800k steps. **Decision
gate:** acceptable (say ≲ a few hours per 8-lunation run at modest lateral resolution) → proceed
to Phase 2; otherwise exercise the mitigations below before integrating.

**Mitigation to test cheaply here (torch is already a dependency via pydisort):** the ADI sweeps
are batched tridiagonal solves — a natural GPU/MPS fit. Prototype a `torch` batched-tridiagonal
variant and time it on MPS (`torch.backends.mps.is_available()` = True on this machine). If CPU
wall-time is marginal, this is the escape hatch, and it costs little to measure now.

**Phase 1 deliverable:** a short report — ADI reproduces the analytic decay to expected order +
projected wall-time + whether MPS is needed. No integration until this passes.

---

## Phase 2 — grid generalization (biggest single piece)

`grid.py`. Add a `VolumeGrid` (or a `mode='3d'` branch on `LayerGrid`) that:
- **reuses `_build_layers` per column** for the vertical τ-grid (keep RTE compatibility), and
- adds lateral node spacings in **meters** plus per-face lateral conductances.
- **Reconcile units:** keep vertical in τ; express lateral conduction in meters and scale lateral
  coefficients by the appropriate `Et` factor so vertical (τ, ×Et²) and lateral (m) terms live in
  one consistent diffusion operator. Write this down explicitly and unit-test a single interior
  node's coefficients by hand.
- Generalize the property arrays (`cond, dens, heat, K`, grid.py:461-464) and depth-dependent /
  temperature-dependent property machinery (already per-node, grid.py:237-348) from 1D to 3D.
  These are per-cell operations and vectorize cleanly.

Regression anchor: with lateral spacing → ∞ (no lateral coupling) the 3D grid must reduce
exactly to the current 1D grid, column by column.

## Phase 3 — stencil + solver + stepping

- `stencils.py`: add lateral tridiagonal builders analogous to
  `fd1d_heat_implicit_diagonal_nonuniform_kieffer` for the x and y directions (non-uniform-safe),
  plus an `adi3d_step(T, sources, matsX, matsY, matsZ)` driver.
- `modelmain.py`: add `_adi3d_step` and switch on a new `cfg.conduction_3d` flag. Replace the
  single solve at **line 873** with the ADI triple *only when the flag is set*; the 1D path stays
  the default → **zero regression risk** for all existing runs. `check_and_update_temperature_
  dependent_properties` (line 869) rebuilds three direction matrices instead of one when
  properties drift.
- `_bc` (line 875): generalize the surface energy balance to the **top face** (a set of surface
  nodes). `_bc_noRTE` / `_T_surf_calc` (modelmain.py:408-442) already use array ops, so this is
  mostly reshaping the surface-node set; the bottom Dirichlet/Neumann generalizes trivially.

**Second-order risk to test HERE, not in Phase 1:** the surface BC is a *nonlinear* radiative
balance solved by Newton (`_T_surf_calc`). Operator-split (ADI) stepping with a strongly
nonlinear boundary can lose order or stability. Phase 1 (linear Dirichlet) will not surface this.
Test with a small 3D box + real nonlinear surface BC (still no RTE) before wiring RTE in.

## Phase 4 — RTE coupling (per column)

For each surface (i,j) column, call the existing 1D DISORT/Hapke solve and scatter its flux
divergence into that column's top cells; assemble `source_term` as `[Nz, Nx, Ny]`. The crater
path already runs N independent columns (`_fd1d_heat_implicit_diag_batch`, modelmain.py:465;
per-facet solvers with `n_cols`, rte_disort.py), so the per-column driver pattern exists to copy.
`rte_disort.py` / `rte_hapke.py` unchanged.

## Phase 5 — integration validation

1. **Reduces to 1D:** a 3D run with laterally-uniform properties and illumination must reproduce
   the 1D result column-by-column (lateral gradient ≈ 0). This is the integration regression test.
2. **Lateral conduction does something physical:** impose a lateral contrast (a shadowed strip,
   or a rock/regolith boundary) and confirm heat flows laterally and the near-boundary columns
   diverge from their isolated-1D counterparts by a sensible amount.
3. Compare against an independent 3D conduction code on a simple case if available (e.g. a
   `heat1d`-style extension, or an analytic wedge).

---

## Reuse vs. change summary

| Component | Fate |
|---|---|
| `rte_disort.py`, `rte_hapke.py`, optical loading, radiance/spectra post-proc | **unchanged** |
| diurnal driver, convergence checks, crater surface-radiation (shadow/view-factor/scatter) | **unchanged** |
| temp/depth-dependent property models | generalize per-cell (vectorized) |
| `grid.py` spatial layout | **major** — add 3D layout + lateral conductances + unit reconciliation |
| `stencils.py` | add lateral builders + ADI driver |
| `modelmain.py` solve call (line 873) + `_bc` (875) | swap 1D solve → ADI triple behind a flag; vectorize surface BC |

## Recommended immediate next step

Write **only** `prototypes/adi3d_prototype.py` (Phase 1). It is self-contained, retires the two
real risks (ADI correctness, wall-time), costs little, and touches nothing in the model. Decide on
the wall-time gate before investing in Phase 2. If wall-time is the problem, resolve it (coarser
lateral grid / matrix-factorization reuse / torch-MPS batched solves) before, not after, building
the integration.
