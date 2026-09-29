# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

ThermoSpec-3D is a thermal + radiative-transfer model for planetary regolith. It solves the
time-dependent heat diffusion equation (implicit Euler) while, at each timestep, computing radiative
flux divergence with a full RTE solver and feeding it back as a source term. The novelty vs. ordinary
thermal models: the regolith surface is treated as **semi-transparent**, not an opaque boundary with a
fixed albedo/emissivity — sunlight penetrates and thermal emission arises from subsurface layers at
different temperatures, which reshapes the emergent spectrum (the "vacuum emissivity" effect).

Read `README.md` for the physics overview and `BATCH_PROCESSING_README.md` for the batch/ML-dataset
system. The model is alpha-stage and under active development; features may be incomplete or broken.

## Environment & running things

Requires a conda env with `pydisort` (the DISORT solver, imported at the top of `rte_disort.py`, so it
is needed even for `use_RTE=False` runs), plus PyTorch, numpy, scipy, h5py, pyyaml, and optionally
rasterio for GeoTIFF DEMs. On this machine that env is **`thermospec`** — **not `sentinel`, which lacks
pydisort.** macOS also needs `KMP_DUPLICATE_LIB_OK=TRUE` to work around a duplicate-libomp crash
(`OMP: Error #15`). Prefix Python invocations with:

```bash
source /Users/phillipsm/anaconda3/etc/profile.d/conda.sh && conda activate thermospec && \
  KMP_DUPLICATE_LIB_OK=TRUE python <script>
```

(The `/sentinel-python` slash command targets the `sentinel` env, which has no pydisort here — do not
use it for this repo; use the explicit `thermospec` activation above.)

- **1D model, quick run + plots:** `python modelmain.py` — instantiates `Simulator()` with `config.py`
  defaults and shows temperature/spectra plots.
- **Batch / single / post-process CLI:** `python run_thermal_batch.py {create-config,single,batch,postprocess}`
  — the YAML-config-driven path for parameter sweeps and ML datasets. See `BATCH_PROCESSING_README.md`
  for full command syntax and the HDF5 output schema.

## Tests

Tests are **standalone scripts, not pytest**. Each has an `if __name__ == "__main__"` block that
discovers its `test_*` functions, runs them, and exits nonzero on failure. Run one directly:

```bash
python prototypes/test_illumination.py
```

`prototypes/test_*.py` covers the 3D/roughness/terrain work; `test_*.py` at the repo root and `Tests/`
cover older 1D/interface pieces. Do not invoke `pytest`. Run tests from the repo root: optical-property
paths in `config.py` are relative to it. `prototypes/audit_*.py` are physics audits (radiosity energy
conservation, view-factor occlusion) — same standalone pattern, run them after touching `crater.py` or
`view_factors.py`.

## Architecture

### Configuration (single source of truth)
- `config.py` — `SimulationConfig` dataclass. **Every** simulation parameter lives here with inline
  documentation; this is the definitive reference for what each knob does.
- `core/config_manager.py` — loads/flattens YAML into a `SimulationConfig` (used by the batch system).
  It filters to valid dataclass fields, so adding a param means adding it to `config.py`.

### 1D core (the foundation everything reuses)
- `modelmain.py` — `Simulator`: ties config + grid + RTE + time-stepping together. `.run()` is the main
  loop; `_T_surf_calc` (surface energy balance Newton solve), `_bc*` (boundary conditions), and
  `_fd1d_heat_implicit_diag*` (the banded conduction solve) are the hot path. Also holds crater setup and
  most plotting/post-processing helpers. **This file is huge (~2000 lines) and TAB-indented** (see
  Conventions).
- `grid.py` — `LayerGrid`: the 1D vertical grid. Its coordinate is **optical depth τ**, and its banded
  conduction operator `diag` carries an `Et²` factor in the conductivity that exactly cancels the `Et²`
  in τ-spacing², so the τ operator is algebraically identical to a physical-metre operator. This identity
  is load-bearing for the 3D extension.
- `rte_disort.py` (`DisortRTESolver`) and `rte_hapke.py` (`RadiativeTransfer`) — the two interchangeable
  RTE cores, selected by `cfg.RTE_solver`. DISORT runs in three spectral modes: `two_wave` (broadband,
  fastest — used for thermal evolution), `multi_wave` (per-wavenumber, needs optical-constant input
  files), and `hybrid` (broadband visible + multi-wave thermal — used for computing output spectra).

### Radiance / spectral post-processing
- `radiance_processor.py` + `postprocessing/radiance_calculator.py` — compute emergent radiance and
  brightness-temperature spectra from saved thermal results. **Radiance always uses DISORT** regardless of
  which solver ran the thermal evolution. `observer_radiance.py` handles observer-geometry radiance.

### Roughness, craters, and real terrain (feature-branch work)
- `crater.py` — `CraterMesh` (hemispherical crater mesh + subdivision), `ShadowTester` (per-facet direct-beam
  illumination via ray casting), `SelfHeatingList` (**reader** of sparse view-factor files),
  `CraterRadiativeTransfer` (facet radiative coupling: shadowing, multiple scattering, self-heating).
  In crater/terrain mode `Simulator` batches DISORT over all facets (`n_cols=n_facets`) and advances all
  columns with one banded solve, so those scale linearly; the per-step term that grows superlinearly with
  facet count is the `ShadowTester` ray cast (O(N_sub²)), amortised by `cfg.illum_freq` and, since
  2026-09-29, compiled by `cfg.shadow_backend='auto'` (numba, ~60–90x on DEM meshes, identical results;
  `prototypes/test_shadow_backend.py`).
- `topography.py` — `DEMMesh(CraterMesh)`: turns a real DEM (GeoTIFF/ASCII/npy) into a mesh with the exact
  `CraterMesh` attribute set, so it drops into `ShadowTester`/`CraterRadiativeTransfer`/`view_factors`.
- `view_factors.py` — `ViewFactorList`: **computes** facet-to-facet view factors (with LOS occlusion) and
  writes them in `SelfHeatingList`'s sparse format, so the crater engine reads them unchanged.
- `terrain_bt.py` — `TerrainObserver`: per-facet observer-geometry brightness-temperature spectra for
  terrain runs.

### 3D conduction (feature/3d-conduction)
- `grid3d.py` (`VolumeGrid`), `sim3d.py` (`Simulator3D`), `radiance3d.py` — a structured [nx,ny,nz]
  conduction model built on the 1D core: vertical sweep reuses `LayerGrid.diag`; lateral (x,y) uses
  Neumann-wall second differences in τ units; one implicit step is LOD/ADI (three batched `solve_banded`
  sweeps). Validated to reduce to the 1D `Simulator` exactly. See `docs/README_3D.md` and
  `docs/3D_conduction_prototype_plan.md`.

Data flow: thermal results → HDF5; radiance is a **separate post-processing stage** off saved thermal data,
which is why the batch system can run thermal-only sweeps fast and post-process spectra later.

## Conventions & gotchas

- **Indentation is not uniform.** `modelmain.py` (and other original Ryan-authored files) use **tabs**;
  `config.py`, `grid3d.py`, `sim3d.py`, and other newer modules use **spaces**. Match the file you're
  editing — do not reformat.
- **Vertical coordinate is optical depth τ, not metres.** The `Et²` cancellation described above is
  intentional; preserve it when touching conduction operators.
- New config parameters must be added to `SimulationConfig` in `config.py` (with a doc comment) to be
  usable by both the direct and YAML/batch paths.

## HANDOFF.md, agenda.json, and the CC/CS collaboration

`HANDOFF.md` is an active message log between **CC** (Claude Code, builds this repo) and **CS** (Claude
Science, does analysis on the same tree). **Newest entries are appended at the TOP**; the header format is
`## YYYY-MM-DD — AUTHOR → RECIPIENT — subject`, and open blockers are tagged `[NEEDS DECISION]`. When the
human says "check HANDOFF," read the top entries — this is the narrative source of truth, and it is
git-versioned. Keep new entries short and link to files/commits.

- **Editing HANDOFF.md:** insert the new entry ABOVE the current top `## ` header — never replace that line
  (an earlier edit pattern silently ate headers). `tools/check_handoff.py` enforces the append-only header
  invariant and is installed as the pre-commit hook.
- **Live state is `agenda.json`, not HANDOFF flags.** Start a session with `python3 tools/agenda.py --for CC`
  (exit 1 = something is waiting on you); close with `--close A<n> --note "..."`, hand over with
  `--open --owner CS --from CC --title "..." [--blocking]`.
- `tools/check_science_gates.py [files.npz]` runs the physics acceptance gates on run outputs;
  `--pair dry_convergence.npz ice_convergence.npz` is the paired matched-season drift gate that decides
  whether an ice–dry differential is converged (needs |drift| < 0.1 K/yr).
- `tools/watch_handoff.py` is a cron-able digest writer (appends to `weekend_progress.md`).

## Production runs (Puma HPC) — the PSR campaign driver is NOT in this repo

`run_psr_floor_puma.py` lives on Puma in `/groups/sbyrne/phillipsm/psr_run/` (repo clone at
`/groups/sbyrne/phillipsm/ThermoSpec`, micromamba env `thermospec`); a stale copy sits in the CS sync tree
(`~/Documents/Research/Publications/artemis-thermal-modeling/claude_session_sync/scripts/`). It builds a
`DEMMesh` + `compute_view_factors(..., occlusion_backend='auto')` + SPICE sun vectors and injects them via
`Simulator(cfg, crater_mesh=, crater_selfheating=, sun_vectors=)`. Outputs (`*_thermal_*.npz`,
`*_convergence_*.npz`, `*_psr_floor_*.npz`) land in `psr_run/` and are rsynced to the sync tree's `data/`
for gating and analysis. SLURM via `ssh junonia`, transfers via `filexfer`; submit with a file-based
`sbatch job.slurm </dev/null` (`--wrap` over ssh hangs). Edit the driver by copying it down, patching,
and copying back with a dated `.bak-cc-*` on Puma.

## Git

`main` is the integration branch. `feature/3d-conduction` and `feature/hybrid-evolution` are fully merged
into `main`. Live work is on **`feature/terrain-viewfactors`** (>130 commits ahead of `main` as of
2026-09-29): the PSR/terrain campaign, `HANDOFF.md`, `agenda.json`, and `tools/` all live there. Branch
before starting unrelated new work; commit only when asked, but **push `feature/terrain-viewfactors` at the
end of every session** even with no new commits (agenda A5 contract — CS reads the tree via origin).
