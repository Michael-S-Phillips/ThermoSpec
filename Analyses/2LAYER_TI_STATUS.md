# Two-layer effective-TI analysis (v2) — status / handoff

**Last updated**: 2026-07-11, end of session. If you are a fresh Claude
session picking this up: read this whole file before doing anything, then
check `git log --oneline main..HEAD` on this branch for the exact commit
history, and read `Analyses/2layer_effective_TI_analysis_v2.py` itself
(well-commented) before making changes.

## Where things live

- **Branch**: `worktree-we-re-working-on-a-glistening-popcorn`, pushed to
  `origin` (repo moved: `https://github.com/andrewryan2727/ThermoSpec.git`,
  old remote name `ThermoRT` still works as an alias). **Not yet merged to
  `main`** — the user needs to do that (or ask a future session to open a
  PR) when ready. All work described below is on this one branch, 4 commits
  on top of the user's own `617dd3b` (`main`):
  1. `9dfa7e8` — new analysis script + `grid.py` geometric two-layer dust
     spacing
  2. `e527658` — rock ghost-node bug fix + root-caused/fixed a major TI
     anomaly (see below)
  3. `a850543` — Mie downsampling + `dust_rte_max_lthick` override params
  4. `7d3ed8f` — `torch.no_grad()` fix + hard memory-safety cap
- **Main script**: `Analyses/2layer_effective_TI_analysis_v2.py`
- **Full design rationale**: `.claude/plans/we-re-working-on-a-glistening-popcorn.md`
  (this session's plan file — may or may not be visible to you depending on
  how you were invoked; this status doc is the fallback if it isn't).
- **Depends on** (pre-existing, from the user's own earlier sessions,
  committed alongside the new script since they were untracked before):
  `Analyses/hybrid_calibration.py`, `Analyses/spectral_calibration.py`, and
  `Optical_props/serpentine_mie_200wns_{1,5}um.txt` + `wn_bounds_200.txt`.

## What this is

A rewrite of `Analyses/2layer_effective_TI_analysis.py` (old script), fitting
an "effective single-layer thermal inertia" (TI) to two-layer dust-over-rock
thermal models, publication-track. Two two-layer models, both fit against
the same one-layer non-RTE lookup table (LUT):

1. **Biele branch** (`build_biele_two_layer_sweep`) — literal reproduction
   of Biele et al. 2019 (opaque dust, pure conduction, fixed material
   properties: `k_dust=0.0025`, `rho_dust=366`, `cp_dust=700`). Fit target:
   `T_surf`.
2. **RTE branch** (`build_rte_two_layer_sweep`) — full spectral
   (hybrid/multi-wave DISORT) treatment, phonon-only `k_dust` from the
   Gundlach & Blum contact-conductivity model as a function of grain
   diameter/porosity. Fit target: `Tb_bol` (directional bolometric
   brightness temperature).

Fitting the RTE model against a non-RTE LUT is **intentional**: mimics how
historical Bennu TI fits (non-RTE) would naively interpret data from the
"real" (RTE) physics. All 6 of Biele's apparent-TI extraction methods are
implemented (`fit_all_methods`): full/day/night chi-squared (headline) +
phase-lag/Tmax/Tmin (diagnostics).

## Key finding this session: a major bug, found and fixed

Initial full-resolution runs showed a wildly anomalous RTE-branch TI at 2cm
dust thickness (~728, vs. the dust's own true TI of ~48-70) — the user was
confident (correctly) this was a bug, not physics. Root cause, found via a
single-layer "infinite dust" diagnostic (ruled out two-layer geometry) and
checking what `auto_dt` would choose: **the RTE branch's near-surface grid
layer (`dust_rte_max_lthick`, in tau/optical-depth units) becomes extremely
thin in physical terms once real Mie-calibrated `Et` is used (~130,714/m,
~18x larger than earlier, non-physically-calibrated validation runs used) —
and the inherited fixed timestep count silently under-resolved it.**

**Fix implemented**: `run_with_convergence_check()` in the script —
independently, per-run, doubles `tsteps_day` until the fit-relevant output
(`T_surf`/`Tb_bol`) stops changing by more than a tolerance (default 0.25 K,
relaxed from an initial 0.1 K per the user for speed), capped at
`max_doublings` (default 8) to avoid unbounded compute. **No cross-run
reuse/interpolation** — deliberately kept simple per the user's explicit
preference (avoiding unverified assumptions about how a converged
resolution transfers across thickness/material combos), at real compute
cost. Convergence diagnostics (`converged`, `n_doublings`, `tsteps_day`,
`max_abs_diff`) are persisted per-run in the HDF5 output, never silently
trusted.

**Validated**: after the fix, the 2cm RTE case's TI dropped from ~728 to
~82 (true value ~70) — confirms the diagnosis was right, even though that
particular run didn't fully converge within the memory-safety cap (see
below). The Biele branch converges cleanly and exactly matches the true
pure-dust TI at thick dust, a strong independent correctness check.

## Second bug found: OOM crash, partially understood

A convergence-checked run got killed by macOS (`jetsam`, 74 GB) once
`tsteps_day` grew into the millions. Investigated: `rte_disort.py`'s
`disort_run` (called once per timestep, sometimes millions of times) had no
`torch.no_grad()` guard, and `Disort` inherits from
`torch.nn.cpp.ModuleWrapper` — every forward call was building an unneeded
autograd graph. **Applied `@torch.no_grad()` to `disort_run`** (safe: this
pipeline never trains/backprops anywhere — calibration uses
`scipy.optimize.brentq`, not gradients).

**This fix was NOT sufficient by itself** — a controlled 1M-total-step test
with the fix applied still peaked at 14.24 GB (a plain per-step history-list
estimate predicts only ~2 GB). The real per-step memory growth rate
(~14 KB/step) is barely different from the original crash's rate
(~11.6 KB/step). **Root cause of the remaining growth is NOT identified** —
candidates include PyTorch allocator caching behavior, or something else in
the per-step DISORT call path not addressed by `no_grad()` alone. This is a
real open item if you want tighter convergence on very-thick-dust RTE cases.

**Mitigation in place**: `run_with_convergence_check` has a hard
`max_tsteps_day` ceiling (default 200,000) — doubling stops (flagged
`hit_memory_cap: True`, `converged: False`) rather than risking another OOM,
regardless of `max_doublings`. This is a safety net, not a fix for the
underlying memory-scaling mystery. The user confirmed (2026-07-11) they're
satisfied with this level of resolution ("good on the memory issue") and do
not currently want further digging into the root cause.

## Current default settings (as of last commit)

Set in `__main__` and as function defaults — these are all real, deliberate
tradeoffs, not arbitrary:
- `convergence_tol_K=0.25` (relaxed from 0.1 for speed, per user)
- `dust_rte_max_lthick=0.10` (nudged up from the `SimulationConfig` default
  of 0.05, per user — modest, deliberately kept in tau/optical-depth units,
  not physical length, since the user recalled needing fine tau-unit
  resolution for physical realism in thin-coating scenarios in prior work
  and didn't want to risk that)
- `n_bands=10` (downsampled from the full ~200-band Mie file via
  `hybrid_calibration.downsample_mie_file` — uniform/unweighted averaging,
  a real accuracy tradeoff. The user flagged that a **smarter,
  property-change-weighted downsampling** would be better but is fine
  deferring that; could be done by visual inspection of the Mie spectrum
  too)
- `max_doublings=8`, `max_tsteps_day=200000` (memory safety cap)
- Real Gundlach & Blum constants: `GB_POISSONS_DEFAULT=0.269`,
  `GB_YOUNGS_DEFAULT=5.625e9` Pa (module-level constants in the script,
  user-supplied, no longer placeholders). `ks=1.0`, `X=0.5` remain
  provisional pending final values from the user.

## What's next (per the user, this session's last message)

**The user wants to pause on long/production-scale runs** and instead run
smaller, cheap test cases to independently validate three things before
trusting a big run:
1. **The fitting routines** (`chi2_fit`, `fit_all_methods`, the six
   apparent-TI methods) — are they finding sensible optima, not just
   "not obviously broken"?
2. **Parameter settings** — do the current defaults (listed above) actually
   produce trustworthy results at small scale, or do they need further
   tuning?
3. **The time-doubling convergence check itself** — is it behaving
   correctly and efficiently across a range of cases (not just the 2 thin/
   thick smoke-test points already exercised)?

No specific test plan has been designed yet for this validation pass — that
is the immediate next task. Good building blocks already exist in the
script: `build_single_layer_lut`, `build_biele_two_layer_sweep`,
`build_rte_two_layer_sweep`, `fit_all_methods`, `plot_summary` /
`_plot_thickness_grid` (diagnostic PDF with per-thickness curve overlays +
relative-RMS annotations), and `run_with_convergence_check` (now printing
per-doubling diagnostics). The most recent full run's output (before being
killed) is at `Analyses/output/` in this worktree — inspect what's there
before deciding what new small tests to run, since a 2-point + partial
12-point sweep's worth of HDF5/plots may already answer some questions.

## Not yet done (known, deferred, not blocking)

- Additional grain sizes beyond 1 µm / 5 µm — needs new Mie tables via the
  `Preprocessing/` pipeline (out of scope for this script).
- Multiple porosities — framework supports it (`porosities` param) but
  never exercised beyond the single `fill_frac=0.37` baseline.
- Smarter (non-uniform) Mie downsampling.
- The full production sweep (all 12 default thicknesses, `DUST_THICKNESS_DEFAULT`)
  was started once and killed partway through by the user for this
  smaller-scale-validation-first reason — not because it was failing.
