# Rate-limiting steps in 3D ThermoSpec, and where PINN surrogates would pay off

Measured on this machine (CPU, single thread, `OMP_NUM_THREADS=1`). Grid: single-layer lunar
column, geometric spacing, depth-dependent properties, nz = 83 vertical nodes; lateral grid as
noted. "step" = one implicit time step.

## Measured cost breakdown (non-RTE 3D conduction)

| component | 16x16x83 | 32x32x83 | share | notes |
|---|---|---|---|---|
| **lateral ADI sweeps (x+y)** | 1.76 ms | 2.40 ms | **~90%** | per-depth Python loop: nz x 2 = 166 `solve_banded` calls/step |
| z-sweep (vertical) | 0.10 ms | 0.32 ms | ~11% | one batched `solve_banded`, reuses `LayerGrid.diag` |
| surface BC (Newton) | negligible | negligible | <1% | vectorized over columns |
| full step | 1.86 ms | 2.71 ms | | 538 / 369 steps/s |
| temp-dependent operator rebuild | 12.6 us / column | | episodic | only for columns whose T drifts > threshold |

Projected wall time (800k-step, 8-lunation run): ~36 min at 32x32x83, non-RTE.

## The rate-limiting steps, ranked

1. **Per-column RTE solve (when `use_RTE=True`) — dominant by orders of magnitude.**
   RTE is 1D per column, so a 3D run calls DISORT/Hapke once per surface column per step. From the
   earlier 1D benchmark: a full 800k-step run costs ~75 s per column with the Hapke core and ~30 min
   per column with DISORT at nstr=16. Times a lateral grid:
   - 32x32 = 1024 columns x Hapke ~75 s  ->  ~21 hours
   - 1024 columns x DISORT nstr=16 ~30 min  ->  ~500+ hours (intractable)
   Conduction (36 min) is a rounding error next to this. **This is the single highest-value
   surrogate target in the whole model.**

2. **Lateral ADI sweeps (dominant *within* non-RTE conduction).** 90% of the conduction step, but
   this is a software/linear-algebra cost, not a physics cost: it is a Python loop of `solve_banded`
   over depth levels. It should be *optimized*, not surrogated (see below).

3. **Per-column vertical operator rebuilds (temperature-dependent properties).** 12.6 us/column when
   a column crosses `temp_change_threshold`. Cheap per event but scales with column count and rebuild
   frequency; episodic.

## What is a PINN/NN-surrogate target and what is not

**Surrogate the RTE, not the conduction.** The distinction is repeated-expensive-nonlinear-map vs.
cheap-linear-algebra:

- **RTE per column = ideal surrogate.** It is an expensive, smooth, repeatedly-evaluated map with a
  low-dimensional interface: `(vertical temperature profile, fixed optical properties, solar mu*F)
  -> (radiative flux-divergence source term over the column, upward flux)`. The optics are fixed
  within a run, illumination is one scalar, and the temperature profile lives on ~80 nodes (or a few
  POD/PCA coefficients). A network that maps profile -> source term + upward flux, trained on DISORT
  outputs, would replace the per-step per-column DISORT call and is exactly what makes **3D + DISORT
  spectra tractable**. This is the natural home for the `refactor/pinn-surrogate-infrastructure`
  work. Physics-informed constraints that fit naturally: energy conservation (integral of the source
  equals net radiative flux divergence), non-negativity of upward flux, monotone response to a
  uniform temperature shift.
  - Two granularities worth prototyping: (a) surrogate the **broadband two_wave source term** used in
    thermal evolution (drop-in for `rte_disort.disort_run` / `rte_hapke.compute_source`); (b)
    surrogate the **hybrid/multi-wave emergent spectrum** used only at output times (fewer calls, but
    the expensive nstr=16 path). (a) removes the per-step cost; (b) removes the per-output cost.

- **Conduction = do NOT surrogate; optimize.** The ADI sweeps are exact, cheap, and unconditionally
  stable linear solves. A surrogate would trade guaranteed accuracy/stability for nothing (they are
  already fast). The lateral-sweep bottleneck is a *loop*, not physics: replace the per-depth
  `solve_banded` loop with (i) a vectorized batched tridiagonal (Thomas) solve over all depth levels
  at once, or (ii) a torch batched tridiagonal solve on MPS/GPU (torch is already a dependency). Both
  keep exactness and should cut the 90% lateral cost several-fold.

- **Surface BC Newton, operator rebuilds = leave alone.** Already negligible / episodic.

## One-line guidance for the PINN branch

Point the surrogate effort at `rte_disort.disort_run` (and `rte_hapke.compute_source`): learn
`T_profile (+ mu*F) -> (source_term, flux_up)` from DISORT training data, with energy-conservation as
the physics constraint. That single substitution is what turns 3D-with-radiative-transfer from
"hundreds of hours" into "feasible," whereas the conduction solver just needs a better tridiagonal
kernel.
