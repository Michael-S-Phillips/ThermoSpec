# 3D conduction in ThermoSpec

A structured-grid, 3D-conduction thermal model built on ThermoSpec's existing 1D core. Heat
conducts in all three directions; each surface column couples to radiation through the existing
1D-per-column machinery (radiation in regolith is effectively vertical). Status: the thermal
model is **fully 3D in conduction** and validated to reduce to the 1D `Simulator` exactly, in the
**non-RTE**, **DISORT**, and **Hapke** radiative-transfer paths (`two_wave`, single-layer), with
**per-column emissivity/brightness-temperature spectra**. Remaining variants (multi_wave/hybrid
thermal *evolution*, two-layer RTE) are follow-ons — see "Not yet" below.

## Modules

| file | role |
|---|---|
| `grid3d.py` | `VolumeGrid`: one LOD-ADI implicit conduction step over an [nx,ny,nz] field. Reuses `LayerGrid.diag` for the vertical sweep; builds conservative Neumann-wall lateral operators. Supports per-column vertical operators for temperature-dependent properties. |
| `sim3d.py` | `Simulator3D`: diurnal driver. Non-RTE path (per-column nonlinear surface energy balance, vectorized copy of `modelmain._T_surf_calc`) and RTE path (batched DISORT over all columns, `n_cols=nx*ny`; Neumann surface BC). Bottom BC + temperature-dependent property updates in both. `run(record_phases=True)` + `phase_spectra()` for noon/pre-dawn spectra. |
| `radiance3d.py` | `compute_spectra`: batched per-column hybrid-thermal DISORT output solve -> emergent thermal spectrum + closed-form per-band brightness temperature. The 3D analogue of `radiance_processor` (hybrid, thermal_only). |
| `prototypes/adi3d.py` | standalone LOD-ADI reference solver (uniform grid, constant props) + its machine-precision eigenmode tests and the wall-time benchmark. |
| `prototypes/test_*.py` | test suites (run directly; no pytest needed). |
| `docs/RATE_LIMITING_AND_PINN.md` | measured cost breakdown + where surrogate models pay off. |

## How it works

- **Vertical (z):** ThermoSpec's coordinate is optical depth tau, and its conduction operator's
  conductivity carries an Et^2 that exactly cancels the Et^2 in tau-spacing^2 -- so the tau operator
  equals a physical-metre operator. The 3D z-sweep therefore reuses `LayerGrid.diag` unchanged.
- **Lateral (x,y):** conservative second differences with insulating (Neumann) walls, in tau units,
  isotropic diffusivity by default. Virtual top/bottom ghost nodes are excluded from lateral
  conduction (they hold BC-enforcement values, not physical temperatures).
- **Time step:** LOD/ADI -- three sequential 1D implicit `solve_banded` sweeps. Unconditionally
  stable, first order in time (the ~25 s diurnal step makes the O(dt) splitting error negligible).
- **Reduction to 1D:** with lateral conduction off, or a laterally-uniform field, the step reproduces
  the 1D banded solve column-for-column.

## Usage

```python
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"; os.environ["OMP_NUM_THREADS"] = "1"
from config import SimulationConfig
from sim3d import Simulator3D

cfg = SimulationConfig(
    use_RTE=False, single_layer=True, diurnal=True, sun=True,
    depth_dependent_properties=True, temperature_dependent_properties=True, temp_dependent_cp=True,
    geometric_spacing=True, auto_dt=False, tsteps_day=100000, ndays=8,
    dust_thickness=0.5, T_bottom=250.0, bottom_bc='dirichlet',
    latitude=0.0, dec=0.0, P=29.5306*24*3600, S=1361.0,
)
sim = Simulator3D(cfg, nx=32, ny=32, dx_m=0.05, dy_m=0.05, lateral_k=None)  # None=isotropic
# lateral illumination structure (optional): shadow mask and per-column incidence factor
sim.F_gate[:, :] = 1.0        # 0 shadows a column; sim.mu_fac scales its incidence cosine
sim.run()                     # sim.T is [nx,ny,nz]; sim.T_surf is [nx,ny]
```

`lateral_k`: `None` isotropic (lateral k = vertical k), `0.0` off (columns independent), or a float
anisotropy factor. `F_gate` (0/1 shadow mask) and `mu_fac` (incidence-cosine multiplier) impose
lateral illumination contrast -- the driver of lateral heat flow.

## Validation (all in `prototypes/test_*.py`, all green)

- ADI reproduces the exact discrete sine-eigenmode decay to 1e-10; unconditionally stable at huge dt.
- Lateral operator conserves heat (Neumann) and relaxes a contrast to uniform.
- A laterally-uniform 3D run reproduces the real 1D `Simulator` to <1e-6 for uniform, geometric, and
  depth-dependent grids, and to <1e-5 with temperature-dependent cp (mid-run operator rebuilds).
- Lateral conduction warms a permanently shadowed half and cools the lit half (physical transport).

## RTE usage (DISORT)

```python
cfg = SimulationConfig(
    use_RTE=True, RTE_solver='disort', thermal_evolution_mode='two_wave',
    single_layer=True, diurnal=True, sun=True, geometric_spacing=True,
    auto_dt=False, tsteps_day=100000, ndays=8, dust_thickness=0.5, T_bottom=250.0,
    latitude=0.0, dec=0.0, P=29.5306*24*3600, S=1361.0, nstr=4, nmom=4,
    ssalb_therm=0.1, ssalb_vis=0.5, eta=1.0,
)
sim = Simulator3D(cfg, nx=8, ny=8, dx_m=0.05, dy_m=0.05, lateral_k=None)
sim.F_gate[...] = ...   # optional shadow mask / sim.mu_fac for facet tilt
sim.run()               # per step: one batched DISORT thermal + one visible solve over all columns
```
Set `RTE_solver='hapke'` instead for the fast broadband core (pure numpy, looped per column, no
optics/torch) — ideal for large-grid scans; `'disort'` for the accurate N-stream solve. Use a dt
that resolves the near-surface layer (dt <~ 20-90 s; Hapke's BVP feedback wants the finer end); at
coarse dt the model's own near-surface response is under-resolved. Per-column RTE dominates cost —
see `docs/RATE_LIMITING_AND_PINN.md`.

## Not yet (next phases)

- **multi_wave thermal evolution** is blocked by the absent solar-spectrum files (same reason as
  multi_wave *output*); use `hybrid` or `two_wave`. **Hybrid** thermal evolution IS supported
  (`thermal_evolution_mode='hybrid'`) but is ~850x two_wave per column (one DISORT per band) --
  practical only for small grids / short runs, or as ground truth for a learned per-column source
  surrogate (see `docs/RATE_LIMITING_AND_PINN.md`). two_wave evolution + hybrid *output*
  (`radiance3d.py`) remains the practical path for production.
- **Native integration.** Currently a companion module reusing the core; a future step could fold a
  `conduction_3d` mode into `modelmain.Simulator` directly.
- **Two-layer RTE** (`single_layer=False`) needs the rock/dust interface source term handled in the
  3D step (`modelmain._fd1d_heat_implicit_diag` has special interface-node handling with no analogue
  here yet).

## Running the tests

```bash
cd <repo root>
for t in adi3d grid3d sim3d; do
  KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 python prototypes/test_$t.py
done
```
