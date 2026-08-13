# 3D conduction in ThermoSpec

A structured-grid, 3D-conduction thermal model built on ThermoSpec's existing 1D core. Heat
conducts in all three directions; each surface column still couples to radiation through the
existing 1D-per-column machinery (radiation in regolith is effectively vertical). Status: the
**non-RTE (traditional) thermal model is fully 3D in conduction** and validated to reduce to the
1D `Simulator` exactly. RTE coupling is the next phase (see "Not yet" below).

## Modules

| file | role |
|---|---|
| `grid3d.py` | `VolumeGrid`: one LOD-ADI implicit conduction step over an [nx,ny,nz] field. Reuses `LayerGrid.diag` for the vertical sweep; builds conservative Neumann-wall lateral operators. Supports per-column vertical operators for temperature-dependent properties. |
| `sim3d.py` | `Simulator3D`: non-RTE diurnal driver. Per-column nonlinear surface energy balance (vectorized copy of `modelmain._T_surf_calc`), bottom BC, temperature-dependent property updates. |
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

## Not yet (next phases)

- **RTE per column.** Wire the existing 1D DISORT/Hapke solve into each surface column and inject its
  flux-divergence source into that column's top cells (`rte_disort.py`/`rte_hapke.py` unchanged). This
  is the expensive part -- see `docs/RATE_LIMITING_AND_PINN.md`; it is also where a learned surrogate
  belongs.
- **Lateral-sweep performance.** The per-depth `solve_banded` loop is ~90% of the conduction step and
  should be replaced by a vectorized/batched tridiagonal (or a torch-MPS batched solve).
- **Native integration.** Currently a companion module reusing the core; a future step could fold a
  `conduction_3d` mode into `modelmain.Simulator` directly.
- **Temperature-dependent k** is supported in the operator but the surface-BC `k_dx` uses the static
  conductivity; refine if `temp_dependent_k=True` is needed (the lunar runs use `temp_dependent_k=False`).

## Running the tests

```bash
cd <repo root>
for t in adi3d grid3d sim3d; do
  KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 python prototypes/test_$t.py
done
```
