"""VolumeGrid: structured 3D conduction on top of ThermoSpec's 1D LayerGrid.

Design (see docs/3D_conduction_prototype_plan.md):
  * The vertical (z) direction reuses the existing LayerGrid banded operator `diag` unchanged.
    ThermoSpec's vertical coordinate is optical depth tau, and its conduction operator carries
    an Et^2 factor in the conductivity that exactly cancels the Et^2 in tau-spacing^2 -- so the
    tau operator is algebraically identical to a physical-metre operator. We therefore drive the
    z-sweep straight from `base_grid.diag`.
  * The lateral (x, y) directions are built here as second-difference operators with insulating
    (Neumann, zero-flux) walls, in tau units for consistency with the vertical operator: the
    lateral spacing in metres is converted to tau via Et, and the lateral diffusivity is the
    grid's own diffusivity K (isotropic material) unless overridden.

One implicit step is LOD/ADI: sweep x, then y, then z, each a batched `solve_banded`. With
lateral conduction off (lateral_k=0) or a laterally-uniform field, the step reduces exactly to
the 1D solve column-by-column.

This handles conduction only. Surface/bottom boundary conditions and the (1D, per-column) RTE
source term are applied by the caller, exactly as in the 1D model.
"""
import numpy as np
from scipy.linalg import solve_banded


def _lateral_banded_neumann(n, r):
    """Banded (I - r D2) for a uniform 1D lateral line of n nodes with zero-flux end walls.

    r = dt K / h^2 (may be a scalar). scipy solve_banded((1,1)) layout: ab[0]=super, ab[1]=main,
    ab[2]=sub. Neumann walls (reflecting) conserve total heat and make the operator identity when
    r = 0.
    """
    ab = np.zeros((3, n))
    ab[1, :] = 1.0 + 2.0 * r
    ab[1, 0] = 1.0 + r          # zero-flux at the two ends: one neighbour only
    ab[1, -1] = 1.0 + r
    ab[0, 1:] = -r              # super-diagonal
    ab[2, :-1] = -r             # sub-diagonal
    return ab


class VolumeGrid:
    """LOD-ADI 3D conduction step reusing a LayerGrid's vertical operator.

    Parameters
    ----------
    base_grid : LayerGrid
        Fully built 1D grid; supplies the vertical operator (`diag`), diffusivity `K`, time step
        `dt`, node count `x_num`, and config (for Et).
    nx, ny : int
        Number of lateral columns in x and y.
    dx_m, dy_m : float
        Lateral node spacing in metres (uniform).
    lateral_k : None | float
        None  -> isotropic: lateral diffusivity equals the vertical grid diffusivity K(z).
        0.0   -> no lateral conduction (operators become identity; step reduces to 1D).
        other -> anisotropy factor multiplying the vertical diffusivity (lateral K = factor*K).
    """

    def __init__(self, base_grid, nx, ny, dx_m, dy_m, lateral_k=None):
        self.g = base_grid
        self.nx, self.ny = int(nx), int(ny)
        self.nz = base_grid.x_num
        self.dt = base_grid.dt
        self.diag = base_grid.diag

        cfg = base_grid.config
        Et = cfg.Et
        Et_arr = Et if np.ndim(Et) else np.full(self.nz, float(Et))
        Kz = np.asarray(base_grid.K, dtype=float)          # tau diffusivity per depth node

        if lateral_k == 0.0:
            rx = np.zeros(self.nz)
            ry = np.zeros(self.nz)
        else:
            fac = 1.0 if lateral_k is None else float(lateral_k)
            dx_tau = dx_m * Et_arr                          # lateral spacing in tau units
            dy_tau = dy_m * Et_arr
            rx = self.dt * fac * Kz / dx_tau**2
            ry = self.dt * fac * Kz / dy_tau**2

        self._abx = [_lateral_banded_neumann(self.nx, float(rx[k])) for k in range(self.nz)]
        self._aby = [_lateral_banded_neumann(self.ny, float(ry[k])) for k in range(self.nz)]
        self._lateral_off = bool(np.all(rx == 0.0) and np.all(ry == 0.0))

    def step(self, T, source=0.0):
        """Advance one implicit time step. T: [nx, ny, nz]. source: scalar or [nx, ny, nz]."""
        b = T + self.dt * source
        if not self._lateral_off:
            for k in range(self.nz):                        # x-sweep: solve along axis 0
                b[:, :, k] = solve_banded((1, 1), self._abx[k], b[:, :, k])
            for k in range(self.nz):                        # y-sweep: solve along axis 1
                b[:, :, k] = solve_banded((1, 1), self._aby[k], b[:, :, k].T).T
        # z-sweep: reuse the 1D vertical operator for every column at once
        B = b.reshape(self.nx * self.ny, self.nz).T         # [nz, ncols]
        B = solve_banded((1, 1), self.diag, B)
        return np.ascontiguousarray(B.T).reshape(self.nx, self.ny, self.nz)
