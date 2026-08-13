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

from stencils import fd1d_heat_implicit_diagonal_nonuniform_kieffer


def build_vertical_diag(lthick, dens, heat, cond, dt):
    """Banded vertical operator for one column, from ThermoSpec's Kieffer coefficients.

    Identical formula to grid.LayerGrid._build_fd_matrix / _update_fd_matrix, exposed so a 3D
    run with temperature-dependent properties can rebuild a per-column operator from that
    column's current cp(T)/k(T). All arrays are length nz (with virtual end nodes).
    """
    A1 = (2 * dt * cond[1:-1]
          / (dens[1:-1] * heat[1:-1] * lthick[1:-1]**2
             * (1 + (lthick[0:-2] * cond[1:-1]) / (lthick[1:-1] * cond[0:-2]))))
    A3 = ((1 + (lthick[0:-2] * cond[1:-1]) / (lthick[1:-1] * cond[0:-2]))
          / (1 + (lthick[2:] * cond[1:-1]) / (lthick[1:-1] * cond[2:])))
    A2 = -1.0 * (1.0 + A3)
    return fd1d_heat_implicit_diagonal_nonuniform_kieffer(len(lthick), A1, A2, A3)


def _lateral_banded_neumann(n, r):
    """Banded (I - r D2) for a uniform 1D lateral line of n nodes with zero-flux end walls.

    r = dt K / h^2 (may be a scalar). scipy solve_banded((1,1)) layout: ab[0]=super, ab[1]=main,
    ab[2]=sub. Neumann walls (reflecting) conserve total heat and make the operator identity when
    r = 0.
    """
    ab = np.zeros((3, n))
    if n == 1:
        ab[1, 0] = 1.0          # a single isolated node has no neighbours: identity
        return ab
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
        self._Et_arr = Et if np.ndim(Et) else np.full(self.nz, float(Et))
        self._dx_m, self._dy_m = float(dx_m), float(dy_m)
        self._lateral_k = lateral_k
        self._build_lateral(np.asarray(base_grid.K, dtype=float))

        # Per-column vertical operators (temperature-dependent case). None -> shared self.diag.
        self._diag_cols = None

    def _build_lateral(self, Kz):
        """(Re)build the per-depth lateral operators from a tau-diffusivity profile Kz[nz]."""
        if self._lateral_k == 0.0:
            rx = np.zeros(self.nz)
            ry = np.zeros(self.nz)
        else:
            fac = 1.0 if self._lateral_k is None else float(self._lateral_k)
            dx_tau = self._dx_m * self._Et_arr              # lateral spacing in tau units
            dy_tau = self._dy_m * self._Et_arr
            rx = self.dt * fac * Kz / dx_tau**2
            ry = self.dt * fac * Kz / dy_tau**2
            # The top (index 0) and bottom (index -1) depth nodes are boundary/virtual ghosts
            # re-set by the caller's BC every step -- the top ghost holds 2*T_surf - T1, not a
            # physical temperature. They must NOT participate in lateral conduction, or the
            # lateral sweep diffuses non-physical ghost values between columns. Only real
            # subsurface nodes (1..nz-2) conduct laterally.
            rx[0] = rx[-1] = 0.0
            ry[0] = ry[-1] = 0.0
        self._abx = [_lateral_banded_neumann(self.nx, float(rx[k])) for k in range(self.nz)]
        self._aby = [_lateral_banded_neumann(self.ny, float(ry[k])) for k in range(self.nz)]
        self._lateral_off = bool(np.all(rx == 0.0) and np.all(ry == 0.0))

    def set_vertical_diag_cols(self, diag_cols):
        """Supply per-column vertical operators [ncols][3, nz] (temperature-dependent case),
        or None to revert to the shared 1D operator. Column order is C-order over (nx, ny)."""
        self._diag_cols = diag_cols

    def set_lateral_diffusivity(self, Kz):
        """Rebuild lateral operators from an updated (laterally-representative) diffusivity Kz."""
        self._build_lateral(np.asarray(Kz, dtype=float))

    def step(self, T, source=0.0):
        """Advance one implicit time step. T: [nx, ny, nz]. source: scalar or [nx, ny, nz].

        The volumetric source (e.g. the RTE flux divergence) is injected into the vertical
        solve, so that for a laterally-uniform field the step reduces exactly to the 1D update
        solve_banded(diag, T + dt*source). Lateral sweeps act on T alone first.
        """
        b = np.array(T, dtype=float)                        # copy; never mutate the caller's T
        if not self._lateral_off:
            for k in range(self.nz):                        # x-sweep: solve along axis 0
                b[:, :, k] = solve_banded((1, 1), self._abx[k], b[:, :, k])
            for k in range(self.nz):                        # y-sweep: solve along axis 1
                b[:, :, k] = solve_banded((1, 1), self._aby[k], b[:, :, k].T).T
        b = b + self.dt * source                            # source enters the vertical solve
        # z-sweep
        flat = b.reshape(self.nx * self.ny, self.nz)        # [ncols, nz], row c = column (i,j)
        if self._diag_cols is None:
            X = solve_banded((1, 1), self.diag, flat.T).T   # shared operator, one batched solve
        else:
            X = np.empty_like(flat)
            for c in range(flat.shape[0]):                  # per-column temperature-dependent op
                X[c] = solve_banded((1, 1), self._diag_cols[c], flat[c])
        return np.ascontiguousarray(X).reshape(self.nx, self.ny, self.nz)
