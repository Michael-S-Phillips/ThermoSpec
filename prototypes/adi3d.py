"""Prototype 3D conduction solver by locally-one-dimensional (LOD) ADI.

One implicit step of dT/dt = K (T_xx + T_yy + T_zz) + S is factored into three sequential 1D
implicit sweeps, one per axis:

    (I - dt K D2_x) T1 = T^n + dt S
    (I - dt K D2_y) T2 = T1
    (I - dt K D2_z) T^{n+1} = T2

Each sweep is a batched tridiagonal solve (scipy.linalg.solve_banded), the exact primitive
ThermoSpec's 1D core already uses. The scheme is unconditionally stable and first order in time
(the O(dt) splitting error is negligible at the ~25 s diurnal step); each 1D operator is the
standard second difference with Dirichlet identity rows at the two boundary planes, so a caller
that leaves boundary planes fixed gets Dirichlet walls for free.

This prototype uses uniform spacing and constant properties on purpose: it isolates ADI
correctness and wall-time from the tau-unit / non-uniform / RTE complications that only appear at
integration. The eigenmode test in test_adi3d.py pins it to machine precision.
"""
import numpy as np
from scipy.linalg import solve_banded


def _backward_euler_banded(n, r):
    """Banded form of (I - r D2) for a uniform 1D grid of n nodes, Dirichlet identity end rows.

    r = dt K / h^2. scipy solve_banded((1,1), ab) layout: ab[0]=super, ab[1]=main, ab[2]=sub.
    """
    ab = np.zeros((3, n))
    ab[1, :] = 1.0 + 2.0 * r     # main diagonal (interior)
    ab[0, 2:] = -r               # super-diagonal for interior rows 1..n-2
    ab[2, :-2] = -r              # sub-diagonal for interior rows 1..n-2
    # Dirichlet identity rows at the two boundary planes
    ab[1, 0] = 1.0
    ab[1, -1] = 1.0
    ab[0, 1] = 0.0               # no coupling out of row 0
    ab[2, -2] = 0.0              # no coupling out of row n-1
    return ab


class ADI3DSolver:
    """LOD-ADI implicit conduction step on a uniform [nx,ny,nz] grid, constant K."""

    def __init__(self, shape, spacing, K, dt):
        self.shape = tuple(int(s) for s in shape)
        self.dt = float(dt)
        (nx, ny, nz) = self.shape
        (hx, hy, hz) = spacing
        self._ab = [
            _backward_euler_banded(nx, self.dt * K / hx**2),
            _backward_euler_banded(ny, self.dt * K / hy**2),
            _backward_euler_banded(nz, self.dt * K / hz**2),
        ]

    def _sweep(self, T, axis):
        """Implicit solve along `axis` for every line, via one batched solve_banded."""
        Tm = np.moveaxis(T, axis, 0)
        n = Tm.shape[0]
        B = np.ascontiguousarray(Tm.reshape(n, -1))
        X = solve_banded((1, 1), self._ab[axis], B)
        return np.moveaxis(X.reshape(Tm.shape), 0, axis)

    def step(self, T, source=0.0):
        """Advance one time step. T: array [nx,ny,nz]. source: scalar or same-shaped array."""
        b = T + self.dt * source
        for axis in (0, 1, 2):
            b = self._sweep(b, axis)
        return b
