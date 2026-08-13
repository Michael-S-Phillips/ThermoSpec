"""Simulator3D: a 3D-conduction diurnal thermal model built from ThermoSpec's own components.

Phase 3 of the 3D-conduction build (see docs/3D_conduction_prototype_plan.md). This is the
non-RTE ("traditional") thermal model made fully 3D in conduction: a structured [nx, ny] field
of vertical columns, each the existing 1D column, coupled laterally by conduction through
VolumeGrid's LOD-ADI step. Every per-step operation mirrors modelmain.Simulator's non-RTE path
exactly, but vectorized over the lateral grid:

  * initial state T = T_bottom everywhere, T_surf = T_bottom  (modelmain.py:74-76);
  * solar geometry mu(t), F(t) from the same formula (modelmain.py:328-331);
  * conduction advance via VolumeGrid.step (reuses LayerGrid.diag for the vertical);
  * surface energy balance solved per column by the same Newton iteration as
    modelmain._T_surf_calc, then the virtual top node and bottom BC set as in _bc_noRTE.

Per-column solar can be modulated with `mu_fac` (multiplies the incidence cosine, e.g. facet
tilt) and `F_gate` (0/1 shadow mask), enabling lateral illumination contrast. With both uniform
the run reduces exactly to the 1D Simulator, column-for-column.

RTE coupling (per-column DISORT/Hapke source injection) is a later phase; this module asserts
non-RTE.
"""
import numpy as np

from config import SimulationConfig
from grid import LayerGrid
from grid3d import VolumeGrid


class Simulator3D:
    def __init__(self, cfg: SimulationConfig, nx, ny, dx_m, dy_m, lateral_k=None):
        if cfg.use_RTE:
            raise NotImplementedError("Simulator3D currently supports non-RTE conduction only.")
        self.cfg = cfg
        self.base = LayerGrid(cfg)
        self.vol = VolumeGrid(self.base, nx, ny, dx_m, dy_m, lateral_k=lateral_k)
        self.nx, self.ny, self.nz = int(nx), int(ny), self.base.x_num

        self.T = np.full((self.nx, self.ny, self.nz), float(cfg.T_bottom))
        self.T_surf = np.full((self.nx, self.ny), float(cfg.T_bottom))

        # per-column solar modulation (defaults: flat, fully lit)
        self.mu_fac = np.ones((self.nx, self.ny))
        self.F_gate = np.ones((self.nx, self.ny))

        self._setup_solar()

        # surface-BC constants, identical to modelmain._T_surf_calc / _bc_noRTE
        g = self.base
        Et = cfg.Et
        Et1 = Et[1] if np.ndim(Et) else Et
        self._dx_surf = ((g.x[1] - g.x[0]) / 2.0) / Et1
        self._k_dx = g.cond[1] / Et1**2 / self._dx_surf
        self._x0, self._x1 = g.x[0], g.x[1]

    def _setup_solar(self):
        cfg = self.cfg
        P = cfg.P
        self.t_num = cfg.tsteps_day * cfg.ndays
        t = np.linspace(0, P * cfg.ndays, self.t_num)
        hour_angle = (np.pi / (P / 2)) * (t - (P / 2))
        mu = (np.sin(cfg.dec) * np.sin(cfg.latitude)
              + np.cos(cfg.latitude) * np.cos(hour_angle) * np.cos(cfg.dec))
        self.mu_array = mu
        self.F_array = (mu > 0.001).astype(float)

    def _surface_bc(self, mu, F):
        """Vectorized copy of modelmain._T_surf_calc + _bc_noRTE + bottom BC, per column."""
        cfg = self.cfg
        Q = (F * self.F_gate) * cfg.J * (mu * self.mu_fac) * (1.0 - cfg.albedo)
        se = cfg.sigma * cfg.em
        T1 = self.T[:, :, 1]
        Ts = self.T_surf
        for _ in range(cfg.T_surf_max_iter):
            W = Q + self._k_dx * (T1 - Ts) - se * Ts**4
            dWdT = -self._k_dx - 4.0 * se * Ts**3
            dT = W / dWdT
            if np.all(np.abs(dT) < cfg.T_surf_tol):
                break
            Ts = Ts - dT
        self.T_surf = Ts
        # virtual top node enforces the surface flux (modelmain._bc_noRTE)
        self.T[:, :, 0] = (T1 - Ts) * (self._x0 / self._x1) + Ts
        # bottom boundary
        if cfg.bottom_bc == 'dirichlet':
            self.T[:, :, -1] = cfg.T_bottom
        elif cfg.bottom_bc == 'neumann':
            self.T[:, :, -1] = self.T[:, :, -2]
        else:
            raise ValueError(f"Invalid bottom_bc: {cfg.bottom_bc}")

    def run(self, record_surf=False):
        surf_hist = []
        for j in range(self.t_num):
            if j > 0:
                self.T = self.vol.step(self.T, 0.0)          # non-RTE: no volumetric source
                self._surface_bc(self.mu_array[j], self.F_array[j])
            if record_surf:
                surf_hist.append(self.T_surf.copy())
        if record_surf:
            self.surf_hist = np.array(surf_hist)
        return self.T, self.T_surf
