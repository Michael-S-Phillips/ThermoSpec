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
from grid3d import VolumeGrid, build_vertical_diag


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
        self._Et1 = Et[1] if np.ndim(Et) else Et
        self._dx_surf = ((g.x[1] - g.x[0]) / 2.0) / self._Et1
        # k_dx is per column so it can track a column's node-1 conductivity when
        # temp_dependent_k rebuilds it (modelmain._T_surf_calc reads grid.cond[1] live).
        self._k_dx = np.full((self.nx, self.ny), g.cond[1] / self._Et1**2 / self._dx_surf)
        self._x0, self._x1 = g.x[0], g.x[1]

        self._setup_temperature_dependence()

    def _setup_temperature_dependence(self):
        """Prepare per-column vertical operators when properties depend on temperature.

        Mirrors grid.check_and_update_temperature_dependent_properties, but per column: each
        column rebuilds its own vertical operator from its own cp(T)/k(T) once its temperature
        has drifted past temp_change_threshold since its last rebuild. For a laterally-uniform
        field every column tracks the 1D operator in lockstep.
        """
        cfg = self.cfg
        g = self.base
        self._temp_dep = (cfg.temperature_dependent_properties
                          and (cfg.temp_dependent_cp or cfg.temp_dependent_k))
        if not self._temp_dep:
            return
        ncols = self.nx * self.ny
        self._lthick = g.l_thick
        self._dens = g.dens
        self._cond_base = np.asarray(g.cond, dtype=float).copy()   # static unless temp_dependent_k
        Et = cfg.Et
        self._Et2 = (np.asarray(Et)**2 if np.ndim(Et) else float(Et)**2)
        self._k_base = (g.k_depth.copy() if getattr(g, 'k_depth', None) is not None
                        else np.full(self.nz, cfg.k_dust))
        # base_grid.diag was already built with cp(T_bottom); start every column from it
        self._diag_cols = [np.array(g.diag, dtype=float) for _ in range(ncols)]
        self._T_last = self.T.reshape(ncols, self.nz).copy()
        self.vol.set_vertical_diag_cols(self._diag_cols)

    def _col_heat_cond(self, Tc):
        cfg = self.cfg
        if cfg.temp_dependent_cp:
            c0, c1, c2, c3, c4 = cfg.cp_coeffs
            heat_c = c0 + c1 * Tc + c2 * Tc**2 + c3 * Tc**3 + c4 * Tc**4
        else:
            heat_c = self.base.heat
        if cfg.temp_dependent_k:
            cond_c = self._k_base * (1.0 + cfg.k_temp_coeff * Tc**3) * self._Et2
        else:
            cond_c = self._cond_base
        return heat_c, cond_c

    def _update_operators(self):
        """Rebuild per-column vertical operators (and the lateral diffusivity) as T drifts."""
        if not self._temp_dep:
            return
        cfg = self.cfg
        flatT = self.T.reshape(self.nx * self.ny, self.nz)
        changed = False
        for c in range(flatT.shape[0]):
            if np.max(np.abs(flatT[c] - self._T_last[c])) > cfg.temp_change_threshold:
                heat_c, cond_c = self._col_heat_cond(flatT[c])
                self._diag_cols[c] = build_vertical_diag(
                    self._lthick, self._dens, heat_c, cond_c, self.base.dt)
                if cfg.temp_dependent_k:
                    # track this column's node-1 conductivity in its surface BC, as the 1D
                    # model does via the live grid.cond[1] (modelmain._T_surf_calc).
                    i, j = divmod(c, self.ny)
                    self._k_dx[i, j] = cond_c[1] / self._Et1**2 / self._dx_surf
                self._T_last[c] = flatT[c].copy()
                changed = True
        if changed and not self.vol._lateral_off:
            # lateral operators use a laterally-representative diffusivity K(z) = cond/(dens*heat)
            if cfg.temp_dependent_cp:
                c0, c1, c2, c3, c4 = cfg.cp_coeffs
                heat_all = c0 + c1 * self.T + c2 * self.T**2 + c3 * self.T**3 + c4 * self.T**4
            else:
                heat_all = self.base.heat[None, None, :]
            if cfg.temp_dependent_k:
                cond_all = self._k_base * (1.0 + cfg.k_temp_coeff * self.T**3) * self._Et2
            else:
                cond_all = self._cond_base[None, None, :]
            K_mean = (cond_all / (self._dens * heat_all)).mean(axis=(0, 1))
            self.vol.set_lateral_diffusivity(K_mean)

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
                self._update_operators()                     # temp-dependent props (no-op if off)
                self.T = self.vol.step(self.T, 0.0)          # non-RTE: no volumetric source
                self._surface_bc(self.mu_array[j], self.F_array[j])
            if record_surf:
                surf_hist.append(self.T_surf.copy())
        if record_surf:
            self.surf_hist = np.array(surf_hist)
        return self.T, self.T_surf
