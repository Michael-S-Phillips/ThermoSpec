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
        self.cfg = cfg
        self._use_rte = cfg.use_RTE
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
        self._setup_rte()

    def _setup_rte(self):
        """Set up per-column radiative transfer (DISORT, batched over all columns).

        RTE in regolith is 1D per column, so all nx*ny columns are solved in one DISORT call
        with n_cols = nx*ny (the same batching the crater model uses over facets). Currently
        supports RTE_solver='disort', thermal_evolution_mode='two_wave', single_layer.
        """
        if not self._use_rte:
            return
        if not self.cfg.single_layer:
            raise NotImplementedError("Simulator3D RTE currently supports single_layer=True.")
        ncols = self.nx * self.ny
        if self.cfg.RTE_solver == 'disort':
            from rte_disort import DisortRTESolver
            mode = self.cfg.thermal_evolution_mode
            if mode == 'two_wave':
                self._rte = DisortRTESolver(self.cfg, self.base, n_cols=ncols,
                                            output_radiance=False, planck=True, solver_mode='two_wave')
                self._rte_vis = DisortRTESolver(self.cfg, self.base, n_cols=ncols,
                                                output_radiance=False, planck=False, solver_mode='two_wave')
            elif mode == 'hybrid':
                # spectral thermal (multi-wave, Planck-only -> no solar file) + broadband visible.
                # EXPENSIVE: the thermal solve is ~850x two_wave (one DISORT per wavenumber band),
                # and cost scales with column count -- practical only for small grids / short runs,
                # or as ground truth for a learned per-column source surrogate.
                self._rte = DisortRTESolver(self.cfg, self.base, n_cols=ncols, output_radiance=False,
                                            planck=True, solver_mode='hybrid',
                                            spectral_component='thermal_only')
                self._rte_vis = DisortRTESolver(self.cfg, self.base, n_cols=ncols, output_radiance=False,
                                                planck=False, solver_mode='hybrid',
                                                spectral_component='visible_only')
            elif mode == 'multi_wave':
                raise NotImplementedError(
                    "multi_wave thermal evolution needs solar-spectrum files, which are absent "
                    "from the repo; use thermal_evolution_mode='hybrid' or 'two_wave'.")
            else:
                raise NotImplementedError(f"Unknown thermal_evolution_mode: {mode}")
        elif self.cfg.RTE_solver == 'hapke':
            # Hapke is scalar-per-column (broadband); loop columns with per-column BVP state,
            # exactly as the 1D model and its crater path do (modelmain.py:841, :992).
            from rte_hapke import RadiativeTransfer
            self._rte_hapke = RadiativeTransfer(self.cfg, self.base)
            nlb = self.base.nlay_dust + 1
            self._phi_vis = np.zeros((ncols, nlb))
            self._phi_therm = np.zeros((ncols, nlb))
        else:
            raise NotImplementedError(f"Unknown RTE_solver: {self.cfg.RTE_solver}")

    def _rte_step(self, j):
        if self.cfg.RTE_solver == 'disort':
            return self._rte_step_disort(j)
        return self._rte_step_hapke(j)

    def _rte_step_disort(self, j):
        """Run DISORT for all columns at this step; return the 3D source and set self.T_surf.

        Mirrors the DISORT branch of modelmain.Simulator.run (two_wave): a thermal (planck)
        solve gives the flux-divergence source and upward flux (-> brightness T_surf), and a
        visible solve adds the absorbed-solar source when any column is sunlit.
        """
        ncols, nz = self.nx * self.ny, self.nz
        mu_col = np.ascontiguousarray((self.mu_array[j] * self.mu_fac).reshape(ncols))
        F_col = np.ascontiguousarray((self.F_array[j] * self.F_gate).reshape(ncols))
        Tflat = np.ascontiguousarray(self.T.reshape(ncols, nz).T)          # [nz, ncols]

        src_th, flup_th = self._rte.disort_run(Tflat, mu_col.copy(), F_col.copy())
        src_th = np.asarray(src_th)                                        # [ncols, nz]
        flup_th = np.asarray(flup_th)
        # two_wave gives upward flux [ncols]; hybrid/multi-wave gives [nwave, ncols] -> sum bands
        flup_col = flup_th.sum(axis=0) if flup_th.ndim == 2 else flup_th.reshape(ncols)
        if np.any(F_col > 0.001):
            src_vis, _ = self._rte_vis.disort_run(Tflat, mu_col.copy(), F_col.copy())
            src = src_th + np.asarray(src_vis)
        else:
            src = src_th
        self.T_surf = (flup_col / self.cfg.sigma).reshape(self.nx, self.ny) ** 0.25
        return src.reshape(self.nx, self.ny, nz)

    def _rte_step_hapke(self, j):
        """Hapke RTE per column (loop; scalar-per-column solver). Returns the 3D source and sets
        self.T_surf. Mirrors modelmain's Hapke branch: compute_source returns a combined
        (visible+thermal) volumetric source and updates the per-column BVP state; the upward
        thermal flux is phi_therm[0]*2*pi -> brightness T_surf."""
        ncols, nz = self.nx * self.ny, self.nz
        mu_col = (self.mu_array[j] * self.mu_fac).reshape(ncols)
        F_col = (self.F_array[j] * self.F_gate).reshape(ncols)
        Tcols = self.T.reshape(ncols, nz)                        # row c = full column profile
        x_RTE, x_b = self.base.x_RTE, self.base.x_boundaries
        source = np.zeros((ncols, nz))
        Tsurf = np.empty(ncols)
        two_pi_sigma = 2.0 * np.pi / self.cfg.sigma
        for c in range(ncols):
            s, self._phi_vis[c], self._phi_therm[c] = self._rte_hapke.compute_source(
                Tcols[c], x_RTE, x_b, self._phi_vis[c], self._phi_therm[c],
                float(mu_col[c]), float(F_col[c]))
            source[c] = s
            Tsurf[c] = (self._phi_therm[c][0] * two_pi_sigma) ** 0.25
        self.T_surf = Tsurf.reshape(self.nx, self.ny)
        return source.reshape(self.nx, self.ny, nz)

    def _rte_bc(self):
        """RTE surface/bottom BC (modelmain._bc, use_RTE branch): Neumann top, then bottom."""
        self.T[:, :, 0] = self.T[:, :, 1]
        if self.cfg.bottom_bc == 'dirichlet':
            self.T[:, :, -1] = self.cfg.T_bottom
        else:
            self.T[:, :, -1] = self.T[:, :, -2]

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

    def run(self, record_surf=False, record_phases=False):
        """Advance the diurnal run.

        record_phases: keep the full temperature field at diurnal noon and pre-dawn (the steps of
        maximum and minimum mean surface brightness temperature) in self.T_noon / self.T_predawn,
        for later spectra via phase_spectra().
        """
        surf_hist = []
        if record_phases:
            self.T_noon = self.T_predawn = None
            self._best_max, self._best_min = -np.inf, np.inf
        for j in range(self.t_num):
            if j > 0:
                self._update_operators()                     # temp-dependent props (no-op if off)
                if self._use_rte:
                    source = self._rte_step(j)               # per-column RTE source + T_surf
                    self.T = self.vol.step(self.T, source)
                    self._rte_bc()
                else:
                    self.T = self.vol.step(self.T, 0.0)      # non-RTE: no volumetric source
                    self._surface_bc(self.mu_array[j], self.F_array[j])
            if record_surf:
                surf_hist.append(self.T_surf.copy())
            if record_phases and j > 0:
                m = float(np.mean(self.T_surf))
                if m > self._best_max:
                    self._best_max, self.T_noon = m, self.T.copy()
                    self.Tsurf_noon = np.asarray(self.T_surf).copy()
                if m < self._best_min:
                    self._best_min, self.T_predawn = m, self.T.copy()
                    self.Tsurf_predawn = np.asarray(self.T_surf).copy()
        if record_surf:
            self.surf_hist = np.array(surf_hist)
        return self.T, self.T_surf

    def phase_spectra(self, phase='noon', observer_mu=1.0):
        """Emergent thermal spectra (wavenumbers, radiance, per-band brightness T) for every
        column, from the recorded noon or pre-dawn field. Requires run(record_phases=True) and
        thermal optics in cfg (mie_file_out / wn_bounds_out). See radiance3d.compute_spectra."""
        from radiance3d import compute_spectra
        field = {'noon': getattr(self, 'T_noon', None),
                 'predawn': getattr(self, 'T_predawn', None)}.get(phase)
        if field is None:
            raise ValueError("No recorded field; run(record_phases=True) first (phase='%s')." % phase)
        return compute_spectra(self.cfg, self.base, field, observer_mu=observer_mu)
