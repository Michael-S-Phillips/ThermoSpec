"""Observer-geometry per-facet brightness-temperature spectra for terrain (DEM) runs.

After a crater/terrain run the facet temperature columns are `sim.T_crater_out[depth, facet, time]`
on the shared vertical grid `sim.grid`. `TerrainObserver` computes, for a fixed observer direction,
each facet's emergent thermal spectrum and per-band brightness temperature over the output times,
evaluated at that facet's LOCAL emission angle (cos = facet_normal . observer) and masked to facets
that face and are not terrain-occluded from the observer.

Reuses the validated spectral path: one batched DISORT hybrid-thermal output solver (built once,
reused across time) + `radiance3d.band_brightness_temperature`. Nadir (observer = mesh +z = up) is
the default, matching a first orbital-geometry comparison; any observer direction is supported.

Primary output (CS-requested): `cube(...)` -> BT[n_facets, n_bands, n_out] on the enstatite band
grid, plus the band wavenumbers/wavelengths and the observer geometry used.
"""
import numpy as np

from radiance3d import band_brightness_temperature


def facet_view_geometry(mesh, observer_vec):
    """(mu_obs, occlusion_fraction, visible) per facet toward `observer_vec`.

    mu_obs = facet_normal . observer (cosine of local emission angle). `visible` = faces the
    observer (mu_obs>0) and not terrain-occluded (ShadowTester ray engine along observer dir)."""
    from crater import ShadowTester
    obs = np.asarray(observer_vec, float)
    obs = obs / np.linalg.norm(obs)
    mu_obs = np.asarray(mesh.normals, float) @ obs
    occl = ShadowTester(mesh).illuminated_facets(obs)
    visible = (mu_obs > 1e-6) & (occl > 1e-6)
    return mu_obs, occl, visible


class TerrainObserver:
    """Per-facet emergent brightness temperature toward a fixed observer, reusable across times."""

    def __init__(self, cfg, base_grid, mesh, observer_vec=None, mu_grid=None,
                 bands=None, band_idx=None):
        """observer_vec: view direction (default nadir = mesh +z). mu_grid: emission-angle cosines
        DISORT solves; each facet is interpolated to its own mu_obs (a single value -> nadir fast
        path, no interpolation). bands / band_idx (Ask 2): restrict the DISORT output solve to a few
        spectral bands -- `band_idx` are explicit indices into the thermal band grid, or `bands` are
        target wavelengths in um mapped to the nearest bands. The BT cost scales with the band
        count, so this is a ~100x speedup for a few-band Diviner comparison."""
        from rte_disort import DisortRTESolver
        self.observer_vec = (np.array([0.0, 0.0, 1.0]) if observer_vec is None
                             else np.asarray(observer_vec, float))
        self.n_facets = len(mesh.normals)
        self.mu_obs, self.occlusion, self.visible = facet_view_geometry(mesh, self.observer_vec)
        self.mu_grid = np.linspace(0.05, 1.0, 12) if mu_grid is None else np.atleast_1d(np.asarray(mu_grid, float))
        self._single_mu = len(self.mu_grid) == 1             # nadir/single-angle fast path (no interp)

        # Solve ONLY the observer-visible facets (skip the NaN ones) -- the DISORT output solve is
        # the dominant cost, so this scales it with the visible count, not the total facet count.
        self._vis_idx = np.where(self.visible)[0]
        nvis = len(self._vis_idx)
        self._pad = nvis == 1                                # DISORT's n_cols=1 path is unsupported
        ncols = 2 if self._pad else nvis

        def _build(ncol):
            s = DisortRTESolver(cfg, base_grid, n_cols=ncol, output_radiance=True, planck=True,
                                observer_mu=self.mu_grid, solver_mode='hybrid',
                                spectral_component='thermal_only')
            if bands is not None or band_idx is not None:
                s.restrict_to_bands(self._resolve_band_idx(s.wavenumbers, bands, band_idx))
            return s

        self._solver = None
        if ncols >= 2:
            self._solver = _build(ncols)
            src = self._solver
        else:
            # nothing visible; still need the band grid for output shapes
            src = _build(2)
        self.wavenumbers = np.asarray(src.wavenumbers, float)   # cm^-1
        self._lo = np.asarray(src.lower_wns, float)
        self._hi = np.asarray(src.upper_wns, float)
        self.wavelengths_um = 1.0e4 / self.wavenumbers

    @staticmethod
    def _resolve_band_idx(wavenumbers, bands, band_idx):
        """Sorted, unique band indices from explicit `band_idx` or target wavelengths `bands` (um)."""
        wn = np.asarray(wavenumbers, float)
        if band_idx is not None:
            idx = np.atleast_1d(np.asarray(band_idx, int))
        else:
            wl = 1.0e4 / wn                                  # um per band
            targets = np.atleast_1d(np.asarray(bands, float))
            idx = np.array([int(np.argmin(np.abs(wl - t))) for t in targets], dtype=int)
        return np.unique(idx)                                # sorted ascending in wavenumber

    def brightness_temperature(self, T_facets):
        """Single time: T_facets [nz, n_facets] -> (radiance, BT) each [n_facets, n_bands],
        NaN on facets not visible to the observer."""
        nb = len(self.wavenumbers)
        radiance = np.full((self.n_facets, nb), np.nan)
        BT = np.full((self.n_facets, nb), np.nan)
        if self._solver is None:
            return radiance, BT
        T = np.ascontiguousarray(np.asarray(T_facets, float)[:, self._vis_idx])   # visible only
        if self._pad:
            T = np.repeat(T, 2, axis=1)                                            # [nz,2]
        rad, _ = self._solver.disort_run(T, 0.0, 0.0)
        rad = np.asarray(rad)[:, :, 0, :]                    # [nwave, ncols, n_mu] (phi=0, all mu)
        for c, f in enumerate(self._vis_idx):
            if self._single_mu:
                # nadir/single-angle fast path: DISORT already solved at the one requested mu,
                # so there is nothing to interpolate (and no zero-width interval to divide by).
                radiance[f] = rad[:, c, 0]
            else:
                k = int(np.clip(np.searchsorted(self.mu_grid, self.mu_obs[f]), 1, len(self.mu_grid) - 1))
                w = (self.mu_obs[f] - self.mu_grid[k - 1]) / (self.mu_grid[k] - self.mu_grid[k - 1])
                radiance[f] = rad[:, c, k - 1] * (1.0 - w) + rad[:, c, k] * w
            BT[f] = band_brightness_temperature(self._lo, self._hi, radiance[f])
        return radiance, BT

    def cube(self, T_crater, time_indices=None):
        """T_crater [nz, n_facets, n_times] (or [nz, n_facets]) -> dict with
        BT [n_facets, n_bands, n_out], radiance [n_facets, n_bands, n_out], wavenumbers,
        wavelengths_um, mu_obs, visible, observer_vec."""
        T = np.asarray(T_crater, float)
        if T.ndim == 2:
            T = T[:, :, None]
        nt = T.shape[2]
        times = range(nt) if time_indices is None else list(time_indices)
        nb = len(self.wavenumbers)
        BT = np.full((self.n_facets, nb, len(times)), np.nan)
        RAD = np.full((self.n_facets, nb, len(times)), np.nan)
        for i, t in enumerate(times):
            RAD[:, :, i], BT[:, :, i] = self.brightness_temperature(T[:, :, t])
        return dict(BT=BT, radiance=RAD, wavenumbers=self.wavenumbers,
                    wavelengths_um=self.wavelengths_um, mu_obs=self.mu_obs, visible=self.visible,
                    observer_vec=self.observer_vec)


def terrain_bt_cube(cfg, base_grid, mesh, T_crater, observer_vec=None, time_indices=None,
                    mu_grid=None, bands=None, band_idx=None):
    """Convenience: build a TerrainObserver and return its BT cube. See TerrainObserver.cube."""
    obs = TerrainObserver(cfg, base_grid, mesh, observer_vec=observer_vec, mu_grid=mu_grid,
                          bands=bands, band_idx=band_idx)
    return obs.cube(T_crater, time_indices=time_indices)
