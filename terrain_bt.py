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

    def __init__(self, cfg, base_grid, mesh, observer_vec=None, mu_grid=None):
        from rte_disort import DisortRTESolver
        self.observer_vec = (np.array([0.0, 0.0, 1.0]) if observer_vec is None
                             else np.asarray(observer_vec, float))
        self.n_facets = len(mesh.normals)
        self.mu_obs, self.occlusion, self.visible = facet_view_geometry(mesh, self.observer_vec)
        self.mu_grid = np.linspace(0.05, 1.0, 12) if mu_grid is None else np.asarray(mu_grid, float)
        self._solver = DisortRTESolver(cfg, base_grid, n_cols=self.n_facets, output_radiance=True,
                                       planck=True, observer_mu=self.mu_grid, solver_mode='hybrid',
                                       spectral_component='thermal_only')
        self.wavenumbers = np.asarray(self._solver.wavenumbers, float)       # cm^-1
        self.wavelengths_um = 1.0e4 / self.wavenumbers
        self._lo = np.asarray(self._solver.lower_wns, float)
        self._hi = np.asarray(self._solver.upper_wns, float)

    def brightness_temperature(self, T_facets):
        """Single time: T_facets [nz, n_facets] -> (radiance, BT) each [n_facets, n_bands],
        NaN on facets not visible to the observer."""
        T_facets = np.ascontiguousarray(np.asarray(T_facets, float))
        rad, _ = self._solver.disort_run(T_facets, 0.0, 0.0)
        rad = np.asarray(rad)[:, :, 0, :]                    # [nwave, ncols, n_mu] (phi=0, all mu)
        nb = len(self.wavenumbers)
        radiance = np.full((self.n_facets, nb), np.nan)
        for f in np.where(self.visible)[0]:
            k = int(np.clip(np.searchsorted(self.mu_grid, self.mu_obs[f]), 1, len(self.mu_grid) - 1))
            w = (self.mu_obs[f] - self.mu_grid[k - 1]) / (self.mu_grid[k] - self.mu_grid[k - 1])
            radiance[f] = rad[:, f, k - 1] * (1.0 - w) + rad[:, f, k] * w
        BT = np.full((self.n_facets, nb), np.nan)
        for f in np.where(self.visible)[0]:
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
                    mu_grid=None):
    """Convenience: build a TerrainObserver and return its BT cube. See TerrainObserver.cube."""
    obs = TerrainObserver(cfg, base_grid, mesh, observer_vec=observer_vec, mu_grid=mu_grid)
    return obs.cube(T_crater, time_indices=time_indices)
