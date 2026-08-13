"""Emergent thermal spectra for the 3D model, batched over all columns.

Radiation in regolith is vertical, so the emergent spectrum is computed one column at a time --
but DISORT batches columns (n_cols), so a single hybrid-thermal output solve over n_cols = nx*ny
gives the per-band radiance for the whole surface at once. This is the 3D analogue of
radiance_processor.calculate_radiances_from_results (hybrid, thermal_only), and reduces to it
column-for-column.

Requires the thermal optical constants (cfg.mie_file_out / cfg.wn_bounds_out, e.g. the enstatite
set). Hybrid-thermal is Planck-only, so no solar-spectrum file is needed.
"""
import numpy as np

_H = 6.62607015e-34   # J s
_C = 2.99792458e8     # m/s
_K = 1.380649e-23     # J/K


def band_brightness_temperature(wn_lo, wn_hi, L_band):
    """Per-band brightness temperature (K) from band-integrated radiance L_band (W/m^2/sr).

    Bands are narrow (~1 cm^-1 for the enstatite set), so L_band / bandwidth is the spectral
    radiance at band centre and the Planck law inverts in closed form. Vectorized over bands;
    non-positive radiances return NaN.
    """
    wn_lo = np.asarray(wn_lo, float)
    wn_hi = np.asarray(wn_hi, float)
    L = np.asarray(L_band, float)
    wn_c_m = 0.5 * (wn_lo + wn_hi) * 100.0            # band-centre wavenumber, m^-1
    dwn_m = (wn_hi - wn_lo) * 100.0                   # band width, m^-1
    B = np.where(L > 0, L / dwn_m, np.nan)            # spectral radiance W/m^2/sr/(m^-1)
    arg = 1.0 + (2.0 * _H * _C**2 * wn_c_m**3) / B
    return (_H * _C * wn_c_m) / (_K * np.log(arg))


def compute_spectra(cfg, base_grid, T_field, observer_mu=1.0):
    """Emergent thermal spectrum for every column of a 3D temperature field.

    Parameters
    ----------
    cfg : SimulationConfig  (output_radiance_mode must permit hybrid; optics via *_out paths)
    base_grid : LayerGrid   (the shared vertical grid)
    T_field : ndarray [nx, ny, nz]
    observer_mu : float     (cosine of the emission angle; 1.0 = nadir)

    Returns
    -------
    wavenumbers : [nwave] band centres (cm^-1)
    radiance    : [nx, ny, nwave] band-integrated thermal radiance (W/m^2/sr)
    BT          : [nx, ny, nwave] per-band brightness temperature (K)
    """
    from rte_disort import DisortRTESolver

    nx, ny, nz = T_field.shape
    ncols = nx * ny
    solver = DisortRTESolver(
        cfg, base_grid, n_cols=ncols, output_radiance=True, planck=True,
        observer_mu=observer_mu, solver_mode='hybrid', spectral_component='thermal_only')

    Tflat = np.ascontiguousarray(T_field.reshape(ncols, nz).T)   # [nz, ncols]
    rad, _ = solver.disort_run(Tflat, 0.0, 0.0)                  # thermal only: no sun
    rad = np.asarray(rad)                                        # [nwave, ncols, n_mu, n_phi]
    rad = rad[:, :, 0, 0]                                        # first observer -> [nwave, ncols]

    wn = np.asarray(solver.wavenumbers, float)                  # [nwave]
    lo = np.asarray(solver.lower_wns, float)                    # band edges the solver used
    hi = np.asarray(solver.upper_wns, float)

    radiance = rad.T                                            # [ncols, nwave]
    BT = np.stack([band_brightness_temperature(lo, hi, radiance[c]) for c in range(ncols)])
    return wn, radiance.reshape(nx, ny, -1), BT.reshape(nx, ny, -1)
