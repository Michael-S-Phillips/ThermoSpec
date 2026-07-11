"""
Tools for setting up hybrid-mode (multi-wave thermal + broadband visible) DISORT
runs consistently with a two-wave (broadband) reference run, without touching
core model code.

Two problems this addresses:

1. Grid spacing is built once from cfg.Et (the "thermal" extinction
   coefficient), in tau units, before any Mie file is read (see grid.py). Even
   in hybrid mode, every wavelength band shares that same grid - DISORT only
   rescales the per-layer optical depth per wavelength via the Mie Cext
   values (rte_disort.py:_setup_optical_properties_advanced), it never
   changes the number of layers or their relative spacing. So a representative
   scalar Et must still be chosen from the (wavelength-dependent) Mie data to
   get sensible grid spacing. set_thermal_Et_from_mie() does this by averaging
   Et(lambda) = n_p * Cext(lambda) * 1e-12 over the wavelengths that will
   actually survive hybrid_wavelength_cutoff filtering at runtime (i.e., the
   subset DisortRTESolver._filter_thermal_wavelengths keeps for the thermal
   band), matching what the simulation will really use.

2. eta (visible/thermal extinction ratio) sets the visible-band grid spacing
   via tau_boundaries *= cfg.eta in rte_disort.py:_setup_optical_properties.
   set_eta_geometric_optics() computes the visible extinction coefficient from
   pure geometric optics (Qext=2 by default): Et_vis = 3*fill_frac*Qext/(4*radius),
   then sets cfg.eta = Et_vis / cfg.Et.

3. To verify two_wave and hybrid modes agree when given equivalent optical
   properties, write_flat_mie_file() builds a companion Mie file with the same
   wavelength/wavenumber grid as a real Mie file but constant g/Cext/ssalb.
   Paired with cfg.scale_Et=True (an existing config.py option), DISORT
   rescales every band's Et(lambda) so their mean - and since the file is
   flat, every individual band - exactly equals cfg.Et. Combined with constant
   ssalb/g, this makes every hybrid-thermal band physically identical to a
   two_wave band, with no core-code changes required.
"""
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import copy
import numpy as np

from modelmain import Simulator
from config import SimulationConfig


def load_mie_file(mie_file):
    """Load a Mie-code output file and return arrays sorted by wavenumber,
    matching the convention used internally by
    rte_disort.DisortRTESolver._load_constants.

    Returns
    -------
    wavenumbers, g, Cext, Csca, ssalb : ndarray
    """
    mie_params = np.loadtxt(mie_file)
    sortidx = np.argsort(10000.0 / mie_params[:, 0])
    wavenumbers = 10000.0 / mie_params[sortidx, 0]
    g = mie_params[sortidx, 1]
    Cext = mie_params[sortidx, 2]
    Csca = mie_params[sortidx, 3]
    ssalb = mie_params[sortidx, 4]
    return wavenumbers, g, Cext, Csca, ssalb


def thermal_mask(wavenumbers, hybrid_wavelength_cutoff):
    """Boolean mask matching DisortRTESolver._filter_thermal_wavelengths:
    keeps wavenumbers <= cutoff_wn (i.e. wavelengths >= hybrid_wavelength_cutoff)."""
    cutoff_wn = 10000.0 / hybrid_wavelength_cutoff
    return wavenumbers <= cutoff_wn


def mie_Et_spectrum(wavenumbers, Cext, fill_frac, radius):
    """Per-wavelength thermal extinction coefficient Et(lambda), using the same
    n_p*Cext formula as DisortRTESolver._setup_optical_properties_advanced
    (uniform fill_frac / particle radius, no depth dependence)."""
    Vp = (4.0 / 3.0) * np.pi * radius ** 3
    n_p = fill_frac / Vp
    return n_p * Cext * 1e-12


def set_thermal_Et_from_mie(cfg: SimulationConfig, mie_file=None, fill_frac=None,
                             radius=None, method='mean'):
    """Set cfg.Et from the Mie file's per-wavelength thermal extinction
    coefficients, restricted to the wavelengths that survive
    hybrid_wavelength_cutoff filtering (the subset actually used as the
    'thermal' band in hybrid mode). Mutates cfg in place and returns it.

    Parameters
    ----------
    cfg : SimulationConfig
    mie_file : str, optional
        Defaults to cfg.mie_file.
    fill_frac, radius : float, optional
        Default to cfg.fill_frac / cfg.radius.
    method : 'mean' or 'max'
    """
    mie_file = mie_file or cfg.mie_file
    fill_frac = fill_frac if fill_frac is not None else cfg.fill_frac
    radius = radius if radius is not None else cfg.radius

    wavenumbers, g, Cext, Csca, ssalb = load_mie_file(mie_file)
    mask = thermal_mask(wavenumbers, cfg.hybrid_wavelength_cutoff)
    if not np.any(mask):
        raise ValueError(
            f"No wavelengths in {mie_file} survive hybrid_wavelength_cutoff="
            f"{cfg.hybrid_wavelength_cutoff} um filtering."
        )
    Et_array = mie_Et_spectrum(wavenumbers[mask], Cext[mask], fill_frac, radius)

    if method == 'mean':
        cfg.Et = float(np.mean(Et_array))
    elif method == 'max':
        cfg.Et = float(np.max(Et_array))
    else:
        raise ValueError("method must be 'mean' or 'max'")

    cfg.__post_init__()
    return cfg


def set_eta_geometric_optics(cfg: SimulationConfig, Qext=2.0, fill_frac=None, radius=None):
    """Set cfg.eta = Et_visible / cfg.Et, where Et_visible is computed from
    pure geometric-optics extinction (Qext=2 by default):

        Et_vis = 3 * fill_frac * Qext / (4 * radius)

    Requires cfg.Et to already reflect the desired thermal reference value
    (e.g. via set_thermal_Et_from_mie, called first). Mutates cfg in place
    and returns (cfg, Et_vis).
    """
    fill_frac = fill_frac if fill_frac is not None else cfg.fill_frac
    radius = radius if radius is not None else cfg.radius

    Et_vis = 3.0 * fill_frac * Qext / (4.0 * radius)
    cfg.eta = Et_vis / cfg.Et
    cfg.__post_init__()
    return cfg, Et_vis


def thermal_mean_properties(mie_file, hybrid_wavelength_cutoff):
    """Simple (unweighted) mean ssalb and g over the thermal-band subset of a
    Mie file, for use as two_wave broadband reference values.

    Returns
    -------
    ssalb_mean, g_mean : float
    """
    wavenumbers, g, Cext, Csca, ssalb = load_mie_file(mie_file)
    mask = thermal_mask(wavenumbers, hybrid_wavelength_cutoff)
    return float(np.mean(ssalb[mask])), float(np.mean(g[mask]))


def write_thermal_wn_bounds_file(wn_bounds_file, mie_file, hybrid_wavelength_cutoff, output_path):
    """Build a wn_bounds file matching the thermal-only wavelength subset that
    DisortRTESolver._filter_thermal_wavelengths keeps for hybrid mode.

    Gotcha this works around: DisortRTESolver._compute_wn_bounds loads
    wn_bounds fresh from disk and requires it to have exactly
    len(thermal_wavenumbers)+1 rows, but the wn_bounds file itself is NOT
    automatically filtered to match hybrid_wavelength_cutoff the way the Mie
    arrays are. If your Mie file spans wavelengths shorter than
    hybrid_wavelength_cutoff (e.g. serpentine_mie_200wns_1um.txt goes down to
    1.0 um vs. the 3.33 um default cutoff), constructing a 'hybrid' /
    'thermal_only' DisortRTESolver with the original, unfiltered wn_bounds
    file raises a ValueError on a bin-count mismatch. Any real hybrid-mode
    run using this Mie file needs a wn_bounds file produced by this function
    (matching whatever hybrid_wavelength_cutoff you're actually using), not
    just the two_wave/hybrid equivalence verification below.
    """
    wavenumbers, g, Cext, Csca, ssalb = load_mie_file(mie_file)
    mask = thermal_mask(wavenumbers, hybrid_wavelength_cutoff)
    if not np.all(np.diff(mask.astype(int)) <= 0):
        raise ValueError(
            "Thermal mask is not a contiguous prefix in sorted-wavenumber "
            "order; cannot safely slice wn_bounds. Check that wn_bounds_file "
            "and mie_file share the same wavelength grid."
        )
    n_keep = int(mask.sum())

    wn_bounds = np.sort(np.loadtxt(wn_bounds_file))
    if len(wn_bounds) != len(wavenumbers) + 1:
        raise ValueError(
            f"{wn_bounds_file} has {len(wn_bounds)} edges, expected "
            f"{len(wavenumbers) + 1} to match {mie_file} ({len(wavenumbers)} rows)."
        )
    filtered_edges = wn_bounds[:n_keep + 1]
    np.savetxt(output_path, filtered_edges, fmt='%.6e')
    return output_path


def write_flat_mie_file(mie_file, output_path, ssalb_const, g_const, Cext_const=1.0):
    """Write a companion Mie file with the identical wavelength grid as
    `mie_file` but constant g/Cext/Csca/ssalb across all rows.

    Cext_const is an arbitrary positive value: as long as cfg.scale_Et=True is
    used downstream, DISORT renormalizes the mean Et(lambda) to cfg.Et
    regardless of the absolute Cext scale, so the specific Cext_const value
    doesn't matter.
    """
    mie_params = np.loadtxt(mie_file)
    wavelengths = mie_params[:, 0]
    Csca_const = ssalb_const * Cext_const
    out = np.column_stack([
        wavelengths,
        np.full_like(wavelengths, g_const),
        np.full_like(wavelengths, Cext_const),
        np.full_like(wavelengths, Csca_const),
        np.full_like(wavelengths, ssalb_const),
    ])
    np.savetxt(output_path, out, fmt='%.6f\t%.6f\t%.6f\t%.6f\t%.6f')
    return output_path


def compare_two_wave_vs_hybrid(cfg: SimulationConfig, flat_mie_file, ssalb_therm,
                                g_therm, T_profile=None, mu=0.0, F=0.0):
    """Run a single RTE flux/source-term calculation in two_wave mode and in
    hybrid mode (with a flattened Mie table matching ssalb_therm/g_therm and
    scale_Et=True), for direct comparison. No time-stepping - this only
    exercises the RTE solve itself, isolated from time-integration effects.

    Parameters
    ----------
    cfg : SimulationConfig
        Base configuration (not mutated). Must already have cfg.Et set to the
        desired reference value (e.g. via set_thermal_Et_from_mie).
    flat_mie_file : str
        Path to a companion Mie file from write_flat_mie_file(), sharing
        cfg.wn_bounds's wavelength/wavenumber grid.
    ssalb_therm, g_therm : float
        Values used both for the two_wave config and as the constants baked
        into flat_mie_file (caller is responsible for consistency).
    T_profile : ndarray, optional
        Full-grid temperature profile (length sim.grid.x_num). Defaults to
        isothermal at cfg.T_bottom.
    mu, F : float
        Solar incidence cosine / flux flag passed to disort_run (0, 0 by
        default: pure thermal self-emission, no solar contribution).

    Returns
    -------
    dict with keys: source_two_wave, source_hybrid, flux_up_two_wave,
    flux_up_hybrid, flux_up_hybrid_total
    """
    cfg_2w = copy.deepcopy(cfg)
    cfg_2w.crater = False
    cfg_2w.thermal_evolution_mode = 'two_wave'
    cfg_2w.ssalb_therm = ssalb_therm
    cfg_2w.g_therm = g_therm
    cfg_2w.__post_init__()
    sim_2w = Simulator(cfg_2w)

    cfg_hy = copy.deepcopy(cfg)
    cfg_hy.crater = False
    cfg_hy.thermal_evolution_mode = 'hybrid'
    cfg_hy.mie_file = flat_mie_file
    cfg_hy.scale_Et = True
    cfg_hy.__post_init__()
    sim_hy = Simulator(cfg_hy)

    if T_profile is None:
        T_profile = np.full(sim_2w.grid.x_num, cfg.T_bottom)

    source_2w, flux_up_2w = sim_2w.rte_disort.disort_run(T_profile, mu, F)
    source_hy, flux_up_hy = sim_hy.rte_disort.disort_run(T_profile, mu, F)

    return {
        'source_two_wave': source_2w,
        'source_hybrid': source_hy,
        'flux_up_two_wave': flux_up_2w,
        'flux_up_hybrid': flux_up_hy,
        'flux_up_hybrid_total': np.sum(flux_up_hy),
    }


def run_diurnal_comparison(cfg: SimulationConfig, hybrid_mie_file, ssalb_therm, g_therm,
                            scale_Et=True, wn_bounds_file=None, run_two_wave=True,
                            observer_angle=0.0):
    """Run a full diurnal thermal simulation in hybrid mode (using
    hybrid_mie_file, any Mie table - flattened, real spectral, downsampled, or
    a different material) and, optionally, a two_wave reference run using
    ssalb_therm/g_therm, then return both complete results for direct
    comparison.

    Unlike compare_two_wave_vs_hybrid (a single no-time-stepping RTE flux
    calculation), this actually drives the full day-stepping loop - including
    the real sun cycling on and off every step - so it also exercises the
    sun-handling code paths in Simulator._setup_thermal_evolution_solvers and
    DisortRTESolver.disort_run, not just the isothermal/no-sun RTE math.

    Parameters
    ----------
    cfg : SimulationConfig
        Base configuration (not mutated), already carrying the real
        diurnal/time-stepping settings (ndays, freq_out, etc.) you want to
        test with. Should already have cfg.Et/cfg.eta set to the calibrated
        reference values (e.g. via set_thermal_Et_from_mie /
        set_eta_geometric_optics) and cfg.wn_bounds set to the pre-filtered
        thermal wn_bounds file for hybrid_mie_file's wavelength grid (unless
        overridden via wn_bounds_file).
    hybrid_mie_file : str
        Path to any Mie file to use for the hybrid thermal solver: a flat
        companion file from write_flat_mie_file(), the real spectral file, a
        downsampled version from downsample_mie_file(), or a different
        material's Mie file entirely.
    ssalb_therm, g_therm : float
        Broadband properties for the two_wave reference config. Only used if
        run_two_wave=True.
    scale_Et : bool
        Whether to set cfg.scale_Et=True on the hybrid config (rescales the
        mean of hybrid_mie_file's per-wavelength Et(lambda) to cfg.Et).
    wn_bounds_file : str, optional
        Overrides cfg.wn_bounds for the hybrid run - needed if hybrid_mie_file
        has a different row count than what cfg.wn_bounds was built for (e.g.
        a downsampled file).
    run_two_wave : bool
        If False, skip the two_wave run and only return hybrid results (e.g.
        when comparing two hybrid runs against each other, like a downsampled
        vs. full-resolution Mie file).
    observer_angle : float
        Emission angle (degrees from nadir) passed to
        bolometric_brightness_temperature_series/max_brightness_temperature_series.

    Returns
    -------
    dict with keys: sim_two_wave, sim_hybrid, T_surf_two_wave, T_surf_hybrid,
    t_out_two_wave, t_out_hybrid, max_abs_T_surf_diff, Tb_bol_two_wave,
    Tb_bol_hybrid, max_abs_Tb_bol_diff, Tb_max_two_wave, Tb_max_hybrid,
    max_abs_Tb_max_diff (two_wave keys are None if run_two_wave=False).

    Three different comparisons, not interchangeable - see the module-level
    comment above bolometric_brightness_temperature_series for the full
    explanation:
    - T_surf_*: Simulator's internal hemispheric-flux-based temperature. Fair,
      apples-to-apples two_wave-vs-hybrid energy-balance check.
    - Tb_bol_*: actual DISORT directional radiance at observer_angle, summed
      across all bands, then inverted. Also a fair, apples-to-apples
      two_wave-vs-hybrid check, but for the directional (not hemispheric)
      quantity - the natural generalization of the original two_wave-only
      `(radiance*pi/sigma)**0.25` notebook calculation.
    - Tb_max_*: max-across-bands "Christiansen feature" radiance-based
      temperature. NOT a fair two_wave-vs-hybrid check (two_wave has no bands
      to pick a peak from) - a real, separate finding about hybrid's spectral
      behavior. Expect Tb_max_hybrid to differ from T_surf/Tb_bol by design,
      not by error.
    """
    sim_2w = None
    T_surf_2w = None
    t_out_2w = None
    if run_two_wave:
        cfg_2w = copy.deepcopy(cfg)
        cfg_2w.crater = False
        cfg_2w.thermal_evolution_mode = 'two_wave'
        cfg_2w.output_radiance_mode = 'two_wave'
        cfg_2w.ssalb_therm = ssalb_therm
        cfg_2w.g_therm = g_therm
        cfg_2w.__post_init__()
        sim_2w = Simulator(cfg_2w)
        _, _, _, T_surf_2w, t_out_2w = sim_2w.run()

    cfg_hy = copy.deepcopy(cfg)
    cfg_hy.crater = False
    cfg_hy.thermal_evolution_mode = 'hybrid'
    cfg_hy.output_radiance_mode = 'hybrid'
    cfg_hy.mie_file = hybrid_mie_file
    cfg_hy.scale_Et = scale_Et
    if wn_bounds_file is not None:
        cfg_hy.wn_bounds = wn_bounds_file
    # Mirror the evolution-mode spectral files for the output/radiance step too,
    # so max_brightness_temperature_series() sees the same spectral resolution
    # that actually produced the temperature field (mie_file_out/wn_bounds_out
    # otherwise stay at whatever cfg originally had, which would silently use a
    # different - possibly mismatched - material/resolution for the observer
    # radiance calculation than what evolved the thermal state).
    cfg_hy.mie_file_out = cfg_hy.mie_file
    cfg_hy.wn_bounds_out = cfg_hy.wn_bounds
    cfg_hy.__post_init__()
    sim_hy = Simulator(cfg_hy)
    _, _, _, T_surf_hy, t_out_hy = sim_hy.run()

    max_abs_diff = (float(np.max(np.abs(T_surf_2w - T_surf_hy)))
                    if run_two_wave else None)

    # Directional radiance-based brightness temperatures - see the module-level
    # comment above bolometric_brightness_temperature_series for what each of
    # these does and doesn't tell you.
    Tb_bol_2w = None
    if run_two_wave:
        _, Tb_bol_2w = bolometric_brightness_temperature_series(sim_2w, observer_angle=observer_angle)
    _, Tb_bol_hy = bolometric_brightness_temperature_series(sim_hy, observer_angle=observer_angle)
    max_abs_Tb_bol_diff = float(np.max(np.abs(Tb_bol_2w - Tb_bol_hy))) if run_two_wave else None

    Tb_max_2w = None
    if run_two_wave:
        _, Tb_max_2w = max_brightness_temperature_series(sim_2w, observer_angle=observer_angle)
    _, Tb_max_hy = max_brightness_temperature_series(sim_hy, observer_angle=observer_angle)
    max_abs_Tb_max_diff = float(np.max(np.abs(Tb_max_2w - Tb_max_hy))) if run_two_wave else None

    return {
        'sim_two_wave': sim_2w,
        'sim_hybrid': sim_hy,
        'T_surf_two_wave': T_surf_2w,
        'T_surf_hybrid': T_surf_hy,
        't_out_two_wave': t_out_2w,
        't_out_hybrid': t_out_hy,
        'max_abs_T_surf_diff': max_abs_diff,
        'Tb_bol_two_wave': Tb_bol_2w,
        'Tb_bol_hybrid': Tb_bol_hy,
        'max_abs_Tb_bol_diff': max_abs_Tb_bol_diff,
        'Tb_max_two_wave': Tb_max_2w,
        'Tb_max_hybrid': Tb_max_hy,
        'max_abs_Tb_max_diff': max_abs_Tb_max_diff,
    }


def downsample_mie_file(mie_file, wn_bounds_file, n_bands, output_mie_path, output_wn_bounds_path):
    """Downsample a Mie file (and its companion wn_bounds file) to n_bands
    coarse bands, for faster hybrid-mode runs.

    Groups the wavenumber-sorted rows into n_bands contiguous, roughly-equal
    chunks (np.array_split) and replaces each chunk with a single row: the
    simple (unweighted) mean of g/Cext/Csca/ssalb and mean wavelength within
    that chunk. The new wn_bounds file's edges are picked directly from the
    original edges at each chunk boundary, so the coarse bins exactly span
    (no gaps/overlaps) the same total wavelength range as the original file.

    This is a simple unweighted-average downsampling, not a physically
    rigorous band-averaging (e.g. Planck- or flux-weighted) - treat results
    from a heavily downsampled file as a speed/accuracy trade-off exploration,
    not a substitute for the full-resolution spectral run.

    Parameters
    ----------
    mie_file, wn_bounds_file : str
        Original full-resolution files (wn_bounds_file must have
        len(mie rows)+1 edges).
    n_bands : int
        Number of coarse output bands.
    output_mie_path, output_wn_bounds_path : str
        Where to write the downsampled files.

    Returns
    -------
    (output_mie_path, output_wn_bounds_path)
    """
    mie_params = np.loadtxt(mie_file)
    sortidx = np.argsort(10000.0 / mie_params[:, 0])
    sorted_params = mie_params[sortidx]  # columns: wavelength, g, Cext, Csca, ssalb

    wn_bounds = np.sort(np.loadtxt(wn_bounds_file))
    if len(wn_bounds) != len(sorted_params) + 1:
        raise ValueError(
            f"{wn_bounds_file} has {len(wn_bounds)} edges, expected "
            f"{len(sorted_params) + 1} to match {mie_file} ({len(sorted_params)} rows)."
        )
    if n_bands < 1 or n_bands > len(sorted_params):
        raise ValueError(f"n_bands must be between 1 and {len(sorted_params)}")

    groups = np.array_split(np.arange(len(sorted_params)), n_bands)

    out_rows = []
    edge_indices = [0]
    for grp in groups:
        chunk = sorted_params[grp]
        out_rows.append([
            chunk[:, 0].mean(),  # wavelength (representative)
            chunk[:, 1].mean(),  # g
            chunk[:, 2].mean(),  # Cext
            chunk[:, 3].mean(),  # Csca
            chunk[:, 4].mean(),  # ssalb
        ])
        edge_indices.append(edge_indices[-1] + len(grp))

    np.savetxt(output_mie_path, np.array(out_rows), fmt='%.6f\t%.6f\t%.6f\t%.6f\t%.6f')
    new_edges = wn_bounds[edge_indices]
    np.savetxt(output_wn_bounds_path, new_edges, fmt='%.6e')
    return output_mie_path, output_wn_bounds_path


# ---------------------------------------------------------------------------
# Brightness temperature from actual DISORT radiance, as distinct from the
# model's internal hemispheric-flux-based T_surf.
#
# Three genuinely different quantities are in play here - don't conflate them:
#
# 1. T_surf (Simulator.run()'s internal state, e.g. T_surf_two_wave/
#    T_surf_hybrid above) = (flux_up/sigma)**0.25, where flux_up is DISORT's
#    actual *hemispherically-integrated* upward diffuse flux at the top of
#    the column (rte_disort.disort_run's fl_up return value, itself
#    result[...,0] from ds.forward() - a genuine, self-consistent RT output).
#    This is NOT a top-node-kinetic-temperature/opaque-surface proxy - it
#    already reflects whatever subsurface, semi-transparent emission DISORT
#    computed, correctly, for the bulk energy budget. (The actual top-node
#    kinetic temperature is the separate quantity self.T[0]/T_out[0,:].)
#    What T_surf lacks is (a) viewing-angle information - it's a hemispheric
#    integral - and (b) spectral information - in hybrid mode it's summed
#    across all bands to one bolometric number. It is the right quantity for
#    checking that two_wave and hybrid agree on overall RT energy balance.
#
# 2. bolometric_brightness_temperature_series() below takes DISORT's actual
#    *directional* radiance (gather_rad(), via
#    radiance_processor.calculate_radiances_from_results) at a chosen
#    observer angle, sums it across all spectral bands (each band's output is
#    already band-integrated, not a density - see the function docstring),
#    and inverts via (total_radiance*pi/sigma)**0.25. This is the direct,
#    fair generalization of the single-band `(radiance*pi/sigma)**0.25`
#    calculation used in the original two_wave-only notebook plotting cell -
#    same physical quantity, extended correctly to multi-wave/hybrid by
#    summing over bands first. Use this to compare two_wave vs. hybrid
#    directional brightness temperature on equal footing.
#
# 3. max_brightness_temperature_series() below inverts *each band separately*
#    to its own brightness temperature and takes the MAXIMUM across bands
#    (within wn_range) - the "Christiansen feature" technique. This
#    deliberately does NOT average over the spectral information the way (2)
#    does - it exploits it, picking out whichever band peers closest to the
#    true near-surface temperature through the semi-transparent medium. It is
#    a genuine, separate physical finding about hybrid's spectral behavior,
#    not a consistency check against two_wave: two_wave has only one band, so
#    there is no spectral peak for it to reveal, and a difference between
#    Tb_max_hybrid and T_surf_two_wave (or Tb_bol_two_wave) is *expected*,
#    not a sign of disagreement between the two modes.
#
# modelmain.max_btemp_blackbody/radiance_processor._fit_brightness_temperature
# already implement a per-band version of (3), but _fit_brightness_temperature
# passes a MinimalSim proxy object lacking T_out into max_btemp_blackbody,
# which raises AttributeError and is silently caught, always falling back to
# the crude Stefan-Boltzmann approximation - so results['brightness_temps'] in
# radiance_processor.py never actually exploits spectral information despite
# appearing to. max_brightness_temperature_series() reimplements the per-band
# inversion directly (as a vectorized Simpson's-rule lookup table, validated
# to <0.0001 K against modelmain.planck_wn_integrated's quad-based
# integration - both far more accurate for wide bands and far faster than
# optimizing each band/timestep individually).
# ---------------------------------------------------------------------------

def _planck_wn_density(wn_cm, T):
    """Spectral radiance density per cm^-1 (vectorized over wn_cm and/or T via
    broadcasting), matching the unit convention of modelmain.planck_wn_integrated
    (i.e. the per-m^-1 Planck function scaled by 100 to get per-cm^-1 density)."""
    h, c, k = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    wn_m = wn_cm * 100.0
    return (2.0 * h * c**2 * wn_m**3) / (np.exp(h * c * wn_m / (k * T)) - 1.0) * 100.0


def _band_radiance_lookup(wn_bounds, T_grid, n_sub=21):
    """Vectorized Simpson's-rule band-integrated Planck radiance for every
    band in wn_bounds (n_bands = len(wn_bounds)-1), across every T in T_grid.

    Returns
    -------
    ndarray, shape [n_bands, len(T_grid)]
    """
    wn_low = wn_bounds[:-1]
    wn_high = wn_bounds[1:]
    frac = np.linspace(0.0, 1.0, n_sub)
    sub_wn = wn_low[:, None] + frac[None, :] * (wn_high - wn_low)[:, None]  # [n_bands, n_sub]
    hgrid = (wn_high - wn_low) / (n_sub - 1)
    w = np.ones(n_sub)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    w = w[None, :] * (hgrid / 3.0)[:, None]  # [n_bands, n_sub]
    B = _planck_wn_density(sub_wn[:, :, None], T_grid[None, None, :])  # [n_bands, n_sub, n_T]
    return (w[:, :, None] * B).sum(axis=1)  # [n_bands, n_T]


def bolometric_brightness_temperature_series(sim, observer_angle=0.0):
    """
    Compute the per-output-time *directional* brightness temperature from
    DISORT's actual upward thermal radiance at the given observer angle,
    integrated (summed) across the full thermal spectrum - the fair,
    apples-to-apples generalization of the single-band
    `(radiance*pi/sigma)**0.25` calculation used for two_wave-only output, now
    extended correctly to hybrid/multi-wave.

    Each output band's radiance from calculate_radiances_from_results is
    already band-integrated (an extensive quantity spanning that band's
    width), not a spectral density - matching the convention used throughout
    this module (e.g. modelmain.planck_wn_integrated's B_bands) and confirmed
    by modelmain.Simulator.run()'s T_surf calculation, which sums
    flux_up_therm directly across bands with no bandwidth weighting. So the
    bolometric directional radiance is simply the sum across bands, with no
    extra integration weights needed.

    Unlike max_brightness_temperature_series (which deliberately picks the
    single most-transparent band to reveal the near-surface temperature),
    this quantity is directly comparable to T_surf and to a two_wave run's
    equivalent number: it's the same "assume the emission is Lambertian, so
    flux = pi*L" inversion, just applied to the true directional radiance at
    observer_angle rather than to the hemispheric flux DISORT also reports.
    It need not exactly equal T_surf even for a self-consistent run, since
    T_surf is hemispheric (angle-integrated) while this is a single
    direction - any difference reflects real (typically modest) deviation
    from Lambertian emission (limb-darkening/brightening), not an error.

    Parameters
    ----------
    sim : Simulator
        A simulator that has already run() to completion.
    observer_angle : float
        Observer emission angle in degrees from nadir (0 = straight overhead).

    Returns
    -------
    times, T_bright : ndarray
        Output time array (s) and bolometric directional brightness
        temperature (K) at each time.
    """
    from radiance_processor import calculate_radiances_from_results

    spec = calculate_radiances_from_results(sim, surface_type='smooth', observer_angles=[observer_angle])
    times = spec['times']
    sigma = 5.670374419e-8

    radiance = spec['radiance_thermal'][0]  # [n_waves, n_times] (multi-wave) or [n_times] (two_wave)
    total_radiance = np.sum(radiance, axis=0) if radiance.ndim == 2 else radiance
    T_bright = (total_radiance * np.pi / sigma) ** 0.25

    return times, T_bright


def max_brightness_temperature_series(sim, observer_angle=0.0, T_min=100.0, T_max=450.0, n_T=3501,
                                       wn_range=(400.0, 2000.0)):
    """
    Compute the per-output-time brightness temperature as actually perceived
    by an observer at the given emission angle, using DISORT's directional
    upward thermal RADIANCE (not Simulator's internal hemispheric-flux-based
    T_surf) - this is what properly captures the semi-transparency of the
    dust layer instead of assuming an opaque blackbody surface.

    For hybrid/multi-wave output: inverts each spectral band's radiance to a
    brightness temperature via an essentially-exact Simpson's-rule lookup
    table, then takes the MAXIMUM across bands (within wn_range) at each time
    - the standard "Christiansen feature" technique for recovering the least-
    biased estimate of the true near-surface temperature.

    IMPORTANT, found by direct investigation while building this: bands deep
    in the Wien tail (e.g. ~1 um, for a ~200-300K surface) are NOT numerically
    unreliable - their tiny radiance values were checked against
    modelmain.planck_wn_integrated at the actual simulated temperature range
    and are physically consistent, real (if minuscule) thermal emission. The
    "excess" brightness temperature they imply over the bulk/broadband value
    is a genuine effect: with a real temperature gradient through the
    (semi-transparent) column, Planck's function is so nonlinear at short
    wavelengths that emission there becomes dominated almost entirely by the
    single hottest layer (Jensen's-inequality-like behavior) - so an
    unrestricted max-across-bands search doesn't converge to a stable
    "brightness temperature", it trends toward the hottest node in the
    profile as you push to shorter wavelengths, limited only by how far the
    Mie file's wavelength coverage extends. That isn't what a real observer
    would measure: no real instrument has sensitivity to a ~300K surface at
    ~1e-16 W/m^2/sr/cm^-1 (1 um). wn_range restricts the search to a genuinely
    observable thermal-IR window (default 5-25 um / 400-2000 cm^-1, comfortably
    spanning the Planck peak for typical planetary surface temperatures and
    matching the spirit - if not the exact silicate-specific band - of
    modelmain.max_btemp_blackbody's hardcoded 900-1700 cm^-1 default).

    For two_wave (broadband) output: there's no spectral information to
    exploit, so this is just the standard blackbody inversion of the single
    broadband radiance value.

    Requires sim.cfg.output_radiance_mode/mie_file_out/wn_bounds_out to
    already be configured consistently with sim.cfg.thermal_evolution_mode/
    mie_file/wn_bounds - run_diurnal_comparison() does this automatically.

    Parameters
    ----------
    sim : Simulator
        A simulator that has already run() to completion.
    observer_angle : float
        Observer emission angle in degrees from nadir (0 = straight overhead).
    T_min, T_max, n_T : float, float, int
        Temperature grid used to build the inversion lookup table. Widen if
        you see values pinned at these bounds (a sign radiance fell outside
        the table's range).
    wn_range : (float, float)
        Wavenumber range (cm^-1) the max-across-bands search is restricted
        to - see note above for why this is necessary, not just a
        convenience. Widen only if you have a specific reason to trust
        radiance further into the Wien tail (e.g. genuinely modeling a
        specific real instrument's broader bandpass).

    Returns
    -------
    times, T_bright : ndarray
        Output time array (s) and brightness temperature (K) at each time.
    """
    from radiance_processor import calculate_radiances_from_results

    spec = calculate_radiances_from_results(sim, surface_type='smooth', observer_angles=[observer_angle])
    times = spec['times']
    sigma = 5.670374419e-8

    if 'wavenumbers' in spec:
        radiance = spec['radiance_thermal'][0]  # [n_waves, n_times]
        wn_bounds = np.sort(np.loadtxt(sim.cfg.wn_bounds_out))
        if len(wn_bounds) != radiance.shape[0] + 1:
            raise ValueError(
                f"{sim.cfg.wn_bounds_out} has {len(wn_bounds)} edges, expected "
                f"{radiance.shape[0] + 1} to match the {radiance.shape[0]} output bands."
            )
        wn_centers = 0.5 * (wn_bounds[:-1] + wn_bounds[1:])
        in_range = (wn_centers >= wn_range[0]) & (wn_centers <= wn_range[1])
        if not np.any(in_range):
            raise ValueError(
                f"No output bands fall within wn_range={wn_range} cm^-1; "
                f"available range is [{wn_centers.min():.1f}, {wn_centers.max():.1f}] cm^-1."
            )

        T_grid = np.linspace(T_min, T_max, n_T)
        lut = _band_radiance_lookup(wn_bounds, T_grid)  # [n_bands, n_T]

        n_bands, n_times = radiance.shape
        T_per_band = np.full((n_bands, n_times), np.nan)
        for b in np.where(in_range)[0]:
            valid = radiance[b] > 0
            T_per_band[b, valid] = np.interp(radiance[b, valid], lut[b], T_grid)
        if np.any((T_per_band[in_range] <= T_min + 1.0) | (T_per_band[in_range] >= T_max - 1.0)):
            print(f"Warning: some per-band brightness temperatures are pinned near the "
                  f"[{T_min}, {T_max}] K lookup-table bounds - widen T_min/T_max.")
        T_bright = np.nanmax(T_per_band, axis=0)
    else:
        radiance = spec['radiance_thermal'][0]
        T_bright = (radiance * np.pi / sigma) ** 0.25

    return times, T_bright
