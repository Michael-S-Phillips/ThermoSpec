"""
Two-layer effective (apparent) thermal-inertia study - v2.

Rewrite of 2layer_effective_TI_analysis.py. Fits an "effective single-layer
thermal inertia" to two-layer dust-over-substrate models, for a sweep of
dust thicknesses, by inverting/fitting against a one-layer, non-RTE,
homogeneous lookup table (LUT) - reproducing the apparent-thermal-inertia
methodology of Biele et al. 2019 ("Effects of dust layers on thermal
emission from airless bodies").

Two two-layer "target" models are fit against the same LUT:

- Biele-style (build_biele_two_layer_sweep): dust treated as thermally
  opaque, pure conduction, fixed material properties taken directly from
  Biele et al. 2019 Table 1. Literal reproduction - only dust_thickness is
  swept, no grain-size/porosity dependence. Fit target: plain kinetic
  T_surf.
- RTE (build_rte_two_layer_sweep): full spectral (hybrid: multi-wave
  thermal + broadband visible) DISORT radiative transfer, using real Mie
  optical data. Phonon-only k_dust is computed from grain diameter and
  porosity via the Gundlach & Blum contact-conductivity model, so this
  branch's true parameter space is (dust thickness x grain size x
  porosity). Fit target: Tb_bol (directional bolometric brightness
  temperature - what a nadir-viewing thermal camera integrating the whole
  blackbody spectrum reports), via
  hybrid_calibration.bolometric_brightness_temperature_series.

Fitting the RTE model against a non-RTE LUT is intentional: historical
Bennu thermal-inertia fits used non-RTE models with no knowledge of dust or
its radiative properties. Comparing the Biele-style and RTE apparent-TI
curves shows what such a naive fit would report under each assumption.

Six apparent-TI extraction methods are computed per (thickness, branch)
pair, matching Biele et al. 2019: chi-squared fits to the full diurnal
curve, day-only, and night-only (the headline result), plus phase-lag,
Tmax, and Tmin (bonus diagnostics).

See the approved plan (.claude/plans/we-re-working-on-a-glistening-popcorn.md)
for full design rationale.
"""
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import copy
import json
import datetime
from dataclasses import asdict

import numpy as np
import h5py
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import CubicSpline, interp1d
from scipy.optimize import minimize_scalar

from config import SimulationConfig
from modelmain import Simulator
from hybrid_calibration import (
    set_thermal_Et_from_mie,
    set_eta_geometric_optics,
    write_thermal_wn_bounds_file,
    bolometric_brightness_temperature_series,
    max_brightness_temperature_series,
)
from spectral_calibration import solve_ssalb_vis

SIGMA = 5.670374419e-8

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
_OPTICAL_PROPS_DIR = os.path.join(_REPO_ROOT, 'Optical_props')
DEFAULT_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')

# -----------------------------------------------------------------------------
# Gundlach & Blum contact-conductivity model for phonon-only k_dust, as a
# function of grain diameter and porosity. ks=1.0, X=0.5 are provisional
# defaults per the user (still pending final values); poissons/youngs
# callers must pass explicitly - GB_POISSONS_DEFAULT/GB_YOUNGS_DEFAULT below
# are the user-supplied real values for the dust material.
# -----------------------------------------------------------------------------
_GB_F1 = 5.18e-2
_GB_F2 = 5.26
_GB_SURFENERGY_SIO = 0.020  # J/m^2
GB_POISSONS_DEFAULT = 0.269    # dimensionless, user-supplied
GB_YOUNGS_DEFAULT = 5.625e9    # Pa, user-supplied


def gundlach_blum_k_dust(sphere_diam, porosity, poissons, youngs, ks=1.0, X=0.5):
    """Phonon-only dust thermal conductivity (W/m/K) from grain diameter (m)
    and porosity, after Gundlach & Blum."""
    buffer = _GB_F1 * np.exp(_GB_F2 * (1.0 - porosity))
    GB_s = (
        (9.0 * np.pi / 4.0) * ((1.0 - poissons ** 2.0) / youngs)
        * (_GB_SURFENERGY_SIO / (sphere_diam / 2.0))
    ) ** (1.0 / 3.0) * buffer * X * ks
    return GB_s


# -----------------------------------------------------------------------------
# Shared assumptions - applied identically across the LUT, Biele-branch, and
# RTE-branch configs. Sharing P/ndays/freq_out/last_day/latitude/dec makes
# every run's t_out bit-identical (see plan: this is what fixes the old
# time-alignment bug at the root, rather than patching it).
# -----------------------------------------------------------------------------
TIMING_DEFAULTS = dict(
    P=15450.0,             # s, Bennu-like rotation period (Two_layer_Bennu_dust.ipynb)
    ndays=5,
    freq_out=96,
    last_day=True,
    latitude=np.radians(50.0),
    dec=np.radians(0.0),
    diurnal=True,
    sun=True,
)

# Bond albedo / thermal emissivity assumptions for the two NON-RTE branches
# (LUT and Biele-style), which have no spectral calibration of their own.
# Bennu-like bond albedo, matching the ssalb_vis target used in
# Two_layer_Bennu_dust.ipynb. Adjust if a different literature value is
# preferred - flagged as an assumption in the plan.
NONRTE_ALBEDO_DEFAULT = 0.025
NONRTE_EM_DEFAULT = 0.95

DUST_THICKNESS_DEFAULT = np.array([
    5.0e-6, 10.0e-6, 30.0e-6, 50.0e-6, 100.0e-6, 200.0e-6,
    500.0e-6, 0.001, 0.002, 0.005, 0.01, 0.02,
])
K_DUST_LUT_DEFAULT = np.logspace(-5, 0.5, 25)

# Available Mie tables keyed by grain diameter (m): (mie_file, radius).
# Only 1um/5um serpentine exist today - more grain sizes require generating
# additional Mie tables via the Preprocessing/ pipeline (see plan).
GRAIN_MIE_FILES = {
    1e-6: dict(
        mie_file=os.path.join(_OPTICAL_PROPS_DIR, 'serpentine_mie_200wns_1um.txt'),
        radius=0.5e-6,
    ),
    5e-6: dict(
        mie_file=os.path.join(_OPTICAL_PROPS_DIR, 'serpentine_mie_200wns_5um.txt'),
        radius=2.5e-6,
    ),
}
WN_BOUNDS_BASE = os.path.join(_OPTICAL_PROPS_DIR, 'wn_bounds_200.txt')

# Placeholder bulk properties for the RTE branch's dust layer. rho_dust is
# derived from porosity (via fill_frac) at run time using GRAIN_DENSITY_DEFAULT;
# cp_dust reuses the serpentine value calibrated in Two_layer_Bennu_dust.ipynb.
GRAIN_DENSITY_DEFAULT = 2500.0  # kg/m^3, approximate serpentine grain density - placeholder
CP_DUST_DEFAULT = 800.0         # J/kg/K

# Biele et al. 2019 Table 1 fixed dust properties - literal reproduction,
# decoupled from the RTE branch's physically-modeled phonon conductivity.
BIELE_K_DUST = 0.0025      # W/m/K (phonon + radiative, non-temperature-dependent)
BIELE_RHO_DUST = 366.0     # kg/m^3
BIELE_CP_DUST = 700.0      # J/kg/K

# Rock/substrate properties (Two_layer_Bennu_dust.ipynb baseline), shared by
# both two-layer branches.
ROCK_DEFAULTS = dict(k_rock=0.5, rho_rock=1800.0, cp_rock=800.0, rock_thickness=1.0)

# Substrate reflectivity for DISORT (RTE branch only - not used by the
# plain-conduction Biele branch).
SUBSTRATE_R_BASE = 0.025


def _apply_timing(cfg, timing=None):
    timing = timing or TIMING_DEFAULTS
    for key, val in timing.items():
        setattr(cfg, key, val)
    return cfg


def _make_common_numerics(cfg):
    """Numerics/boundary settings shared across all branches, matching the
    validated Two_layer_Bennu_dust.ipynb baseline."""
    cfg.auto_dt = False
    cfg.tsteps_day = 5000
    cfg.dtfac = 80000.0
    cfg.minsteps = 50000
    cfg.T_bottom = 260.0
    cfg.bottom_bc = 'dirichlet'
    cfg.disort_space_temp = 0.0
    cfg.enable_diurnal_convergence = False
    cfg.crater = False
    return cfg


def _refine_day_peak(t_out, T, day_mask):
    """Quintic-polynomial refinement of the daytime temperature peak (time,
    value), restricted to the (contiguous) day window."""
    idx = np.where(day_mask)[0]
    t_day = t_out[idx]
    T_day = T[idx]
    poly = np.polyfit(t_day, T_day, 5)
    p = np.poly1d(poly)
    res = minimize_scalar(lambda t: -p(t), bounds=(t_day.min(), t_day.max()), method='bounded')
    return float(p(res.x)), float(res.x)


# -----------------------------------------------------------------------------
# Timestep convergence check: repeatedly double tsteps_day (starting from
# whatever the caller's base config already has, e.g. 5000 from
# _make_common_numerics) until the fit-relevant output stops changing by
# more than tol_K, fully independent per run (no cross-run reuse or
# interpolation - deliberately simple per the user's preference, to avoid
# baking in untested assumptions about which parameters a converged
# resolution transfers across).
#
# This directly fixes a real bug found via diagnosis: the RTE branch's
# near-surface grid layer (dust_rte_max_lthick, in tau units) becomes very
# thin in physical terms once Et is large (real Mie-calibrated Et is ~18x
# the placeholder value used in prior validation), and the fixed
# tsteps_day=5000 inherited from that prior validation silently
# under-resolves it - confirmed empirically (5000 steps/day gave a 67 K
# diurnal amplitude for a materially "infinite dust" case; 100000 steps/day
# gave 157.5 K, still not converged). auto_dt was considered and rejected -
# the user's own experience is that it overshoots unreliably and cannot be
# trusted as a stand-in for actual verified convergence.
#
# Capped (max_doublings) to avoid unbounded runaway compute for pathological
# cases; a capped-out run is flagged as non-converged in its output, never
# silently trusted.
# -----------------------------------------------------------------------------

def run_with_convergence_check(build_and_run_fn, extract_fn, base_tsteps_day=5000,
                                tol_K=0.1, max_doublings=8, verbose=True, tag=''):
    """build_and_run_fn(tsteps_day) -> completed Simulator.
    extract_fn(sim) -> 1D ndarray, the fit-relevant series to check for
    convergence (e.g. sim.T_surf_out or a Tb_bol array).

    Returns a dict: sim, converged (bool), n_doublings, tsteps_day,
    max_abs_diff (K, from the final comparison)."""
    label = f' {tag}' if tag else ''
    tsteps_day = base_tsteps_day
    sim_prev = build_and_run_fn(tsteps_day)
    series_prev = extract_fn(sim_prev)
    max_abs_diff = None
    for n in range(1, max_doublings + 1):
        tsteps_day *= 2
        sim_cur = build_and_run_fn(tsteps_day)
        series_cur = extract_fn(sim_cur)
        max_abs_diff = float(np.max(np.abs(series_cur - series_prev)))
        if verbose:
            print(f"  [convergence{label}] tsteps_day={tsteps_day}: "
                  f"max_abs_diff={max_abs_diff:.4f} K")
        if max_abs_diff < tol_K:
            return {'sim': sim_cur, 'converged': True, 'n_doublings': n,
                    'tsteps_day': tsteps_day, 'max_abs_diff': max_abs_diff}
        sim_prev, series_prev = sim_cur, series_cur
    print(f"  [convergence{label}] WARNING: did not converge within "
          f"{max_doublings} doublings (final tsteps_day={tsteps_day}, "
          f"max_abs_diff={max_abs_diff:.4f} K >= {tol_K} K) - result may still "
          f"contain significant numerical error.")
    return {'sim': sim_cur, 'converged': False, 'n_doublings': max_doublings,
            'tsteps_day': tsteps_day, 'max_abs_diff': max_abs_diff}


# -----------------------------------------------------------------------------
# 1. One-layer non-RTE LUT
# -----------------------------------------------------------------------------

def build_single_layer_lut(k_dust_values=None, timing=None, rho_dust_ref=None,
                            cp_dust_ref=None, albedo=NONRTE_ALBEDO_DEFAULT,
                            em=NONRTE_EM_DEFAULT, enable_convergence_check=False,
                            convergence_tol_K=0.1, max_doublings=8, verbose=True):
    """One-layer, non-RTE homogeneous lookup table swept over k_dust - the
    common reference both two-layer branches are fit against."""
    if k_dust_values is None:
        k_dust_values = K_DUST_LUT_DEFAULT
    timing = timing or TIMING_DEFAULTS
    if rho_dust_ref is None:
        rho_dust_ref = ROCK_DEFAULTS['rho_rock']
    if cp_dust_ref is None:
        cp_dust_ref = ROCK_DEFAULTS['cp_rock']

    canonical_t_out = None
    canonical_mu_out = None
    T_surf_list = []
    convergence_list = []

    for k in k_dust_values:
        cfg = SimulationConfig()
        cfg.single_layer = True
        cfg.use_RTE = False
        cfg.k_dust_auto = False
        _apply_timing(cfg, timing)
        _make_common_numerics(cfg)
        cfg.albedo = albedo
        cfg.em = em
        cfg.k_dust = float(k)
        cfg.rho_dust = rho_dust_ref
        cfg.cp_dust = cp_dust_ref
        cfg.dust_thickness = 1.0  # deep column -> semi-infinite-slab approximation
        cfg.__post_init__()

        if verbose:
            print(f"[LUT] k_dust={k:.4e} W/m/K")

        def _build_and_run(tsteps_day, _cfg=cfg):
            c = copy.deepcopy(_cfg)
            c.tsteps_day = tsteps_day
            c.__post_init__()
            s = Simulator(c)
            s.run()
            return s

        if enable_convergence_check:
            conv = run_with_convergence_check(
                _build_and_run, lambda s: s.T_surf_out, base_tsteps_day=cfg.tsteps_day,
                tol_K=convergence_tol_K, max_doublings=max_doublings, verbose=verbose,
                tag=f'LUT k={k:.2e}')
            sim = conv['sim']
            convergence_list.append({ck: cv for ck, cv in conv.items() if ck != 'sim'})
        else:
            sim = _build_and_run(cfg.tsteps_day)
            convergence_list.append(None)

        if canonical_t_out is None:
            canonical_t_out = sim.t_out.copy()
            canonical_mu_out = sim.mu_out.copy()
        else:
            assert np.allclose(sim.t_out, canonical_t_out), (
                "LUT run produced a different t_out than earlier LUT runs - "
                "shared timing config invariant violated."
            )
        T_surf_list.append(sim.T_surf_out.copy())

    T_surf_arr = np.array(T_surf_list)  # [n_k, n_t]
    log_k = np.log10(k_dust_values)
    lut_k_interp = CubicSpline(log_k, T_surf_arr, axis=0)
    day_mask = canonical_mu_out > 0.001

    n_k = len(k_dust_values)
    max_T = np.zeros(n_k)
    min_T = np.zeros(n_k)
    max_time = np.zeros(n_k)
    for i in range(n_k):
        max_T[i], max_time[i] = _refine_day_peak(canonical_t_out, T_surf_arr[i], day_mask)
        night_vals = T_surf_arr[i][~day_mask]
        min_T[i] = np.min(night_vals) if night_vals.size else np.min(T_surf_arr[i])

    t_noon = float(canonical_t_out[np.argmax(canonical_mu_out)])
    period = float(timing['P'])
    lag_deg = (max_time - t_noon) * (360.0 / period)

    return {
        'k_dust_values': np.asarray(k_dust_values), 'log_k': log_k,
        't_out': canonical_t_out, 'mu_out': canonical_mu_out, 'day_mask': day_mask,
        'T_surf': T_surf_arr, 'lut_k_interp': lut_k_interp,
        'max_T': max_T, 'min_T': min_T, 'max_time': max_time,
        't_noon': t_noon, 'period': period, 'lag_deg': lag_deg,
        'rho_dust_ref': rho_dust_ref, 'cp_dust_ref': cp_dust_ref,
        'albedo': albedo, 'em': em, 'convergence': convergence_list,
    }


# -----------------------------------------------------------------------------
# 2. Biele-style two-layer branch (non-RTE, opaque dust, fixed properties)
# -----------------------------------------------------------------------------

def build_biele_two_layer_sweep(dust_thickness_values=None, timing=None,
                                 rock=None, albedo=NONRTE_ALBEDO_DEFAULT,
                                 em=NONRTE_EM_DEFAULT, canonical_t_out=None,
                                 enable_convergence_check=False,
                                 convergence_tol_K=0.1, max_doublings=8,
                                 verbose=True):
    """Literal reproduction of Biele et al. 2019's two-layer opaque-dust
    conduction model: fixed dust material properties (their Table 1 values),
    only dust_thickness swept. No grain-size/porosity dependence, and no
    core-model changes needed - grid.py/modelmain.py already build a genuine
    step-function two-material FD grid and boundary solve for
    single_layer=False, use_RTE=False."""
    if dust_thickness_values is None:
        dust_thickness_values = DUST_THICKNESS_DEFAULT
    timing = timing or TIMING_DEFAULTS
    rock = rock or ROCK_DEFAULTS

    results = {}
    for thickness in dust_thickness_values:
        cfg = SimulationConfig()
        cfg.single_layer = False
        cfg.use_RTE = False
        cfg.k_dust_auto = False
        _apply_timing(cfg, timing)
        _make_common_numerics(cfg)
        cfg.albedo = albedo
        cfg.em = em
        cfg.k_dust = BIELE_K_DUST
        cfg.rho_dust = BIELE_RHO_DUST
        cfg.cp_dust = BIELE_CP_DUST
        for key, val in rock.items():
            setattr(cfg, key, val)
        cfg.dust_thickness = float(thickness)
        cfg.__post_init__()

        if verbose:
            print(f"[Biele branch] dust_thickness={thickness:.3e} m")

        def _build_and_run(tsteps_day, _cfg=cfg):
            c = copy.deepcopy(_cfg)
            c.tsteps_day = tsteps_day
            c.__post_init__()
            s = Simulator(c)
            s.run()
            return s

        if enable_convergence_check:
            conv = run_with_convergence_check(
                _build_and_run, lambda s: s.T_surf_out, base_tsteps_day=cfg.tsteps_day,
                tol_K=convergence_tol_K, max_doublings=max_doublings, verbose=verbose,
                tag=f'Biele d={thickness:.2e}')
            sim = conv['sim']
            convergence = {ck: cv for ck, cv in conv.items() if ck != 'sim'}
        else:
            sim = _build_and_run(cfg.tsteps_day)
            convergence = None

        if canonical_t_out is not None:
            assert np.allclose(sim.t_out, canonical_t_out), (
                "Biele-branch run produced a different t_out than the LUT - "
                "shared timing config invariant violated."
            )

        results[float(thickness)] = {
            'config': sim.cfg, 't_out': sim.t_out.copy(), 'mu_out': sim.mu_out.copy(),
            'T_surf': sim.T_surf_out.copy(), 'convergence': convergence,
        }
    return results


# -----------------------------------------------------------------------------
# 3. RTE (hybrid) two-layer branch
# -----------------------------------------------------------------------------

def build_rte_two_layer_sweep(dust_thickness_values=None, grain_diameters=None,
                               porosities=None, poissons=None, youngs=None,
                               ks=1.0, X=0.5, timing=None, rock=None,
                               target_bond_albedo=NONRTE_ALBEDO_DEFAULT,
                               observer_angle=0.0, compute_tb_max=True,
                               canonical_t_out=None, output_dir=DEFAULT_OUTPUT_DIR,
                               enable_convergence_check=False,
                               convergence_tol_K=0.1, max_doublings=8,
                               verbose=True):
    """RTE (hybrid/multi-wave) two-layer sweep over dust thickness x grain
    diameter x porosity. Phonon-only k_dust from gundlach_blum_k_dust;
    Et/eta calibrated per (grain, porosity) combo from the matching Mie
    file; ssalb_vis calibrated once per combo against the thickest ("most
    opaque") dust thickness in the sweep, since thin dust has little
    leverage over bond albedo (substrate-reflectivity-dominated). Fit
    target: Tb_bol (directional bolometric brightness temperature)."""
    if dust_thickness_values is None:
        dust_thickness_values = DUST_THICKNESS_DEFAULT
    if grain_diameters is None:
        grain_diameters = sorted(GRAIN_MIE_FILES.keys())
    if porosities is None:
        porosities = [1.0 - 0.37]  # matches Two_layer_Bennu_dust.ipynb's fill_frac=0.37 baseline
    if poissons is None or youngs is None:
        raise ValueError(
            "poissons (Poisson's ratio) and youngs (Young's modulus, Pa) are "
            "required for the Gundlach & Blum phonon-conductivity model and "
            "have not yet been finalized by the user (see plan 'Open "
            "dependencies'). Pass provisional values explicitly to run a "
            "Phase-A smoke test."
        )
    timing = timing or TIMING_DEFAULTS
    rock = rock or ROCK_DEFAULTS
    os.makedirs(output_dir, exist_ok=True)
    thickness_max = float(np.max(dust_thickness_values))

    results = {}
    for grain_diam in grain_diameters:
        if grain_diam not in GRAIN_MIE_FILES:
            raise ValueError(
                f"No Mie table available for grain diameter {grain_diam} m. "
                f"Available: {sorted(GRAIN_MIE_FILES.keys())}. Generate more "
                f"via the Preprocessing/ pipeline (see plan)."
            )
        mie_info = GRAIN_MIE_FILES[grain_diam]
        mie_file = mie_info['mie_file']
        radius = mie_info['radius']

        for porosity in porosities:
            fill_frac = 1.0 - porosity
            phonon_k_dust = gundlach_blum_k_dust(grain_diam, porosity, poissons, youngs, ks=ks, X=X)
            combo_tag = f"grain{grain_diam * 1e6:.1f}um_poro{porosity:.2f}"

            base_cfg = SimulationConfig()
            base_cfg.single_layer = False
            base_cfg.use_RTE = True
            base_cfg.RTE_solver = 'disort'
            base_cfg.thermal_evolution_mode = 'hybrid'
            base_cfg.output_radiance_mode = 'hybrid'
            _apply_timing(base_cfg, timing)
            _make_common_numerics(base_cfg)
            for key, val in rock.items():
                setattr(base_cfg, key, val)
            base_cfg.R_base = SUBSTRATE_R_BASE
            base_cfg.mie_file = mie_file
            base_cfg.radius = radius
            base_cfg.fill_frac = fill_frac
            base_cfg.k_dust = float(phonon_k_dust)
            base_cfg.rho_dust = GRAIN_DENSITY_DEFAULT * fill_frac
            base_cfg.cp_dust = CP_DUST_DEFAULT
            base_cfg.k_dust_auto = False  # phonon-only; RTE solves radiative transport explicitly
            base_cfg.g_vis = 0.75

            set_thermal_Et_from_mie(base_cfg, mie_file=mie_file, fill_frac=fill_frac, radius=radius)
            set_eta_geometric_optics(base_cfg, fill_frac=fill_frac, radius=radius)

            wn_bounds_path = os.path.join(output_dir, f"wn_bounds_thermal_{combo_tag}.txt")
            write_thermal_wn_bounds_file(WN_BOUNDS_BASE, mie_file, base_cfg.hybrid_wavelength_cutoff, wn_bounds_path)
            base_cfg.wn_bounds = wn_bounds_path
            base_cfg.mie_file_out = mie_file
            base_cfg.wn_bounds_out = wn_bounds_path
            base_cfg.__post_init__()

            thick_cfg = copy.deepcopy(base_cfg)
            thick_cfg.dust_thickness = thickness_max
            thick_cfg.__post_init__()
            if verbose:
                print(f"[RTE branch {combo_tag}] calibrating ssalb_vis against "
                      f"dust_thickness={thickness_max:.3e} m (target bond albedo="
                      f"{target_bond_albedo}) ...")
            best_ssalb_vis = solve_ssalb_vis(thick_cfg, target_bond_albedo)
            base_cfg.ssalb_vis = best_ssalb_vis
            base_cfg.__post_init__()

            combo_results = {}
            for thickness in dust_thickness_values:
                cfg = copy.deepcopy(base_cfg)
                cfg.dust_thickness = float(thickness)
                cfg.__post_init__()
                if verbose:
                    print(f"[RTE branch {combo_tag}] dust_thickness={thickness:.3e} m")

                def _build_and_run(tsteps_day, _cfg=cfg):
                    c = copy.deepcopy(_cfg)
                    c.tsteps_day = tsteps_day
                    c.__post_init__()
                    s = Simulator(c)
                    s.run()
                    return s

                def _extract_tb_bol(s, _angle=observer_angle):
                    _, tb = bolometric_brightness_temperature_series(s, observer_angle=_angle)
                    return tb

                if enable_convergence_check:
                    conv = run_with_convergence_check(
                        _build_and_run, _extract_tb_bol, base_tsteps_day=cfg.tsteps_day,
                        tol_K=convergence_tol_K, max_doublings=max_doublings, verbose=verbose,
                        tag=f'RTE {combo_tag} d={thickness:.2e}')
                    sim = conv['sim']
                    convergence = {ck: cv for ck, cv in conv.items() if ck != 'sim'}
                else:
                    sim = _build_and_run(cfg.tsteps_day)
                    convergence = None

                if canonical_t_out is not None:
                    assert np.allclose(sim.t_out, canonical_t_out), (
                        "RTE-branch run produced a different t_out than the "
                        "LUT - shared timing config invariant violated."
                    )

                _, Tb_bol = bolometric_brightness_temperature_series(sim, observer_angle=observer_angle)
                Tb_max = None
                if compute_tb_max:
                    _, Tb_max = max_brightness_temperature_series(sim, observer_angle=observer_angle)

                combo_results[float(thickness)] = {
                    'config': sim.cfg, 't_out': sim.t_out.copy(), 'mu_out': sim.mu_out.copy(),
                    'T_surf': sim.T_surf_out.copy(), 'Tb_bol': Tb_bol, 'Tb_max': Tb_max,
                    'convergence': convergence,
                }

            results[(grain_diam, porosity)] = {
                'ssalb_vis': best_ssalb_vis, 'Et': base_cfg.Et, 'eta': base_cfg.eta,
                'phonon_k_dust': phonon_k_dust, 'rho_dust': base_cfg.rho_dust,
                'mie_file': mie_file, 'wn_bounds': wn_bounds_path,
                'thicknesses': combo_results,
            }
    return results


# -----------------------------------------------------------------------------
# 4. Fitting methods
# -----------------------------------------------------------------------------

def _global_then_local_argmin(objective, lo, hi, n_coarse=300):
    """Coarse log-spaced scan (cheap - just spline evals) to bracket the
    minimum, then local refinement - avoids biasing the fit toward a
    possibly-inaccurate bracketing heuristic."""
    grid = np.linspace(lo, hi, n_coarse)
    vals = np.array([objective(x) for x in grid])
    i0 = int(np.argmin(vals))
    hit_boundary = (i0 == 0) or (i0 == n_coarse - 1)
    lo_ref = grid[max(i0 - 1, 0)]
    hi_ref = grid[min(i0 + 1, n_coarse - 1)]
    if lo_ref >= hi_ref:
        return float(grid[i0]), float(vals[i0]), hit_boundary
    res = minimize_scalar(objective, bounds=(lo_ref, hi_ref), method='bounded')
    return float(res.x), float(res.fun), hit_boundary


def chi2_fit(modelT, lut, mask=None):
    """Fit log10(k_dust) minimizing MSE against the LUT's spline, optionally
    restricted to a boolean time mask (day/night)."""
    log_k = lut['log_k']
    if mask is None:
        mask = np.ones_like(modelT, dtype=bool)
    n_pts = int(np.sum(mask))

    def objective(logk):
        lut_curve = lut['lut_k_interp'](logk)
        diff = modelT[mask] - lut_curve[mask]
        return float(np.mean(diff ** 2))

    best_logk, mse, hit_boundary = _global_then_local_argmin(objective, log_k.min(), log_k.max())
    best_k = 10 ** best_logk
    TI = float(np.sqrt(best_k * lut['rho_dust_ref'] * lut['cp_dust_ref']))
    return {'best_k_dust': best_k, 'best_TI': TI, 'mse': mse, 'n_points': n_pts,
            'hit_boundary': hit_boundary}


def _k_from_lut_extremum(value, lut_values, log_k):
    interp = interp1d(lut_values, log_k, fill_value='extrapolate', kind='quadratic')
    return 10 ** float(interp(value))


def tmax_fit(modelT, lut):
    T_peak, t_peak = _refine_day_peak(lut['t_out'], modelT, lut['day_mask'])
    best_k = _k_from_lut_extremum(T_peak, lut['max_T'], lut['log_k'])
    TI = float(np.sqrt(best_k * lut['rho_dust_ref'] * lut['cp_dust_ref']))
    return {'T_peak': T_peak, 't_peak': t_peak, 'best_k_dust': best_k, 'best_TI': TI}


def tmin_fit(modelT, lut):
    night_mask = ~lut['day_mask']
    night_vals = modelT[night_mask] if np.any(night_mask) else modelT
    T_min = float(np.min(night_vals))
    best_k = _k_from_lut_extremum(T_min, lut['min_T'], lut['log_k'])
    TI = float(np.sqrt(best_k * lut['rho_dust_ref'] * lut['cp_dust_ref']))
    return {'T_min': T_min, 'best_k_dust': best_k, 'best_TI': TI}


def phase_lag_fit(modelT, lut):
    _, t_peak = _refine_day_peak(lut['t_out'], modelT, lut['day_mask'])
    lag_deg = (t_peak - lut['t_noon']) * (360.0 / lut['period'])
    best_k = _k_from_lut_extremum(lag_deg, lut['lag_deg'], lut['log_k'])
    TI = float(np.sqrt(best_k * lut['rho_dust_ref'] * lut['cp_dust_ref']))
    return {'lag_deg': lag_deg, 'best_k_dust': best_k, 'best_TI': TI}


def fit_all_methods(modelT, lut):
    """All six Biele et al. 2019 apparent-TI methods. 'full'/'day'/'night'
    (chi-squared) are the headline result; 'tmax'/'tmin'/'phase_lag' are
    supplementary diagnostics."""
    day_mask = lut['day_mask']
    night_mask = ~day_mask
    return {
        'full': chi2_fit(modelT, lut, mask=None),
        'day': chi2_fit(modelT, lut, mask=day_mask),
        'night': chi2_fit(modelT, lut, mask=night_mask),
        'tmax': tmax_fit(modelT, lut),
        'tmin': tmin_fit(modelT, lut),
        'phase_lag': phase_lag_fit(modelT, lut),
    }


def attach_fits(lut, biele_results, rte_results):
    """Compute and attach fit results (res['fits']) in place, once, so
    save_results/plot_summary only ever read them."""
    for res in biele_results.values():
        res['fits'] = fit_all_methods(res['T_surf'], lut)
    for combo in rte_results.values():
        for res in combo['thicknesses'].values():
            res['fits'] = fit_all_methods(res['Tb_bol'], lut)


# -----------------------------------------------------------------------------
# 5. Output persistence (HDF5) and diagnostic plots
# -----------------------------------------------------------------------------

def _h5_write_config(group, cfg):
    group.attrs['config_yaml'] = yaml.dump(asdict(cfg), default_flow_style=False)


def _h5_write_fits(group, fits):
    for method, res in fits.items():
        g = group.create_group(method)
        for key, val in res.items():
            g.attrs[key] = val


def _h5_write_convergence(group, convergence):
    """convergence: None (check disabled) or the dict from
    run_with_convergence_check (minus 'sim')."""
    if convergence is None:
        group.attrs['enabled'] = False
        return
    group.attrs['enabled'] = True
    for key, val in convergence.items():
        group.attrs[key] = val


def save_results(output_path, lut, biele_results, rte_results, meta=None):
    """Structured HDF5 output: /meta, /lut, /two_layer_biele, /two_layer_rte,
    /summary. Requires attach_fits() to have already been called."""
    meta = meta or {}
    with h5py.File(output_path, 'w') as f:
        m = f.create_group('meta')
        for key, val in meta.items():
            if val is None:
                continue
            try:
                m.attrs[key] = val
            except TypeError:
                m.attrs[key] = json.dumps(val)

        g_lut = f.create_group('lut')
        g_lut.create_dataset('k_dust_values', data=lut['k_dust_values'])
        g_lut.create_dataset('t_out', data=lut['t_out'])
        g_lut.create_dataset('mu_out', data=lut['mu_out'])
        g_lut.create_dataset('day_mask', data=lut['day_mask'])
        g_lut.create_dataset('T_surf', data=lut['T_surf'])
        g_lut.create_dataset('max_T', data=lut['max_T'])
        g_lut.create_dataset('min_T', data=lut['min_T'])
        g_lut.create_dataset('max_time', data=lut['max_time'])
        g_lut.create_dataset('lag_deg', data=lut['lag_deg'])
        g_lut.attrs['rho_dust_ref'] = lut['rho_dust_ref']
        g_lut.attrs['cp_dust_ref'] = lut['cp_dust_ref']
        g_lut.attrs['albedo'] = lut['albedo']
        g_lut.attrs['em'] = lut['em']
        g_lut.attrs['t_noon'] = lut['t_noon']
        g_lut.attrs['period'] = lut['period']
        g_lut_conv = g_lut.create_group('convergence')
        for i, conv in enumerate(lut['convergence']):
            _h5_write_convergence(g_lut_conv.create_group(f'k_{i:03d}'), conv)

        g_biele = f.create_group('two_layer_biele')
        for i, thickness in enumerate(sorted(biele_results.keys())):
            res = biele_results[thickness]
            g = g_biele.create_group(f'thickness_{i:03d}')
            g.attrs['dust_thickness'] = thickness
            _h5_write_config(g, res['config'])
            g.create_dataset('T_surf', data=res['T_surf'])
            _h5_write_fits(g.create_group('fits'), res['fits'])
            _h5_write_convergence(g.create_group('convergence'), res['convergence'])

        g_rte = f.create_group('two_layer_rte')
        for (grain, porosity), combo in rte_results.items():
            combo_tag = f"grain{grain * 1e6:.1f}um_poro{porosity:.2f}"
            g_combo = g_rte.create_group(combo_tag)
            g_combo.attrs['grain_diameter'] = grain
            g_combo.attrs['porosity'] = porosity
            g_combo.attrs['ssalb_vis'] = combo['ssalb_vis']
            g_combo.attrs['Et'] = combo['Et']
            g_combo.attrs['eta'] = combo['eta']
            g_combo.attrs['phonon_k_dust'] = combo['phonon_k_dust']
            g_combo.attrs['rho_dust'] = combo['rho_dust']
            g_combo.attrs['mie_file'] = combo['mie_file']
            g_combo.attrs['wn_bounds'] = combo['wn_bounds']
            for i, thickness in enumerate(sorted(combo['thicknesses'].keys())):
                res = combo['thicknesses'][thickness]
                g = g_combo.create_group(f'thickness_{i:03d}')
                g.attrs['dust_thickness'] = thickness
                _h5_write_config(g, res['config'])
                g.create_dataset('T_surf', data=res['T_surf'])
                g.create_dataset('Tb_bol', data=res['Tb_bol'])
                if res['Tb_max'] is not None:
                    g.create_dataset('Tb_max', data=res['Tb_max'])
                _h5_write_fits(g.create_group('fits'), res['fits'])
                _h5_write_convergence(g.create_group('convergence'), res['convergence'])

        g_sum = f.create_group('summary')
        biele_thick = np.array(sorted(biele_results.keys()))
        g_sum.create_dataset('biele_dust_thickness', data=biele_thick)
        for method in ('full', 'day', 'night'):
            g_sum.create_dataset(
                f'biele_TI_{method}',
                data=np.array([biele_results[t]['fits'][method]['best_TI'] for t in biele_thick]))
        for (grain, porosity), combo in rte_results.items():
            combo_tag = f"grain{grain * 1e6:.1f}um_poro{porosity:.2f}"
            rte_thick = np.array(sorted(combo['thicknesses'].keys()))
            g_sum.create_dataset(f'rte_{combo_tag}_dust_thickness', data=rte_thick)
            for method in ('full', 'day', 'night'):
                g_sum.create_dataset(
                    f'rte_{combo_tag}_TI_{method}',
                    data=np.array([combo['thicknesses'][t]['fits'][method]['best_TI'] for t in rte_thick]))
    print(f"Saved results to {output_path}")


def _plot_thickness_grid(pdf, items, get_curve_key, title_prefix):
    fig = plt.figure(figsize=(15, 10))
    for i, (thickness, res) in enumerate(items):
        ax = plt.subplot(3, 4, (i % 12) + 1)
        best_k = res['fits']['full']['best_k_dust']
        lut_curve = res['_lut']['lut_k_interp'](np.log10(best_k))
        model_curve = res[get_curve_key]
        amplitude = float(np.max(model_curve) - np.min(model_curve))
        mse = res['fits']['full']['mse']
        rel_rms_pct = 100.0 * np.sqrt(mse) / amplitude if amplitude > 0 else float('nan')
        ax.plot(res['_lut']['t_out'] / 3600.0, model_curve, 'b-', label=title_prefix)
        ax.plot(res['_lut']['t_out'] / 3600.0, lut_curve, 'r--', label=f'Single-layer fit (k={best_k:.2e})')
        ax.set_title(f'd={thickness:.2e} m, TI={res["fits"]["full"]["best_TI"]:.1f}\n'
                     f'full-curve rel. RMS={rel_rms_pct:.1f}%', fontsize=8)
        if i % 12 == 0:
            ax.legend(fontsize=7)
        if (i + 1) % 12 == 0 or i == len(items) - 1:
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            if i < len(items) - 1:
                fig = plt.figure(figsize=(15, 10))


def _relative_rms_pct(thicknesses, curve_key, results_by_thickness):
    """sqrt(MSE)/diurnal-amplitude, as a percentage - large values mean the
    single-parameter homogeneous fit is a poor match to the curve shape,
    not just a surprising TI number."""
    out = []
    for t in thicknesses:
        res = results_by_thickness[t]
        curve = res[curve_key]
        amp = float(np.max(curve) - np.min(curve))
        mse = res['fits']['full']['mse']
        out.append(100.0 * np.sqrt(mse) / amp if amp > 0 else np.nan)
    return np.array(out)


def plot_summary(output_pdf, lut, biele_results, rte_results):
    """Diagnostic PDF: headline TI-vs-thickness summary (with true
    pure-material TI reference lines), a companion fit-quality (relative
    RMS) panel, then per-thickness two-layer-vs-best-fit overlays (with
    relative RMS annotated) for each branch. Requires attach_fits() to have
    already been called."""
    rock_TI = float(np.sqrt(ROCK_DEFAULTS['k_rock'] * ROCK_DEFAULTS['rho_rock'] * ROCK_DEFAULTS['cp_rock']))
    biele_dust_TI = float(np.sqrt(BIELE_K_DUST * BIELE_RHO_DUST * BIELE_CP_DUST))
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(rte_results), 1)))

    with PdfPages(output_pdf) as pdf:
        # --- Headline TI-vs-thickness summary, with true-TI reference lines ---
        fig, ax = plt.subplots(figsize=(8, 6))
        biele_thick = np.array(sorted(biele_results.keys()))
        biele_TI = np.array([biele_results[t]['fits']['full']['best_TI'] for t in biele_thick])
        ax.semilogx(biele_thick, biele_TI, 'o-', label='Biele-style (non-RTE)', color='tab:blue')
        ax.axhline(rock_TI, color='gray', linestyle=':', linewidth=1,
                    label=f'Rock TI ({rock_TI:.0f})')
        ax.axhline(biele_dust_TI, color='tab:blue', linestyle=':', linewidth=1,
                    label=f'Biele pure-dust TI ({biele_dust_TI:.1f})')
        for idx, ((grain, porosity), combo) in enumerate(rte_results.items()):
            rte_thick = np.array(sorted(combo['thicknesses'].keys()))
            rte_TI = np.array([combo['thicknesses'][t]['fits']['full']['best_TI'] for t in rte_thick])
            label = f'RTE (grain={grain * 1e6:.1f}um, poro={porosity:.2f})'
            color = colors[idx]
            ax.semilogx(rte_thick, rte_TI, 's--', label=label, color=color)
            rte_dust_TI = float(np.sqrt(combo['phonon_k_dust'] * combo['rho_dust'] * CP_DUST_DEFAULT))
            ax.axhline(rte_dust_TI, color=color, linestyle=':', linewidth=1,
                        label=f'{label} pure-dust TI ({rte_dust_TI:.1f})')
        ax.set_xlabel('Dust thickness (m)')
        ax.set_ylabel(r'Effective single-layer TI (J m$^{-2}$ K$^{-1}$ s$^{-1/2}$)')
        ax.set_title('Full-curve chi-squared fit (dotted: true pure-material TI)')
        ax.legend(fontsize=6, loc='best')
        ax.grid(True, which='both', alpha=0.3)
        pdf.savefig(fig)
        plt.close(fig)

        # --- Fit-quality (relative RMS) companion plot ---
        fig, ax = plt.subplots(figsize=(8, 6))
        biele_rel_rms = _relative_rms_pct(biele_thick, 'T_surf', biele_results)
        ax.semilogx(biele_thick, biele_rel_rms, 'o-', label='Biele-style (non-RTE)', color='tab:blue')
        for idx, ((grain, porosity), combo) in enumerate(rte_results.items()):
            rte_thick = np.array(sorted(combo['thicknesses'].keys()))
            rte_rel_rms = _relative_rms_pct(rte_thick, 'Tb_bol', combo['thicknesses'])
            label = f'RTE (grain={grain * 1e6:.1f}um, poro={porosity:.2f})'
            ax.semilogx(rte_thick, rte_rel_rms, 's--', label=label, color=colors[idx])
        ax.axhline(10.0, color='red', linestyle=':', linewidth=1, label='10% (caution threshold)')
        ax.set_xlabel('Dust thickness (m)')
        ax.set_ylabel('Full-curve fit relative RMS residual (%)')
        ax.set_title('Fit quality - large values mean the single-parameter homogeneous\n'
                      'model is a poor match to the curve shape, not just a surprising TI')
        ax.legend(fontsize=7, loc='best')
        ax.grid(True, which='both', alpha=0.3)
        pdf.savefig(fig)
        plt.close(fig)

        biele_items = []
        for thickness in sorted(biele_results.keys()):
            res = dict(biele_results[thickness])
            res['_lut'] = lut
            biele_items.append((thickness, res))
        _plot_thickness_grid(pdf, biele_items, 'T_surf', 'Two-layer (Biele)')

        for (grain, porosity), combo in rte_results.items():
            rte_items = []
            for thickness in sorted(combo['thicknesses'].keys()):
                res = dict(combo['thicknesses'][thickness])
                res['_lut'] = lut
                rte_items.append((thickness, res))
            _plot_thickness_grid(pdf, rte_items, 'Tb_bol',
                                  f'Two-layer RTE (grain={grain * 1e6:.1f}um, poro={porosity:.2f})')
    print(f"Saved diagnostic plots to {output_pdf}")


# -----------------------------------------------------------------------------
# 6. Driver
# -----------------------------------------------------------------------------

def run_full_analysis(dust_thickness_values=None, k_dust_lut_values=None,
                       grain_diameters=None, porosities=None,
                       poissons=None, youngs=None, ks=1.0, X=0.5,
                       target_bond_albedo=NONRTE_ALBEDO_DEFAULT, observer_angle=0.0,
                       compute_tb_max=True, output_dir=DEFAULT_OUTPUT_DIR,
                       enable_convergence_check=False, convergence_tol_K=0.1,
                       max_doublings=8, run_tag=None):
    os.makedirs(output_dir, exist_ok=True)
    if run_tag is None:
        run_tag = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

    timing = TIMING_DEFAULTS
    conv_kwargs = dict(enable_convergence_check=enable_convergence_check,
                        convergence_tol_K=convergence_tol_K, max_doublings=max_doublings)

    print("=== Building one-layer non-RTE LUT ===")
    lut = build_single_layer_lut(k_dust_lut_values, timing, **conv_kwargs)

    print("=== Building Biele-style two-layer (non-RTE) sweep ===")
    biele_results = build_biele_two_layer_sweep(
        dust_thickness_values, timing, canonical_t_out=lut['t_out'], **conv_kwargs)

    print("=== Building RTE (hybrid) two-layer sweep ===")
    rte_results = build_rte_two_layer_sweep(
        dust_thickness_values, grain_diameters, porosities, poissons, youngs,
        ks=ks, X=X, timing=timing, target_bond_albedo=target_bond_albedo,
        observer_angle=observer_angle, compute_tb_max=compute_tb_max,
        canonical_t_out=lut['t_out'], output_dir=output_dir, **conv_kwargs)

    print("=== Fitting all six apparent-TI methods ===")
    attach_fits(lut, biele_results, rte_results)

    h5_path = os.path.join(output_dir, f'2layer_effective_TI_v2_{run_tag}.h5')
    pdf_path = os.path.join(output_dir, f'2layer_effective_TI_v2_{run_tag}.pdf')

    meta = dict(
        script_version='v2', run_tag=run_tag, P=timing['P'], ndays=timing['ndays'],
        freq_out=timing['freq_out'], latitude=timing['latitude'], dec=timing['dec'],
        target_bond_albedo=target_bond_albedo, observer_angle=observer_angle,
        ks=ks, X=X, poissons=poissons, youngs=youngs,
    )
    save_results(h5_path, lut, biele_results, rte_results, meta=meta)
    plot_summary(pdf_path, lut, biele_results, rte_results)
    return lut, biele_results, rte_results


if __name__ == "__main__":
    # Phase-A smoke test: reduced sweep (2 thicknesses spanning thin<->thick,
    # one grain size), to validate the pipeline end-to-end. poissons/youngs
    # are the user-supplied real values (GB_POISSONS_DEFAULT/GB_YOUNGS_DEFAULT);
    # ks/X remain provisional (1.0/0.5) pending final values.
    #
    # enable_convergence_check=True is load-bearing, not optional, for the
    # RTE branch: the real Mie-calibrated Et (~130714/m) makes the RTE
    # branch's near-surface grid layer (dust_rte_max_lthick, tau units) only
    # ~0.4-0.8 microns physically, and the fixed tsteps_day inherited from
    # earlier (smaller-Et) validation silently under-resolves it - confirmed
    # via diagnosis (5000 steps/day gave a 67 K diurnal amplitude for an
    # "infinite dust" case that should show ~150-200+ K; still climbing at
    # 100000 steps/day). Only 2 thicknesses here because full RTE-branch
    # convergence checking is expensive (each doubling re-runs the full
    # diurnal cycle) - this is a real, accepted compute-vs-correctness
    # tradeoff, not a shortcut (see plan).
    smoke_thicknesses = np.array([10.0e-6, 0.02])

    lut, biele_results, rte_results = run_full_analysis(
        dust_thickness_values=smoke_thicknesses,
        k_dust_lut_values=K_DUST_LUT_DEFAULT,
        grain_diameters=[1e-6],
        porosities=[1.0 - 0.37],
        poissons=GB_POISSONS_DEFAULT,
        youngs=GB_YOUNGS_DEFAULT,
        enable_convergence_check=True,
        convergence_tol_K=0.1,
        max_doublings=8,
        run_tag='phaseA_smoketest_convcheck',
    )

    # --- Full production sweep (uncomment once Poisson's ratio / Young's
    # modulus are finalized and additional grain-size Mie tables exist): ---
    # lut, biele_results, rte_results = run_full_analysis(
    #     dust_thickness_values=DUST_THICKNESS_DEFAULT,
    #     k_dust_lut_values=K_DUST_LUT_DEFAULT,
    #     grain_diameters=sorted(GRAIN_MIE_FILES.keys()),
    #     porosities=[1.0 - 0.37],
    #     poissons=<final value>,
    #     youngs=<final value>,
    #     run_tag='production',
    # )
