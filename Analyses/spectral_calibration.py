"""
Lightweight tools for calibrating broadband single-scattering albedo parameters
(ssalb_vis, ssalb_therm) against known bulk radiometric properties (bond albedo,
thermal emissivity), without running a full diurnal thermal simulation.

Reuses Simulator.compute_spectral_properties (modelmain.py), which computes the
upward flux from a single isothermal RTE solve with the sun straight overhead
(mu=1, F=1) for albedo, and a single no-sun solve (mu=0, F=0) for emissivity.
No time-stepping is required: constructing a Simulator already builds the
DISORT/Hapke solver instances needed for this calculation.
"""
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import copy
import numpy as np
from scipy.optimize import brentq

from modelmain import Simulator
from config import SimulationConfig


def effective_albedo_emissivity(cfg: SimulationConfig, T=None, solver_mode=None):
    """
    Compute the isothermal, sun-overhead bond albedo and thermal emissivity of
    the regolith column described by `cfg`.

    Parameters
    ----------
    cfg : SimulationConfig
        Base configuration. A deep copy is used internally, so `cfg` itself is
        never mutated.
    T : float, optional
        Isothermal temperature (K) to evaluate at. Defaults to cfg.T_bottom.
    solver_mode : str, optional
        Overrides cfg.thermal_evolution_mode for this calculation only, e.g.
        'two_wave' (broadband, matches typical diurnal-run settings) or
        'hybrid' (spectral thermal emissivity + broadband visible albedo).
        Defaults to cfg.thermal_evolution_mode.

    Returns
    -------
    albedo : float or ndarray
        Bond albedo. Scalar for 'two_wave'/'hybrid', spectral array for
        'multi_wave'.
    emissivity : float or ndarray
        Thermal emissivity. Scalar for 'two_wave', spectral array (per
        wavenumber) for 'hybrid'/'multi_wave'.
    nwaves : int
        Number of wavelength bands represented in `emissivity` (1 if scalar).
    """
    cfg = copy.deepcopy(cfg)
    cfg.crater = False  # crater geometry is never needed for this calculation
    if solver_mode is not None:
        cfg.thermal_evolution_mode = solver_mode
    if T is not None:
        cfg.T_bottom = T
    cfg.__post_init__()

    sim = Simulator(cfg)
    T_profile = np.full(sim.grid.x_num, cfg.T_bottom)
    albedo, emissivity, _, nwaves = sim.compute_spectral_properties(
        T_profile, cfg.thermal_evolution_mode
    )
    return albedo, emissivity, nwaves


def solve_ssalb_vis(cfg: SimulationConfig, target_albedo: float,
                     bracket=(1e-6, 0.999999), T=None, solver_mode=None,
                     xtol=1e-6):
    """
    Solve for the ssalb_vis value that reproduces `target_albedo`, holding all
    other config parameters (including g_vis) fixed.

    Parameters
    ----------
    cfg : SimulationConfig
        Base configuration (not mutated).
    target_albedo : float
        Target bond albedo to match (e.g. a literature value for Bennu).
    bracket : (float, float)
        Search bracket for ssalb_vis. Must bound the root (residual must
        change sign across the bracket) or scipy raises a ValueError.
    T, solver_mode :
        Passed through to effective_albedo_emissivity.
    xtol : float
        Absolute tolerance on ssalb_vis passed to scipy.optimize.brentq.

    Returns
    -------
    ssalb_vis : float
        Best-fit value.
    """
    def residual(ssalb_vis):
        trial_cfg = copy.deepcopy(cfg)
        trial_cfg.ssalb_vis = ssalb_vis
        albedo, _, _ = effective_albedo_emissivity(trial_cfg, T=T, solver_mode=solver_mode)
        return float(albedo) - target_albedo

    r_lo, r_hi = residual(bracket[0]), residual(bracket[1])
    if r_lo * r_hi > 0:
        albedo_lo, albedo_hi = r_lo + target_albedo, r_hi + target_albedo
        raise ValueError(
            f"target_albedo={target_albedo} is not reachable by varying ssalb_vis "
            f"alone over bracket {bracket}: achievable albedo range is "
            f"[{min(albedo_lo, albedo_hi):.4f}, {max(albedo_lo, albedo_hi):.4f}]. "
            "For optically thin dust layers, the substrate reflectivity (R_base / "
            "substrate_spectrum) can dominate the bond albedo and set a floor/ceiling "
            "that ssalb_vis alone cannot cross - check R_base and dust_thickness/Et "
            "(optical depth) before assuming ssalb_vis is the free parameter to blame."
        )

    return brentq(residual, *bracket, xtol=xtol)


if __name__ == "__main__":
    # Minimal smoke test / usage example.
    cfg = SimulationConfig(use_RTE=True, RTE_solver='disort',
                            thermal_evolution_mode='two_wave')
    albedo, emissivity, nwaves = effective_albedo_emissivity(cfg)
    print(f"ssalb_vis={cfg.ssalb_vis}: bond albedo={albedo:.4f}, "
          f"emissivity={emissivity:.4f} (nwaves={nwaves})")

    # TODO: fill in a literature Bennu bond albedo value before using this.
    target_albedo = None
    if target_albedo is not None:
        best_ssalb_vis = solve_ssalb_vis(cfg, target_albedo)
        print(f"Best-fit ssalb_vis for target albedo {target_albedo}: "
              f"{best_ssalb_vis:.6f}")
