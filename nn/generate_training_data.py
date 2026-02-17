"""
Generate training data for the DISORT surrogate model.

Runs thermal simulations with DISORT recording enabled to collect
(T_profile, mu, F, optical_props) -> (source_term, flux_up) pairs.

Covers three target regimes:
    - Lunar: P=2360592s, R=1AU, varied thermal inertia
    - Asteroid (Bennu-like): P=15450s, R=0.9-1.4AU
    - Lab samples: P=60-3600s, controlled T_bottom

Usage:
    python -m nn.generate_training_data --output_dir training_data
    python -m nn.generate_training_data --regime lunar --output_dir training_data/lunar
    python -m nn.generate_training_data --regime all --n_sims 200 --output_dir training_data
"""

import argparse
import logging
import sys
import time
from pathlib import Path
from itertools import product

import numpy as np
import h5py

sys.path.insert(0, str(Path(__file__).parent.parent))

from config import SimulationConfig
from modelmain import Simulator

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
#  Parameter grids for each regime
# --------------------------------------------------------------------------- #

def _lunar_params(n: int) -> list:
    """Lunar regime: slow rotator, 1 AU, varied thermal inertia and latitude."""
    rng = np.random.default_rng(42)
    params = []
    for _ in range(n):
        params.append({
            'P': 2360592.0,        # Lunar synodic day (s)
            'R': 1.0,              # AU
            'k_dust': rng.uniform(5e-4, 5e-3),
            'rho_dust': rng.uniform(1000, 1800),
            'cp_dust': rng.uniform(600, 900),
            'Et': rng.uniform(500, 5000),
            'ssalb_therm': rng.uniform(0.0, 0.3),
            'g_therm': rng.uniform(-0.2, 0.2),
            'latitude': rng.uniform(0, np.pi / 2),
            'T_bottom': rng.uniform(200, 300),
            'ndays': 2,
            'tsteps_day': 16000,
        })
    return params


def _asteroid_params(n: int) -> list:
    """Asteroid regime: fast rotator, varied distance and TI."""
    rng = np.random.default_rng(123)
    params = []
    for _ in range(n):
        params.append({
            'P': rng.uniform(10000, 40000),    # 3-11 hours
            'R': rng.uniform(0.9, 1.5),        # AU
            'k_dust': rng.uniform(1e-4, 1e-2),
            'rho_dust': rng.uniform(800, 2500),
            'cp_dust': rng.uniform(500, 900),
            'Et': rng.uniform(100, 10000),
            'ssalb_therm': rng.uniform(0.0, 0.25),
            'g_therm': rng.uniform(-0.1, 0.1),
            'latitude': rng.uniform(0, np.pi / 3),
            'T_bottom': rng.uniform(100, 350),
            'ndays': 3,
            'tsteps_day': 12000,
        })
    return params


def _lab_params(n: int) -> list:
    """Lab sample regime: very fast rotation equivalent, controlled bottom T."""
    rng = np.random.default_rng(456)
    params = []
    for _ in range(n):
        params.append({
            'P': rng.uniform(60, 3600),        # 1 min to 1 hour
            'R': rng.uniform(0.5, 2.0),        # Varied equivalent flux
            'k_dust': rng.uniform(1e-4, 1e-2),
            'rho_dust': rng.uniform(800, 2500),
            'cp_dust': rng.uniform(500, 1000),
            'Et': rng.uniform(200, 8000),
            'ssalb_therm': rng.uniform(0.0, 0.3),
            'g_therm': rng.uniform(-0.2, 0.2),
            'latitude': 0.0,
            'T_bottom': rng.uniform(200, 500),
            'ndays': 5,
            'tsteps_day': 8000,
        })
    return params


REGIME_GENERATORS = {
    'lunar': _lunar_params,
    'asteroid': _asteroid_params,
    'lab': _lab_params,
}


# --------------------------------------------------------------------------- #
#  Core runner
# --------------------------------------------------------------------------- #

def run_single_recording(sim_id: int, param_overrides: dict,
                         output_file: str) -> dict:
    """
    Run one thermal simulation with DISORT recording enabled.

    Args:
        sim_id: Simulation identifier
        param_overrides: Dict of SimulationConfig fields to override
        output_file: HDF5 file to append DISORT records to

    Returns:
        dict with run metadata
    """
    t0 = time.time()

    # Build config with recording enabled
    cfg = SimulationConfig()
    cfg.use_RTE = True
    cfg.RTE_solver = 'disort'
    cfg.thermal_evolution_mode = 'two_wave'
    cfg.diurnal = True
    cfg.sun = True
    cfg.single_layer = True
    cfg.record_disort_data = True
    cfg.disort_record_file = output_file
    cfg.compute_observer_radiance = False
    cfg.last_day = True

    # Apply parameter overrides
    for key, val in param_overrides.items():
        if hasattr(cfg, key):
            setattr(cfg, key, val)
        else:
            logger.warning(f"Sim {sim_id}: unknown config key '{key}', skipping")

    try:
        sim = Simulator(cfg)
        sim.run()

        # Flush any remaining buffered records
        if hasattr(sim, 'rte_solver') and hasattr(sim.rte_solver, 'flush_recording'):
            sim.rte_solver.flush_recording()

        dt = time.time() - t0
        logger.info(f"Sim {sim_id} complete ({dt:.1f}s)")
        return {'sim_id': sim_id, 'status': 'success', 'time_s': dt,
                'params': param_overrides}

    except Exception as e:
        dt = time.time() - t0
        logger.error(f"Sim {sim_id} failed ({dt:.1f}s): {e}")
        return {'sim_id': sim_id, 'status': 'failed', 'time_s': dt,
                'error': str(e), 'params': param_overrides}


def generate_training_data(regime: str, n_sims: int, output_dir: str):
    """
    Generate training data for a given regime.

    Args:
        regime: 'lunar', 'asteroid', 'lab', or 'all'
        n_sims: Number of simulations per regime
        output_dir: Directory to write HDF5 files
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if regime == 'all':
        regimes = list(REGIME_GENERATORS.keys())
    else:
        regimes = [regime]

    total_t0 = time.time()
    all_results = []

    for reg in regimes:
        logger.info(f"=== Generating {reg} regime: {n_sims} simulations ===")
        generator = REGIME_GENERATORS[reg]
        param_list = generator(n_sims)
        h5_file = str(output_path / f'disort_records_{reg}.h5')

        for i, params in enumerate(param_list):
            result = run_single_recording(i, params, h5_file)
            all_results.append(result)

        # Report dataset size
        if Path(h5_file).exists():
            with h5py.File(h5_file, 'r') as f:
                n_records = f['inputs/mu'].shape[0] if 'inputs/mu' in f else 0
            logger.info(f"{reg}: {n_records} total DISORT records in {h5_file}")

    # Merge regime files into a single combined file
    if len(regimes) > 1:
        combined_file = str(output_path / 'disort_records_combined.h5')
        _merge_h5_files(
            [str(output_path / f'disort_records_{r}.h5') for r in regimes],
            combined_file
        )
        logger.info(f"Combined dataset: {combined_file}")

    total_dt = time.time() - total_t0
    n_success = sum(1 for r in all_results if r['status'] == 'success')
    n_failed = sum(1 for r in all_results if r['status'] == 'failed')
    logger.info(f"Done: {n_success} succeeded, {n_failed} failed, {total_dt:.0f}s total")


def _merge_h5_files(input_files: list, output_file: str):
    """Merge multiple DISORT record HDF5 files into one."""
    datasets_to_merge = [
        'inputs/T_profile', 'inputs/mu', 'inputs/F',
        'inputs/ssalb', 'inputs/g', 'inputs/tau_total',
        'outputs/source_term', 'outputs/flux_up',
    ]

    # Collect all data
    all_data = {name: [] for name in datasets_to_merge}
    for fpath in input_files:
        if not Path(fpath).exists():
            continue
        with h5py.File(fpath, 'r') as f:
            for name in datasets_to_merge:
                if name in f:
                    all_data[name].append(f[name][:])

    # Write combined
    with h5py.File(output_file, 'w') as f:
        for name, arrays in all_data.items():
            if arrays:
                combined = np.concatenate(arrays, axis=0)
                f.create_dataset(name, data=combined, compression='gzip')

    # Report
    with h5py.File(output_file, 'r') as f:
        n = f['inputs/mu'].shape[0] if 'inputs/mu' in f else 0
    logger.info(f"Merged {len(input_files)} files -> {output_file} ({n} records)")


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(
        description='Generate DISORT training data for NN surrogate')
    parser.add_argument('--regime', type=str, default='all',
                        choices=['lunar', 'asteroid', 'lab', 'all'],
                        help='Parameter regime (default: all)')
    parser.add_argument('--n_sims', type=int, default=100,
                        help='Number of simulations per regime (default: 100)')
    parser.add_argument('--output_dir', type=str, default='training_data',
                        help='Output directory for HDF5 files (default: training_data)')
    parser.add_argument('--log_level', type=str, default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'])
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s %(name)s %(levelname)s: %(message)s',
    )

    generate_training_data(args.regime, args.n_sims, args.output_dir)


if __name__ == '__main__':
    main()
