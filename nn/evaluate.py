"""
Evaluation tools for DISORT surrogate accuracy.

Provides metrics for:
    - Per-call accuracy (source term RMSE, flux error, T_surf error)
    - Full-simulation comparison (diurnal T_surf curves, T(z) profiles)
    - Speed benchmarks (DISORT vs NN wall-clock)
    - Error distribution analysis (percentiles)
"""

import logging
import time
from typing import Dict, List, Optional

import numpy as np
import torch
from torch.utils.data import DataLoader

logger = logging.getLogger(__name__)

SIGMA_SB = 5.670373e-8


def evaluate_per_call_accuracy(
    model, test_loader: DataLoader, norm_params=None, device: str = "cpu"
) -> Dict[str, float]:
    """
    Evaluate surrogate accuracy on test dataset.

    Args:
        model: Trained surrogate model
        test_loader: Test DataLoader
        norm_params: Normalization parameters for denormalization
        device: Torch device

    Returns:
        Dict with metrics: source_rmse, source_rel_rmse, flux_rel_error,
        flux_abs_error_mean, T_surf_abs_error_K, percentiles
    """
    model.eval()
    all_source_errors = []
    all_flux_errors = []
    all_flux_rel_errors = []
    all_Tsurf_errors = []

    with torch.no_grad():
        for batch in test_loader:
            T_profile = batch["T_profile"].to(device)
            params = batch["params"].to(device)
            target_source = batch["source_term"].numpy()
            target_flux = batch["flux_up"].numpy()

            pred_source, pred_flux = model(T_profile, params)
            pred_source = pred_source.cpu().numpy()
            pred_flux = pred_flux.cpu().numpy().squeeze()

            # Denormalize if needed
            if norm_params is not None:
                pred_source = pred_source * norm_params.source_std + norm_params.source_mean
                pred_flux = pred_flux * norm_params.flux_std + norm_params.flux_mean
                target_source = target_source * norm_params.source_std + norm_params.source_mean
                target_flux = target_flux * norm_params.flux_std + norm_params.flux_mean

            # Source term errors
            source_err = np.abs(pred_source - target_source)
            all_source_errors.append(source_err)

            # Flux errors
            flux_err = np.abs(pred_flux - target_flux)
            all_flux_errors.append(flux_err)

            flux_ref = np.abs(target_flux) + 1e-10
            flux_rel_err = flux_err / flux_ref
            all_flux_rel_errors.append(flux_rel_err)

            # T_surf errors (from flux via Stefan-Boltzmann)
            T_surf_pred = (np.maximum(pred_flux, 0) / SIGMA_SB) ** 0.25
            T_surf_target = (np.maximum(target_flux, 0) / SIGMA_SB) ** 0.25
            all_Tsurf_errors.append(np.abs(T_surf_pred - T_surf_target))

    # Aggregate
    source_errors = np.concatenate(all_source_errors, axis=0)
    flux_errors = np.concatenate(all_flux_errors)
    flux_rel_errors = np.concatenate(all_flux_rel_errors)
    Tsurf_errors = np.concatenate(all_Tsurf_errors)

    metrics = {
        "source_rmse": float(np.sqrt(np.mean(source_errors ** 2))),
        "source_rel_rmse": float(
            np.sqrt(np.mean(source_errors ** 2))
            / (np.std(np.concatenate([b["source_term"].numpy() for b in test_loader.dataset])) + 1e-10)
        ) if hasattr(test_loader.dataset, '__getitem__') else 0.0,
        "flux_abs_error_mean": float(np.mean(flux_errors)),
        "flux_rel_error_mean": float(np.mean(flux_rel_errors)),
        "flux_rel_error_p50": float(np.percentile(flux_rel_errors, 50)),
        "flux_rel_error_p95": float(np.percentile(flux_rel_errors, 95)),
        "flux_rel_error_p99": float(np.percentile(flux_rel_errors, 99)),
        "T_surf_abs_error_mean_K": float(np.mean(Tsurf_errors)),
        "T_surf_abs_error_max_K": float(np.max(Tsurf_errors)),
        "T_surf_abs_error_p95_K": float(np.percentile(Tsurf_errors, 95)),
        "T_surf_abs_error_p99_K": float(np.percentile(Tsurf_errors, 99)),
    }

    logger.info("Per-call accuracy metrics:")
    for k, v in metrics.items():
        logger.info(f"  {k}: {v:.6f}")

    return metrics


def benchmark_speed(
    model, test_loader: DataLoader, disort_solver=None,
    n_calls: int = 1000, device: str = "cpu"
) -> Dict[str, float]:
    """
    Benchmark NN surrogate inference speed vs DISORT.

    Args:
        model: Trained surrogate model
        test_loader: Test DataLoader for sample inputs
        disort_solver: Original DISORT solver for comparison (optional)
        n_calls: Number of inference calls to benchmark
        device: Torch device

    Returns:
        Dict with timing metrics
    """
    model.eval()

    # Get a sample batch
    sample = next(iter(test_loader))
    T_profile = sample["T_profile"][:1].to(device)
    params = sample["params"][:1].to(device)

    # Warmup
    with torch.no_grad():
        for _ in range(10):
            model(T_profile, params)

    # Benchmark NN
    t0 = time.perf_counter()
    with torch.no_grad():
        for _ in range(n_calls):
            model(T_profile, params)
    nn_total = time.perf_counter() - t0
    nn_per_call = nn_total / n_calls

    metrics = {
        "nn_total_time_s": nn_total,
        "nn_per_call_ms": nn_per_call * 1000,
        "nn_calls_per_sec": n_calls / nn_total,
    }

    logger.info(f"NN benchmark: {nn_per_call*1000:.3f} ms/call ({n_calls} calls)")

    # Benchmark DISORT if available
    if disort_solver is not None:
        n_disort = min(n_calls, 100)  # Fewer calls for DISORT (much slower)

        # Create synthetic inputs matching DisortRTESolver.disort_run() interface
        T_dummy = np.full(disort_solver.grid.x_num, 300.0)
        mu_dummy = 0.5
        F_dummy = 1.0

        t0 = time.perf_counter()
        for _ in range(n_disort):
            disort_solver.disort_run(T_dummy, mu_dummy, F_dummy)
        disort_total = time.perf_counter() - t0
        disort_per_call = disort_total / n_disort

        speedup = disort_per_call / nn_per_call

        metrics["disort_per_call_ms"] = disort_per_call * 1000
        metrics["speedup"] = speedup

        logger.info(
            f"DISORT benchmark: {disort_per_call*1000:.3f} ms/call ({n_disort} calls)"
        )
        logger.info(f"Speedup: {speedup:.1f}x")

    return metrics


def compare_full_simulation(
    sim_disort, sim_surrogate, time_indices: Optional[List[int]] = None
) -> Dict[str, float]:
    """
    Compare full diurnal simulation results between DISORT and surrogate.

    Args:
        sim_disort: Simulator results using DISORT
        sim_surrogate: Simulator results using NN surrogate
        time_indices: Specific output time indices to compare (None = all)

    Returns:
        Dict with comparison metrics
    """
    T_surf_disort = sim_disort.T_surf_out
    T_surf_nn = sim_surrogate.T_surf_out

    if time_indices is not None:
        T_surf_disort = T_surf_disort[time_indices]
        T_surf_nn = T_surf_nn[time_indices]

    # Surface temperature comparison
    T_surf_diff = np.abs(T_surf_nn - T_surf_disort)

    # Depth profile comparison (at all output times)
    T_profile_disort = sim_disort.T_out
    T_profile_nn = sim_surrogate.T_out
    T_profile_diff = np.abs(T_profile_nn - T_profile_disort)

    metrics = {
        "T_surf_max_error_K": float(np.max(T_surf_diff)),
        "T_surf_mean_error_K": float(np.mean(T_surf_diff)),
        "T_surf_rms_error_K": float(np.sqrt(np.mean(T_surf_diff ** 2))),
        "T_profile_max_error_K": float(np.max(T_profile_diff)),
        "T_profile_mean_error_K": float(np.mean(T_profile_diff)),
    }

    logger.info("Full simulation comparison:")
    for k, v in metrics.items():
        logger.info(f"  {k}: {v:.4f}")

    return metrics
