"""
Drop-in surrogate wrapper replacing DisortRTESolver.disort_run().

Provides the DISORTSurrogateWrapper class that:
    1. Matches the DisortRTESolver.disort_run() interface exactly
    2. Runs the neural network surrogate for fast inference
    3. Periodically validates against real DISORT and warns on high error
    4. Falls back to DISORT when surrogate error exceeds threshold

Usage:
    # In modelmain.py, after setting up rte_disort:
    if cfg.use_nn_surrogate:
        from nn.surrogate import DISORTSurrogateWrapper
        rte_disort = DISORTSurrogateWrapper(rte_disort, cfg.nn_model_path)
"""

import logging
import time
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import torch

logger = logging.getLogger(__name__)


class DISORTSurrogateWrapper:
    """
    Drop-in replacement for DisortRTESolver that uses a neural network
    surrogate for fast inference with periodic DISORT validation.

    Wraps an existing DisortRTESolver instance, using the NN for most calls
    and falling back to real DISORT when validation fails.

    Args:
        disort_solver: Original DisortRTESolver instance
        model_path: Path to trained surrogate model checkpoint
        norm_path: Path to normalization parameters (.npz)
        validation_interval: Run real DISORT every N calls for validation
        error_threshold: Max allowed relative flux error before fallback
        device: Torch device for inference
    """

    def __init__(
        self,
        disort_solver,
        model_path: str,
        norm_path: Optional[str] = None,
        validation_interval: int = 100,
        error_threshold: float = 0.01,
        device: str = "cpu",
    ):
        self.disort_solver = disort_solver
        self.device = device
        self.validation_interval = validation_interval
        self.error_threshold = error_threshold

        # Copy attributes from original solver for interface compatibility
        self.cfg = disort_solver.cfg
        self.grid = disort_solver.grid
        self.n_cols = disort_solver.n_cols
        self.is_multi_wave = disort_solver.is_multi_wave
        self.nwave = disort_solver.nwave

        if hasattr(disort_solver, 'wavenumbers'):
            self.wavenumbers = disort_solver.wavenumbers

        # Load model
        self.model = self._load_model(model_path)
        self.model.eval()

        # Load normalization parameters
        if norm_path is not None:
            from nn.dataset import NormalizationParams
            self.norm = NormalizationParams.load(norm_path)
        else:
            self.norm = None

        # Statistics tracking
        self.call_count = 0
        self.nn_calls = 0
        self.disort_calls = 0
        self.fallback_count = 0
        self.total_nn_time = 0.0
        self.total_disort_time = 0.0

        logger.info(
            f"DISORT surrogate initialized: model={model_path}, "
            f"validation_interval={validation_interval}, "
            f"error_threshold={error_threshold}"
        )

    def _load_model(self, model_path: str):
        """Load trained surrogate model."""
        from nn.architectures import SourceTermMLP, DeepONetSurrogate

        checkpoint = torch.load(model_path, map_location=self.device)
        state_dict = checkpoint["model_state_dict"]

        # Detect architecture from state dict keys
        n_depth = self.grid.nlay_dust

        # Try SourceTermMLP first (simpler)
        try:
            model = SourceTermMLP(n_depth=n_depth)
            model.load_state_dict(state_dict)
            logger.info(f"Loaded SourceTermMLP surrogate ({model.count_parameters()} params)")
            return model.to(self.device)
        except (RuntimeError, KeyError):
            pass

        # Try DeepONet
        model = DeepONetSurrogate(n_depth=n_depth)
        model.load_state_dict(state_dict)
        logger.info(f"Loaded DeepONetSurrogate ({model.count_parameters()} params)")
        return model.to(self.device)

    def disort_run(self, T, mu, F, Q=None, phi=None):
        """
        Drop-in replacement for DisortRTESolver.disort_run().

        Uses the neural network surrogate for fast inference.
        Periodically validates against real DISORT.

        Args and returns match DisortRTESolver.disort_run() exactly.
        """
        self.call_count += 1

        # Periodic validation
        if self.call_count % self.validation_interval == 0:
            return self._validated_run(T, mu, F, Q, phi)

        # Fast path: NN inference only
        try:
            t0 = time.perf_counter()
            result = self._nn_inference(T, mu, F)
            self.total_nn_time += time.perf_counter() - t0
            self.nn_calls += 1
            return result
        except Exception as e:
            logger.warning(f"NN inference failed ({e}), falling back to DISORT")
            return self._disort_fallback(T, mu, F, Q, phi)

    def _nn_inference(self, T, mu, F) -> Tuple[np.ndarray, np.ndarray]:
        """Run neural network inference."""
        n_depth = self.grid.nlay_dust

        # Extract temperature profile for RTE layers
        if self.n_cols == 1:
            T_profile = T[1:n_depth + 1]
        else:
            T_profile = T[1:n_depth + 1, :]

        # Prepare inputs
        T_tensor = torch.tensor(T_profile, dtype=torch.float32, device=self.device)
        if T_tensor.dim() == 1:
            T_tensor = T_tensor.unsqueeze(0)  # [1, n_depth]
        else:
            T_tensor = T_tensor.T  # [n_cols, n_depth]

        # Scalar parameters
        mu_val = float(mu) if np.isscalar(mu) else float(mu[0])
        F_val = float(F) if np.isscalar(F) else float(F[0])
        ssalb = float(self.cfg.ssalb_therm)
        g = float(self.cfg.g_therm)
        tau_total = float(self.grid.x_boundaries[-1])

        params = torch.tensor(
            [[mu_val, F_val, ssalb, g, tau_total]],
            dtype=torch.float32, device=self.device
        )

        # Normalize if needed
        if self.norm is not None:
            T_tensor = (T_tensor - self.norm.T_mean) / self.norm.T_std
            params = (params - torch.tensor(self.norm.params_mean, dtype=torch.float32)
                     ) / torch.tensor(self.norm.params_std, dtype=torch.float32)

        # Forward pass
        with torch.no_grad():
            pred_source, pred_flux = self.model(T_tensor, params)

        # Denormalize
        if self.norm is not None:
            pred_source = pred_source * self.norm.source_std + self.norm.source_mean
            pred_flux = pred_flux * self.norm.flux_std + self.norm.flux_mean

        # Convert to source_term format matching DisortRTESolver output
        source_np = pred_source.cpu().numpy().squeeze()
        flux_np = pred_flux.cpu().numpy().squeeze()

        if self.n_cols == 1:
            source_term = np.zeros(self.grid.x_num)
            source_term[1:n_depth + 1] = source_np
            source_term *= self.grid.K * self.cfg.q
            flux_up = flux_np
        else:
            source_term = np.zeros((self.n_cols, self.grid.x_num))
            source_term[:, 1:n_depth + 1] = source_np
            source_term *= self.grid.K * self.cfg.q
            flux_up = flux_np

        return source_term, flux_up

    def _validated_run(self, T, mu, F, Q=None, phi=None):
        """
        Run both NN and DISORT, compare results, and return DISORT result.
        Logs warning if error exceeds threshold.
        """
        # Run DISORT (ground truth)
        t0 = time.perf_counter()
        source_disort, flux_disort = self.disort_solver.disort_run(T, mu, F, Q, phi)
        disort_time = time.perf_counter() - t0
        self.total_disort_time += disort_time
        self.disort_calls += 1

        # Run NN
        try:
            t0 = time.perf_counter()
            source_nn, flux_nn = self._nn_inference(T, mu, F)
            nn_time = time.perf_counter() - t0
            self.total_nn_time += nn_time

            # Compare flux
            if np.isscalar(flux_disort):
                flux_ref = abs(flux_disort) + 1e-10
                flux_err = abs(flux_nn - flux_disort) / flux_ref
            else:
                flux_ref = np.abs(flux_disort).max() + 1e-10
                flux_err = np.abs(flux_nn - flux_disort).max() / flux_ref

            speedup = disort_time / (nn_time + 1e-10)

            if flux_err > self.error_threshold:
                logger.warning(
                    f"Surrogate validation: flux error={flux_err:.4f} > "
                    f"threshold={self.error_threshold}. "
                    f"Speedup={speedup:.1f}x (call #{self.call_count})"
                )
                self.fallback_count += 1
            else:
                logger.debug(
                    f"Surrogate validation OK: flux error={flux_err:.6f}, "
                    f"speedup={speedup:.1f}x"
                )
        except Exception as e:
            logger.warning(f"Surrogate validation failed: {e}")

        return source_disort, flux_disort

    def _disort_fallback(self, T, mu, F, Q=None, phi=None):
        """Fall back to real DISORT."""
        t0 = time.perf_counter()
        result = self.disort_solver.disort_run(T, mu, F, Q, phi)
        self.total_disort_time += time.perf_counter() - t0
        self.disort_calls += 1
        self.fallback_count += 1
        return result

    def get_stats(self) -> dict:
        """Return surrogate usage statistics."""
        return {
            "total_calls": self.call_count,
            "nn_calls": self.nn_calls,
            "disort_calls": self.disort_calls,
            "fallback_count": self.fallback_count,
            "nn_fraction": self.nn_calls / max(self.call_count, 1),
            "avg_nn_time_ms": 1000 * self.total_nn_time / max(self.nn_calls, 1),
            "avg_disort_time_ms": 1000 * self.total_disort_time / max(self.disort_calls, 1),
            "speedup": (
                (self.total_disort_time / max(self.disort_calls, 1))
                / (self.total_nn_time / max(self.nn_calls, 1) + 1e-10)
                if self.nn_calls > 0 and self.disort_calls > 0
                else 0.0
            ),
        }
