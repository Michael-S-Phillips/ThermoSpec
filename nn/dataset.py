"""
Dataset handling for DISORT surrogate training.

Loads DISORT input/output pairs from HDF5, applies normalization,
and provides train/val/test splits via PyTorch DataLoader.

HDF5 schema (written by DisortRTESolver in record_mode):
    /inputs/T_profile      [N, n_depth]   - Temperature profiles (K)
    /inputs/mu             [N]            - Solar cosine angles
    /inputs/F              [N]            - Solar illumination flags
    /inputs/ssalb          [N]            - Single-scattering albedo
    /inputs/g              [N]            - Asymmetry parameter
    /inputs/tau_total      [N]            - Total optical depth
    /outputs/source_term   [N, n_depth]   - Volumetric heating rates
    /outputs/flux_up       [N]            - Upward thermal flux at surface
"""

import logging
import numpy as np
import h5py
import torch
from torch.utils.data import Dataset, DataLoader, random_split
from pathlib import Path
from typing import Tuple, Dict, Optional

logger = logging.getLogger(__name__)


class NormalizationParams:
    """Stores normalization parameters for inputs and outputs."""

    def __init__(self):
        self.T_mean = 0.0
        self.T_std = 1.0
        self.params_mean = None  # [n_params]
        self.params_std = None
        self.source_mean = 0.0
        self.source_std = 1.0
        self.flux_mean = 0.0
        self.flux_std = 1.0

    def save(self, path: str):
        """Save normalization parameters to npz file."""
        np.savez(
            path,
            T_mean=self.T_mean, T_std=self.T_std,
            params_mean=self.params_mean, params_std=self.params_std,
            source_mean=self.source_mean, source_std=self.source_std,
            flux_mean=self.flux_mean, flux_std=self.flux_std,
        )

    @classmethod
    def load(cls, path: str):
        """Load normalization parameters from npz file."""
        data = np.load(path)
        params = cls()
        params.T_mean = float(data['T_mean'])
        params.T_std = float(data['T_std'])
        params.params_mean = data['params_mean']
        params.params_std = data['params_std']
        params.source_mean = float(data['source_mean'])
        params.source_std = float(data['source_std'])
        params.flux_mean = float(data['flux_mean'])
        params.flux_std = float(data['flux_std'])
        return params


class DisortDataset(Dataset):
    """
    PyTorch Dataset for DISORT surrogate training.

    Loads data from HDF5 file and applies optional normalization.
    """

    def __init__(self, h5_path: str, normalize: bool = True,
                 norm_params: Optional[NormalizationParams] = None):
        """
        Args:
            h5_path: Path to HDF5 data file
            normalize: Whether to normalize inputs/outputs
            norm_params: Pre-computed normalization params (if None, computed from data)
        """
        self.h5_path = Path(h5_path)
        self.normalize = normalize

        # Load data into memory
        logger.info(f"Loading dataset from {h5_path}")
        with h5py.File(h5_path, 'r') as f:
            self.T_profiles = f['inputs/T_profile'][:]
            self.mu = f['inputs/mu'][:]
            self.F = f['inputs/F'][:]
            self.ssalb = f['inputs/ssalb'][:]
            self.g = f['inputs/g'][:]
            self.tau_total = f['inputs/tau_total'][:]
            self.source_terms = f['outputs/source_term'][:]
            self.flux_up = f['outputs/flux_up'][:]

        self.n_samples = len(self.mu)
        self.n_depth = self.T_profiles.shape[1]
        logger.info(f"Loaded {self.n_samples} samples, {self.n_depth} depth points")

        # Stack scalar params: [N, 5]
        self.params = np.stack([
            self.mu, self.F, self.ssalb, self.g, self.tau_total
        ], axis=1)

        # Compute or use provided normalization
        if normalize:
            if norm_params is None:
                self.norm = self._compute_normalization()
            else:
                self.norm = norm_params
            self._apply_normalization()
        else:
            self.norm = NormalizationParams()

    def _compute_normalization(self) -> NormalizationParams:
        """Compute normalization statistics from data."""
        norm = NormalizationParams()
        norm.T_mean = float(np.mean(self.T_profiles))
        norm.T_std = float(np.std(self.T_profiles)) + 1e-8
        norm.params_mean = np.mean(self.params, axis=0)
        norm.params_std = np.std(self.params, axis=0) + 1e-8
        norm.source_mean = float(np.mean(self.source_terms))
        norm.source_std = float(np.std(self.source_terms)) + 1e-8
        norm.flux_mean = float(np.mean(self.flux_up))
        norm.flux_std = float(np.std(self.flux_up)) + 1e-8
        return norm

    def _apply_normalization(self):
        """Apply normalization to stored data."""
        self.T_profiles = (self.T_profiles - self.norm.T_mean) / self.norm.T_std
        self.params = (self.params - self.norm.params_mean) / self.norm.params_std
        self.source_terms = (self.source_terms - self.norm.source_mean) / self.norm.source_std
        self.flux_up = (self.flux_up - self.norm.flux_mean) / self.norm.flux_std

    def __len__(self):
        return self.n_samples

    def __getitem__(self, idx):
        return {
            'T_profile': torch.tensor(self.T_profiles[idx], dtype=torch.float32),
            'params': torch.tensor(self.params[idx], dtype=torch.float32),
            'source_term': torch.tensor(self.source_terms[idx], dtype=torch.float32),
            'flux_up': torch.tensor(self.flux_up[idx], dtype=torch.float32),
        }


def create_data_loaders(
    h5_path: str,
    batch_size: int = 256,
    train_frac: float = 0.8,
    val_frac: float = 0.1,
    num_workers: int = 4,
    seed: int = 42,
) -> Tuple[DataLoader, DataLoader, DataLoader, NormalizationParams]:
    """
    Create train/val/test DataLoaders from HDF5 file.

    Args:
        h5_path: Path to HDF5 data file
        batch_size: Batch size for training
        train_frac: Fraction of data for training
        val_frac: Fraction of data for validation (rest is test)
        num_workers: DataLoader workers
        seed: Random seed for reproducible splits

    Returns:
        (train_loader, val_loader, test_loader, norm_params)
    """
    # Load full dataset (normalized)
    full_dataset = DisortDataset(h5_path, normalize=True)
    norm_params = full_dataset.norm

    # Split into train/val/test
    n = len(full_dataset)
    n_train = int(n * train_frac)
    n_val = int(n * val_frac)
    n_test = n - n_train - n_val

    generator = torch.Generator().manual_seed(seed)
    train_ds, val_ds, test_ds = random_split(
        full_dataset, [n_train, n_val, n_test], generator=generator
    )

    logger.info(f"Dataset split: train={n_train}, val={n_val}, test={n_test}")

    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True,
                              num_workers=num_workers, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False,
                            num_workers=num_workers, pin_memory=True)
    test_loader = DataLoader(test_ds, batch_size=batch_size, shuffle=False,
                             num_workers=num_workers, pin_memory=True)

    return train_loader, val_loader, test_loader, norm_params
