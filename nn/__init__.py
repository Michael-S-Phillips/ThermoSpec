"""
ThermoSpec Neural Network Surrogate Module.

Provides a physics-informed neural network (PINN) surrogate for the DISORT
radiative transfer solver, targeting 10-50x speedup while maintaining
sub-Kelvin accuracy.

Submodules:
    architectures  - Model definitions (DeepONet, MLP variants)
    dataset        - HDF5 dataset loading, normalization, train/val/test splits
    training       - Training loop with W&B logging
    surrogate      - Drop-in wrapper replacing DisortRTESolver
    physics_constraints - Physics-informed loss terms
    evaluate       - Accuracy metrics, comparison tools
"""
