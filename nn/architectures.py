"""
Neural network architectures for DISORT surrogate modeling.

Provides a Modified DeepONet architecture that maps:
    (T_profile, depth, mu, F, ssalb, g, tau_total) -> (source_term, flux_up)

The architecture consists of:
    - Branch net: Encodes temperature profile T(z) -> latent basis coefficients
    - Trunk net: Encodes (depth z, mu, F, ssalb, g, tau_total) -> latent basis functions
    - Flux head: Predicts scalar flux_up from branch encoding
"""

import torch
import torch.nn as nn
import numpy as np


class MLP(nn.Module):
    """Simple multi-layer perceptron with configurable depth and width."""

    def __init__(self, input_dim, output_dim, hidden_dims, activation=nn.GELU):
        super().__init__()
        layers = []
        prev_dim = input_dim
        for h in hidden_dims:
            layers.append(nn.Linear(prev_dim, h))
            layers.append(activation())
            prev_dim = h
        layers.append(nn.Linear(prev_dim, output_dim))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


class DeepONetSurrogate(nn.Module):
    """
    Modified DeepONet for DISORT surrogate.

    Branch net encodes the temperature profile (function input).
    Trunk net encodes the query coordinates and physical parameters.
    Output is the dot product of branch and trunk latent representations,
    plus a separate flux head.

    Args:
        n_depth: Number of depth points in temperature profile
        trunk_input_dim: Dimension of trunk inputs (depth, mu, F, ssalb, g, tau)
        latent_dim: Dimension of latent basis (p in DeepONet)
        branch_hidden: List of hidden layer sizes for branch net
        trunk_hidden: List of hidden layer sizes for trunk net
        flux_hidden: List of hidden layer sizes for flux prediction head
    """

    def __init__(
        self,
        n_depth: int = 50,
        trunk_input_dim: int = 6,
        latent_dim: int = 32,
        branch_hidden: list = None,
        trunk_hidden: list = None,
        flux_hidden: list = None,
    ):
        super().__init__()

        if branch_hidden is None:
            branch_hidden = [128, 128]
        if trunk_hidden is None:
            trunk_hidden = [64, 64]
        if flux_hidden is None:
            flux_hidden = [64, 32]

        # Branch net: T(z) -> latent basis coefficients [batch, latent_dim]
        self.branch_net = MLP(n_depth, latent_dim, branch_hidden)

        # Trunk net: (z, mu, F, ssalb, g, tau) -> latent basis functions [batch, n_query, latent_dim]
        self.trunk_net = MLP(trunk_input_dim, latent_dim, trunk_hidden)

        # Flux head: branch encoding -> scalar flux_up
        self.flux_head = MLP(latent_dim, 1, flux_hidden)

        # Bias for source term output
        self.source_bias = nn.Parameter(torch.zeros(1))

        self._n_depth = n_depth
        self._latent_dim = latent_dim

    def forward(self, T_profile, trunk_inputs):
        """
        Forward pass.

        Args:
            T_profile: Temperature profile [batch, n_depth]
            trunk_inputs: Query parameters [batch, n_query, trunk_input_dim]
                          Each query point contains (z, mu, F, ssalb, g, tau_total)

        Returns:
            source_term: Volumetric heating rate [batch, n_query]
            flux_up: Upward thermal flux [batch, 1]
        """
        # Branch encoding
        branch_out = self.branch_net(T_profile)  # [batch, latent_dim]

        # Trunk encoding
        trunk_out = self.trunk_net(trunk_inputs)  # [batch, n_query, latent_dim]

        # Source term: dot product of branch and trunk
        # branch_out: [batch, latent_dim] -> [batch, 1, latent_dim]
        source_term = torch.sum(
            branch_out.unsqueeze(1) * trunk_out, dim=-1
        ) + self.source_bias  # [batch, n_query]

        # Flux prediction from branch encoding
        flux_up = self.flux_head(branch_out)  # [batch, 1]

        return source_term, flux_up

    def count_parameters(self):
        """Return total number of trainable parameters."""
        return sum(p.numel() for p in self.parameters() if p.requires_grad)


class SourceTermMLP(nn.Module):
    """
    Simpler MLP-based surrogate that directly maps
    (T_profile, mu, F, ssalb, g, tau) -> (source_term_profile, flux_up).

    Useful as a baseline comparison for the DeepONet.

    Args:
        n_depth: Number of depth points
        n_params: Number of scalar physical parameters (mu, F, ssalb, g, tau)
        hidden_dims: List of hidden layer sizes
    """

    def __init__(self, n_depth=50, n_params=5, hidden_dims=None):
        super().__init__()

        if hidden_dims is None:
            hidden_dims = [256, 256, 128]

        input_dim = n_depth + n_params
        output_dim = n_depth + 1  # source_term at each depth + flux_up

        self.net = MLP(input_dim, output_dim, hidden_dims)
        self._n_depth = n_depth

    def forward(self, T_profile, params):
        """
        Args:
            T_profile: [batch, n_depth]
            params: [batch, n_params] containing (mu, F, ssalb, g, tau_total)

        Returns:
            source_term: [batch, n_depth]
            flux_up: [batch, 1]
        """
        x = torch.cat([T_profile, params], dim=-1)
        out = self.net(x)
        source_term = out[:, :self._n_depth]
        flux_up = out[:, self._n_depth:]
        return source_term, flux_up

    def count_parameters(self):
        return sum(p.numel() for p in self.parameters() if p.requires_grad)
