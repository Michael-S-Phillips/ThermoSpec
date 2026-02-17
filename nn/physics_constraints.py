"""
Physics-informed loss terms for DISORT surrogate training.

These constraints encode known physical laws that the surrogate must satisfy,
improving generalization and physical consistency.

Loss terms:
    1. Energy conservation: integral(source_term * dz) ~ boundary flux balance
    2. Flux positivity: upward thermal flux >= 0
    3. Deep-layer decay: source_term -> 0 beyond ~3 optical depths (nighttime)
    4. Stefan-Boltzmann consistency: flux_up ~ sigma * T_surf^4
"""

import torch
import numpy as np


SIGMA_SB = 5.670373e-8  # Stefan-Boltzmann constant (W/m^2/K^4)


def energy_conservation_loss(source_term, flux_up, dz, flux_down_top=0.0):
    """
    Energy conservation: the integral of the source term over depth
    should balance the net boundary fluxes.

    integral(S dz) = F_up_surface - F_down_top - F_up_bottom

    For a semi-infinite medium with no bottom flux:
        integral(S dz) ~ F_up_surface - F_down_top

    Args:
        source_term: Predicted heating rates [batch, n_depth]
        flux_up: Predicted upward flux at surface [batch, 1]
        dz: Layer thicknesses [n_depth] (in physical units)
        flux_down_top: Downward flux at top boundary [batch, 1] or scalar

    Returns:
        Scalar loss (MSE of energy imbalance)
    """
    # Integrate source term over depth
    dz_tensor = torch.as_tensor(dz, dtype=source_term.dtype, device=source_term.device)
    integrated_source = torch.sum(source_term * dz_tensor.unsqueeze(0), dim=1, keepdim=True)

    # Energy balance residual
    residual = integrated_source - (flux_up - flux_down_top)

    return torch.mean(residual ** 2)


def flux_positivity_loss(flux_up):
    """
    Enforce that upward thermal flux is non-negative.
    Uses a soft penalty (ReLU on negative values).

    Args:
        flux_up: Predicted upward flux [batch, 1]

    Returns:
        Scalar loss
    """
    return torch.mean(torch.relu(-flux_up) ** 2)


def deep_layer_decay_loss(source_term, tau_cumulative, tau_threshold=3.0, night_mask=None):
    """
    Enforce that source term decays to near-zero beyond a few optical depths,
    especially at nighttime when there is no solar illumination.

    Args:
        source_term: Predicted heating rates [batch, n_depth]
        tau_cumulative: Cumulative optical depth at each layer [n_depth]
        tau_threshold: Optical depth beyond which source should be negligible
        night_mask: Boolean mask [batch] indicating nighttime samples (F=0)

    Returns:
        Scalar loss
    """
    tau_tensor = torch.as_tensor(tau_cumulative, dtype=source_term.dtype,
                                  device=source_term.device)

    # Mask for deep layers
    deep_mask = (tau_tensor > tau_threshold).unsqueeze(0)  # [1, n_depth]

    if night_mask is not None:
        # Apply only to nighttime samples
        night_mask = night_mask.unsqueeze(1)  # [batch, 1]
        deep_mask = deep_mask & night_mask

    # Penalize non-zero source terms in deep layers
    deep_source = source_term * deep_mask.float()
    return torch.mean(deep_source ** 2)


def stefan_boltzmann_consistency_loss(flux_up, T_surface, emissivity=1.0, weight=0.1):
    """
    Soft constraint: predicted flux_up should be approximately consistent
    with Stefan-Boltzmann law for the surface temperature.

    This is NOT an exact constraint (RTE flux != sigma*T^4 in general),
    but provides a useful regularization signal.

    Args:
        flux_up: Predicted upward flux [batch, 1]
        T_surface: Surface temperature (first depth layer) [batch]
        emissivity: Effective emissivity (scalar or [batch])
        weight: Relative weight of this soft constraint

    Returns:
        Scalar loss
    """
    sb_flux = emissivity * SIGMA_SB * T_surface.unsqueeze(1) ** 4
    # Relative error (soft constraint)
    relative_error = (flux_up - sb_flux) / (sb_flux + 1e-10)
    return weight * torch.mean(relative_error ** 2)


class PhysicsInformedLoss(torch.nn.Module):
    """
    Combined physics-informed loss for DISORT surrogate training.

    Combines data-driven MSE loss with physics constraint terms.

    Args:
        dz: Layer thicknesses [n_depth]
        tau_cumulative: Cumulative optical depth [n_depth]
        lambda_energy: Weight for energy conservation loss
        lambda_positivity: Weight for flux positivity loss
        lambda_decay: Weight for deep-layer decay loss
        lambda_sb: Weight for Stefan-Boltzmann consistency loss
    """

    def __init__(self, dz, tau_cumulative,
                 lambda_energy=0.1,
                 lambda_positivity=1.0,
                 lambda_decay=0.01,
                 lambda_sb=0.01):
        super().__init__()
        self.register_buffer('dz', torch.tensor(dz, dtype=torch.float32))
        self.register_buffer('tau_cumulative',
                             torch.tensor(tau_cumulative, dtype=torch.float32))
        self.lambda_energy = lambda_energy
        self.lambda_positivity = lambda_positivity
        self.lambda_decay = lambda_decay
        self.lambda_sb = lambda_sb

    def forward(self, pred_source, pred_flux, target_source, target_flux,
                T_profile=None, F_values=None):
        """
        Compute combined loss.

        Args:
            pred_source: Predicted source terms [batch, n_depth]
            pred_flux: Predicted flux_up [batch, 1]
            target_source: Ground truth source terms [batch, n_depth]
            target_flux: Ground truth flux_up [batch, 1]
            T_profile: Temperature profiles [batch, n_depth] (for S-B constraint)
            F_values: Solar illumination flags [batch] (for nighttime detection)

        Returns:
            total_loss, loss_dict
        """
        # Data-driven losses
        source_mse = torch.nn.functional.mse_loss(pred_source, target_source)
        flux_mse = torch.nn.functional.mse_loss(pred_flux, target_flux)

        # Physics constraints
        energy_loss = energy_conservation_loss(pred_source, pred_flux, self.dz)
        positivity_loss = flux_positivity_loss(pred_flux)

        night_mask = (F_values < 0.01) if F_values is not None else None
        decay_loss = deep_layer_decay_loss(
            pred_source, self.tau_cumulative, night_mask=night_mask
        )

        total = (source_mse + flux_mse
                 + self.lambda_energy * energy_loss
                 + self.lambda_positivity * positivity_loss
                 + self.lambda_decay * decay_loss)

        loss_dict = {
            'source_mse': source_mse.item(),
            'flux_mse': flux_mse.item(),
            'energy_conservation': energy_loss.item(),
            'flux_positivity': positivity_loss.item(),
            'deep_decay': decay_loss.item(),
        }

        # Optional Stefan-Boltzmann consistency
        if T_profile is not None and self.lambda_sb > 0:
            sb_loss = stefan_boltzmann_consistency_loss(
                pred_flux, T_profile[:, 0]
            )
            total += self.lambda_sb * sb_loss
            loss_dict['stefan_boltzmann'] = sb_loss.item()

        loss_dict['total'] = total.item()

        return total, loss_dict
