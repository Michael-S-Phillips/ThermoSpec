"""
End-to-end training script for the DISORT surrogate model.

Loads HDF5 training data, sets up the model and physics-informed loss,
trains with optional W&B logging, and saves the best checkpoint.

Usage:
    python -m nn.train_surrogate --data training_data/disort_records_combined.h5
    python -m nn.train_surrogate --data data.h5 --arch deeponet --wandb --epochs 500
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).parent.parent))

from nn.architectures import SourceTermMLP, DeepONetSurrogate
from nn.dataset import create_data_loaders
from nn.physics_constraints import PhysicsInformedLoss
from nn.training import TrainerConfig, SurrogateTrainer

logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description='Train DISORT surrogate model')

    # Data
    parser.add_argument('--data', type=str, required=True,
                        help='Path to HDF5 training data file')
    parser.add_argument('--batch_size', type=int, default=256)
    parser.add_argument('--num_workers', type=int, default=4)

    # Architecture
    parser.add_argument('--arch', type=str, default='mlp',
                        choices=['mlp', 'deeponet'],
                        help='Model architecture (default: mlp)')
    parser.add_argument('--n_depth', type=int, default=None,
                        help='Number of depth points (auto-detected from data)')
    parser.add_argument('--latent_dim', type=int, default=32,
                        help='Latent dimension for DeepONet')

    # Training
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--weight_decay', type=float, default=1e-5)
    parser.add_argument('--epochs', type=int, default=200)
    parser.add_argument('--patience', type=int, default=20)
    parser.add_argument('--checkpoint_dir', type=str, default='checkpoints')

    # Physics loss weights
    parser.add_argument('--lambda_energy', type=float, default=0.1)
    parser.add_argument('--lambda_positivity', type=float, default=1.0)
    parser.add_argument('--lambda_decay', type=float, default=0.01)
    parser.add_argument('--lambda_sb', type=float, default=0.01)

    # W&B
    parser.add_argument('--wandb', action='store_true',
                        help='Enable Weights & Biases logging')
    parser.add_argument('--wandb_project', type=str,
                        default='thermospec-surrogate')
    parser.add_argument('--wandb_run_name', type=str, default=None)

    # Device
    parser.add_argument('--device', type=str, default='cpu',
                        choices=['cpu', 'cuda', 'mps'])

    parser.add_argument('--log_level', type=str, default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'])

    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s %(name)s %(levelname)s: %(message)s',
    )

    # --- Load data ---
    logger.info(f"Loading data from {args.data}")
    train_loader, val_loader, test_loader, norm_params = create_data_loaders(
        args.data,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
    )

    # Detect n_depth from data
    sample = next(iter(train_loader))
    n_depth = sample['T_profile'].shape[1]
    if args.n_depth is not None and args.n_depth != n_depth:
        logger.warning(f"--n_depth={args.n_depth} overridden by data n_depth={n_depth}")
    args.n_depth = n_depth
    logger.info(f"Detected n_depth={n_depth} from training data")

    # --- Build model ---
    if args.arch == 'mlp':
        model = SourceTermMLP(n_depth=n_depth, n_params=5)
        logger.info(f"SourceTermMLP: {model.count_parameters()} parameters")
    elif args.arch == 'deeponet':
        model = DeepONetSurrogate(n_depth=n_depth, latent_dim=args.latent_dim)
        logger.info(f"DeepONetSurrogate: {model.count_parameters()} parameters")

    # --- Build physics-informed loss ---
    # Placeholder dz and tau_cumulative — will be refined with real grid params
    dz = np.ones(n_depth) / n_depth
    tau_cumulative = np.cumsum(dz)

    loss_fn = PhysicsInformedLoss(
        dz=dz,
        tau_cumulative=tau_cumulative,
        lambda_energy=args.lambda_energy,
        lambda_positivity=args.lambda_positivity,
        lambda_decay=args.lambda_decay,
        lambda_sb=args.lambda_sb,
    )

    # --- Trainer config ---
    trainer_cfg = TrainerConfig(
        lr=args.lr,
        weight_decay=args.weight_decay,
        epochs=args.epochs,
        patience=args.patience,
        checkpoint_dir=args.checkpoint_dir,
        use_wandb=args.wandb,
        wandb_project=args.wandb_project,
        wandb_run_name=args.wandb_run_name,
    )

    # --- Train ---
    trainer = SurrogateTrainer(
        model=model,
        loss_fn=loss_fn,
        train_loader=train_loader,
        val_loader=val_loader,
        config=trainer_cfg,
        device=args.device,
    )

    history = trainer.train()

    logger.info(f"Training complete. Best val loss: {history['best_val_loss']:.6f} "
                f"at epoch {history['best_epoch']}")
    logger.info(f"Best model saved to: {history['best_model_path']}")

    # --- Save normalization params alongside model ---
    norm_path = Path(args.checkpoint_dir) / 'norm_params.npz'
    norm_params.save(str(norm_path))
    logger.info(f"Normalization params saved to: {norm_path}")

    # --- Quick eval on test set ---
    from nn.evaluate import evaluate_per_call_accuracy
    logger.info("Evaluating on test set...")
    model.load_state_dict(
        torch.load(history['best_model_path'], map_location=args.device)['model_state_dict']
    )
    metrics = evaluate_per_call_accuracy(model, test_loader, device=args.device)
    logger.info(f"Test T_surf error: {metrics['T_surf_abs_error_mean_K']:.4f} K (mean), "
                f"{metrics['T_surf_abs_error_p95_K']:.4f} K (p95)")


if __name__ == '__main__':
    main()
