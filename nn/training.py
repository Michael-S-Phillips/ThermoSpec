"""
Training loop for DISORT surrogate model with optional W&B logging.

Supports:
    - Physics-informed loss with configurable constraint weights
    - Learning rate scheduling (cosine annealing)
    - Early stopping based on validation loss
    - Checkpoint saving/loading
    - W&B metric tracking (if enabled)
"""

import logging
import time
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

logger = logging.getLogger(__name__)


class TrainerConfig:
    """Training hyperparameters."""

    def __init__(
        self,
        lr: float = 1e-3,
        weight_decay: float = 1e-5,
        epochs: int = 200,
        patience: int = 20,
        min_delta: float = 1e-6,
        checkpoint_dir: str = "checkpoints",
        use_wandb: bool = False,
        wandb_project: str = "thermospec-surrogate",
        wandb_run_name: Optional[str] = None,
    ):
        self.lr = lr
        self.weight_decay = weight_decay
        self.epochs = epochs
        self.patience = patience
        self.min_delta = min_delta
        self.checkpoint_dir = Path(checkpoint_dir)
        self.use_wandb = use_wandb
        self.wandb_project = wandb_project
        self.wandb_run_name = wandb_run_name


class SurrogateTrainer:
    """
    Trains a DISORT surrogate model with physics-informed constraints.
    """

    def __init__(
        self,
        model: nn.Module,
        loss_fn: nn.Module,
        train_loader: DataLoader,
        val_loader: DataLoader,
        config: TrainerConfig,
        device: str = "cpu",
    ):
        self.model = model.to(device)
        self.loss_fn = loss_fn.to(device)
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.config = config
        self.device = device

        self.optimizer = torch.optim.AdamW(
            model.parameters(), lr=config.lr, weight_decay=config.weight_decay
        )
        self.scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            self.optimizer, T_max=config.epochs, eta_min=config.lr * 0.01
        )

        # Early stopping state
        self.best_val_loss = float("inf")
        self.patience_counter = 0
        self.best_epoch = 0

        # W&B setup
        self.wandb_run = None
        if config.use_wandb:
            self._init_wandb()

        # Ensure checkpoint directory exists
        config.checkpoint_dir.mkdir(parents=True, exist_ok=True)

    def _init_wandb(self):
        """Initialize Weights & Biases logging."""
        try:
            import wandb

            self.wandb_run = wandb.init(
                project=self.config.wandb_project,
                name=self.config.wandb_run_name,
                config={
                    "lr": self.config.lr,
                    "weight_decay": self.config.weight_decay,
                    "epochs": self.config.epochs,
                    "model_params": self.model.count_parameters(),
                    "architecture": self.model.__class__.__name__,
                },
            )
            wandb.watch(self.model, log="gradients", log_freq=100)
            logger.info(f"W&B initialized: {wandb.run.url}")
        except ImportError:
            logger.warning("wandb not installed. Logging disabled.")
            self.config.use_wandb = False

    def train(self) -> Dict:
        """
        Run full training loop.

        Returns:
            Dict with training history and best model path.
        """
        logger.info(
            f"Starting training: {self.config.epochs} epochs, "
            f"{self.model.count_parameters()} parameters"
        )

        history = {"train_loss": [], "val_loss": [], "val_metrics": []}

        for epoch in range(self.config.epochs):
            t0 = time.time()

            # Training epoch
            train_loss, train_metrics = self._train_epoch()

            # Validation epoch
            val_loss, val_metrics = self._validate()

            # Update learning rate
            self.scheduler.step()

            dt = time.time() - t0
            history["train_loss"].append(train_loss)
            history["val_loss"].append(val_loss)
            history["val_metrics"].append(val_metrics)

            # Log to W&B
            if self.config.use_wandb and self.wandb_run:
                import wandb

                log_dict = {
                    "epoch": epoch,
                    "train/loss": train_loss,
                    "val/loss": val_loss,
                    "lr": self.optimizer.param_groups[0]["lr"],
                    "epoch_time_s": dt,
                }
                for k, v in train_metrics.items():
                    log_dict[f"train/{k}"] = v
                for k, v in val_metrics.items():
                    log_dict[f"val/{k}"] = v
                wandb.log(log_dict)

            # Checkpoint and early stopping
            if val_loss < self.best_val_loss - self.config.min_delta:
                self.best_val_loss = val_loss
                self.best_epoch = epoch
                self.patience_counter = 0
                self._save_checkpoint(epoch, is_best=True)
            else:
                self.patience_counter += 1

            if epoch % 10 == 0 or self.patience_counter == 0:
                logger.info(
                    f"Epoch {epoch:4d}/{self.config.epochs}: "
                    f"train={train_loss:.6f}, val={val_loss:.6f}, "
                    f"lr={self.optimizer.param_groups[0]['lr']:.2e}, "
                    f"patience={self.patience_counter}/{self.config.patience}, "
                    f"dt={dt:.1f}s"
                )

            if self.patience_counter >= self.config.patience:
                logger.info(
                    f"Early stopping at epoch {epoch}. "
                    f"Best val loss: {self.best_val_loss:.6f} at epoch {self.best_epoch}"
                )
                break

        if self.config.use_wandb and self.wandb_run:
            import wandb
            wandb.finish()

        best_path = self.config.checkpoint_dir / "best_model.pt"
        history["best_model_path"] = str(best_path)
        history["best_val_loss"] = self.best_val_loss
        history["best_epoch"] = self.best_epoch

        return history

    def _train_epoch(self):
        """Run one training epoch."""
        self.model.train()
        total_loss = 0.0
        total_metrics = {}
        n_batches = 0

        for batch in self.train_loader:
            T_profile = batch["T_profile"].to(self.device)
            params = batch["params"].to(self.device)
            target_source = batch["source_term"].to(self.device)
            target_flux = batch["flux_up"].to(self.device).unsqueeze(1)

            # Forward pass
            pred_source, pred_flux = self.model(T_profile, params)

            # Compute loss
            loss, metrics = self.loss_fn(
                pred_source, pred_flux, target_source, target_flux,
                T_profile=T_profile, F_values=params[:, 1]
            )

            # Backward pass
            self.optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
            self.optimizer.step()

            total_loss += loss.item()
            for k, v in metrics.items():
                total_metrics[k] = total_metrics.get(k, 0.0) + v
            n_batches += 1

        avg_loss = total_loss / max(n_batches, 1)
        avg_metrics = {k: v / max(n_batches, 1) for k, v in total_metrics.items()}

        return avg_loss, avg_metrics

    @torch.no_grad()
    def _validate(self):
        """Run validation."""
        self.model.eval()
        total_loss = 0.0
        total_metrics = {}
        n_batches = 0

        for batch in self.val_loader:
            T_profile = batch["T_profile"].to(self.device)
            params = batch["params"].to(self.device)
            target_source = batch["source_term"].to(self.device)
            target_flux = batch["flux_up"].to(self.device).unsqueeze(1)

            pred_source, pred_flux = self.model(T_profile, params)

            loss, metrics = self.loss_fn(
                pred_source, pred_flux, target_source, target_flux,
                T_profile=T_profile, F_values=params[:, 1]
            )

            total_loss += loss.item()
            for k, v in metrics.items():
                total_metrics[k] = total_metrics.get(k, 0.0) + v
            n_batches += 1

        avg_loss = total_loss / max(n_batches, 1)
        avg_metrics = {k: v / max(n_batches, 1) for k, v in total_metrics.items()}

        return avg_loss, avg_metrics

    def _save_checkpoint(self, epoch: int, is_best: bool = False):
        """Save model checkpoint."""
        checkpoint = {
            "epoch": epoch,
            "model_state_dict": self.model.state_dict(),
            "optimizer_state_dict": self.optimizer.state_dict(),
            "scheduler_state_dict": self.scheduler.state_dict(),
            "val_loss": self.best_val_loss,
        }

        if is_best:
            path = self.config.checkpoint_dir / "best_model.pt"
        else:
            path = self.config.checkpoint_dir / f"checkpoint_epoch_{epoch}.pt"

        torch.save(checkpoint, path)

    def load_checkpoint(self, path: str):
        """Load model from checkpoint."""
        checkpoint = torch.load(path, map_location=self.device)
        self.model.load_state_dict(checkpoint["model_state_dict"])
        self.optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
        self.scheduler.load_state_dict(checkpoint["scheduler_state_dict"])
        self.best_val_loss = checkpoint["val_loss"]
        logger.info(f"Loaded checkpoint from {path} (epoch {checkpoint['epoch']})")
