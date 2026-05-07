from __future__ import annotations

import shutil
from pathlib import Path

import torch


class AverageMeter:
    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0.0
        self.avg = 0.0
        self.sum = 0.0
        self.count = 0

    def update(self, val, n=1):
        self.val = float(val)
        self.sum += float(val) * n
        self.count += n
        self.avg = self.sum / max(self.count, 1)


def save_checkpoint(state, is_best: bool, checkpoint_path: str | Path, best_path: str | Path):
    checkpoint_path = Path(checkpoint_path)
    best_path = Path(best_path)
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    best_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(state, checkpoint_path)
    if is_best:
        shutil.copyfile(checkpoint_path, best_path)


def binary_metrics_from_logits(logits: torch.Tensor, target: torch.Tensor, threshold: float = 0.5):
    probs = torch.sigmoid(logits.detach())
    pred = probs >= threshold
    truth = target.detach() >= 0.5
    tp = torch.logical_and(pred, truth).sum().item()
    tn = torch.logical_and(~pred, ~truth).sum().item()
    fp = torch.logical_and(pred, ~truth).sum().item()
    fn = torch.logical_and(~pred, truth).sum().item()

    total = max(tp + tn + fp + fn, 1)
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-12)
    return {
        "accuracy": (tp + tn) / total,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
    }

