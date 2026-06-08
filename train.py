from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.optim.lr_scheduler import ExponentialLR

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from src.GCN import SolventAtomClassifier
from src.data import SolventCIFData, collate_pool, get_data_loader
from src.utils import AverageMeter, metric_functions, save_checkpoint, load_checkpoint

_ROC_TOOLS = None

def compute_pos_weight(label_dir: Path, ids: list[str]) -> float:
    positive = 0.0
    total = 0.0
    for cif_id in ids:
        labels = np.load(label_dir / f"{cif_id}.npy")
        positive += float(labels.sum())
        total += float(labels.size)
    negative = max(total - positive, 0.0)
    if positive <= 0:
        return 1.0
    return max(negative / positive, 1.0)


def move_model_input(batch_input, device: torch.device):
    atom_fea = torch.cat((batch_input[0], batch_input[6]), dim=-1).to(device)
    return (
        atom_fea,
        batch_input[1].to(device),
        batch_input[2].to(device),
        batch_input[3].to(device),
        batch_input[4].to(device),
    )


def run_epoch(loader, model, criterion, optimizer, device, threshold: float, train: bool, collect_scores: bool = False):
    losses = AverageMeter()
    totals = {"tp": 0, "tn": 0, "fp": 0, "fn": 0}
    y_true_parts = []
    y_score_parts = []

    model.train(train)
    for batch_input, target, _ in loader:
        target = target.to(device)
        model_input = move_model_input(batch_input, device)
        with torch.set_grad_enabled(train):
            output = model(*model_input)
            loss = criterion(output, target)
            if train:
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

        metrics = metric_functions(output, target, threshold=threshold)
        for key in totals:
            totals[key] += metrics[key]
        losses.update(loss.item(), target.numel())
        if collect_scores:
            y_true_parts.append(target.detach().cpu().numpy())
            y_score_parts.append(torch.sigmoid(output.detach()).cpu().numpy())

    tp, tn, fp, fn = totals["tp"], totals["tn"], totals["fp"], totals["fn"]
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-12)
    accuracy = (tp + tn) / max(tp + tn + fp + fn, 1)
    result = {"loss": losses.avg, "accuracy": accuracy, "precision": precision, "recall": recall, "f1": f1, **totals}
    if collect_scores:
        result["y_true"] = np.concatenate(y_true_parts)
        result["y_score"] = np.concatenate(y_score_parts)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Train a GTsR3 atom-level solvent classifier.")
    parser.add_argument("--task", choices=("free", "all"), required=True)
    parser.add_argument("--data-dir", type=Path, default=SCRIPT_DIR / "data")
    parser.add_argument("--model-dir", type=Path, default=SCRIPT_DIR / "pth")
    parser.add_argument("--radius", type=float, default=8.0)
    parser.add_argument("--dmin", type=float, default=0.0)
    parser.add_argument("--step", type=float, default=0.2)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--atom-fea-len", type=int, default=128)
    parser.add_argument("--n-conv", type=int, default=6)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--lr", type=float, default=5e-4)
    parser.add_argument("--lr-decay-rate", type=float, default=0.99)
    parser.add_argument("--weight-decay", type=float, default=0.0)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--resume", type=Path, default=None)
    args = parser.parse_args()

    task_data_dir = args.data_dir / args.task
    graph_dir = task_data_dir / "json"
    pos_dir = task_data_dir / "npy" / "pos"
    label_dir = task_data_dir / "npy" / "label"
    split_dir = args.data_dir / "splits"

    train_dataset = SolventCIFData(
        graph_dir, pos_dir, label_dir, split_dir / f"{args.task}_train.csv", args.radius, args.dmin, args.step
    )
    val_dataset = SolventCIFData(
        graph_dir, pos_dir, label_dir, split_dir / f"{args.task}_val.csv", args.radius, args.dmin, args.step
    )
    test_dataset = SolventCIFData(
        graph_dir, pos_dir, label_dir, split_dir / f"{args.task}_test.csv", args.radius, args.dmin, args.step
    )
    train_loader = get_data_loader(train_dataset, collate_pool, args.batch_size, args.num_workers, test=False)
    val_loader = get_data_loader(val_dataset, collate_pool, args.batch_size, args.num_workers, test=True)
    test_loader = get_data_loader(test_dataset, collate_pool, args.batch_size, args.num_workers, test=True)

    sample_input, _, _ = collate_pool([train_dataset[0]])
    orig_atom_fea_len = sample_input[0].shape[-1] + sample_input[6].shape[-1]
    nbr_fea_len = sample_input[1].shape[-1]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = SolventAtomClassifier(
        orig_atom_fea_len=orig_atom_fea_len,
        nbr_fea_len=nbr_fea_len,
        atom_fea_len=args.atom_fea_len,
        n_conv=args.n_conv,
        dropout=args.dropout,
    ).to(device)

    pos_weight = compute_pos_weight(label_dir, train_dataset.id_prop_data)
    criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([pos_weight], device=device))
    optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    scheduler = ExponentialLR(optimizer, gamma=args.lr_decay_rate)

    start_epoch = 0
    best_f1 = -1.0
    if args.resume is not None:
        checkpoint = load_checkpoint(args.resume, device)
        model.load_state_dict(checkpoint["state_dict"])
        optimizer.load_state_dict(checkpoint["optimizer"])
        start_epoch = int(checkpoint["epoch"]) + 1
        best_f1 = float(checkpoint.get("best_f1", best_f1))

    task_dir = args.model_dir / args.task
    checkpoint_path = task_dir / "checkpoint.pth"
    best_path = task_dir / "best.pth"
    roc_dir = args.roc_dir if args.roc_dir is not None else task_dir / "roc"

    print(f"task: {args.task}")
    print(f"train/val/test: {len(train_dataset)}/{len(val_dataset)}/{len(test_dataset)}")
    print(f"device: {device}")
    print(f"pos_weight: {pos_weight:.4f}")

    model_config = {
        "orig_atom_fea_len": orig_atom_fea_len,
        "nbr_fea_len": nbr_fea_len,
        "atom_fea_len": args.atom_fea_len,
        "n_conv": args.n_conv,
        "dropout": args.dropout,
    }

    for epoch in range(start_epoch, args.epochs):
        train_metrics = run_epoch(train_loader, model, criterion, optimizer, device, args.threshold, train=True)
        val_metrics = run_epoch(
            val_loader,
            model,
            criterion,
            optimizer,
            device,
            args.threshold,
            train=False,
            collect_scores=args.save_roc,
        )
        scheduler.step()

        is_best = val_metrics["f1"] > best_f1
        best_f1 = max(best_f1, val_metrics["f1"])
        save_checkpoint(
            {
                "epoch": epoch,
                "state_dict": model.state_dict(),
                "optimizer": optimizer.state_dict(),
                "best_f1": best_f1,
                "task": args.task,
                "threshold": args.threshold,
                "radius": args.radius,
                "dmin": args.dmin,
                "step": args.step,
                "model_config": model_config,
                "args": vars(args),
            },
            is_best,
            checkpoint_path,
            best_path,
        )

        auc_text = f" auc {val_metrics['auc']:.4f}" if args.save_roc else ""
        print(
            f"epoch {epoch:03d} train loss {train_metrics['loss']:.4f} f1 {train_metrics['f1']:.4f} "
            f"val loss {val_metrics['loss']:.4f} f1 {val_metrics['f1']:.4f}{auc_text} "
            f"p/r {val_metrics['precision']:.4f}/{val_metrics['recall']:.4f}"
        )

    best_checkpoint = load_checkpoint(best_path, device)
    model.load_state_dict(best_checkpoint["state_dict"])
    test_metrics = run_epoch(
        test_loader,
        model,
        criterion,
        optimizer,
        device,
        args.threshold,
        train=False,
        collect_scores=args.save_roc,
    )
    
    test_auc_text = f" auc {test_metrics['auc']:.4f}" if args.save_roc else ""
    print(
        f"test loss {test_metrics['loss']:.4f} f1 {test_metrics['f1']:.4f}{test_auc_text} "
        f"accuracy {test_metrics['accuracy']:.4f} p/r {test_metrics['precision']:.4f}/{test_metrics['recall']:.4f}"
    )
    print(f"best model: {best_path}")
    if args.save_roc:
        print(f"roc curves: {roc_dir}")


if __name__ == "__main__":
    main()
