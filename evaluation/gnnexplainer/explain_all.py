#!/usr/bin/env python3
import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import torch
from ase.data import chemical_symbols
from ase.io import read, write

HERE = Path(__file__).resolve().parent
GTSR = HERE.parents[1]
sys.path.insert(0, str(GTSR))
from src.GCN import SolventAtomClassifier
from src.cif_utils import cif2graph
from src.data import GaussianDistance


def scatter_mean(src, index, n, count):
    out = torch.zeros((n, src.shape[1]), dtype=src.dtype, device=src.device)
    return out.scatter_add(0, index.view(-1, 1).expand_as(src), src / count)


def conv_masked(layer, atom_in_fea, nbr_fea, i1, i2, crys, emask):
    n_atom, width = atom_in_fea.shape
    counts = torch.zeros(n_atom, dtype=atom_in_fea.dtype, device=atom_in_fea.device).scatter_add(
        0, i1, torch.ones_like(i1, dtype=atom_in_fea.dtype))
    cnt = counts[i1].clamp(min=1.0).view(-1, 1)
    m = emask.view(-1, 1)
    ek = layer.phi_e(torch.cat([atom_in_fea[i1], atom_in_fea[i2], nbr_fea], dim=1))
    rho_e_v = scatter_mean(ek * m, i1, n_atom, cnt)
    ek = nbr_fea + ek
    ek_sum = scatter_mean(ek * m, i1, n_atom, cnt)
    vi = atom_in_fea + layer.norm_v(layer.phi_v(torch.cat([atom_in_fea, rho_e_v], dim=1)))
    atom_nbr_fea = torch.cat([vi, ek_sum], dim=1)
    n_crys = int(crys.max()) + 1
    global_fea = torch.zeros((n_crys, 2 * width), dtype=vi.dtype, device=vi.device).scatter_add(
        0, crys.view(-1, 1).expand_as(atom_nbr_fea), atom_nbr_fea)
    return ek, vi, global_fea, atom_nbr_fea


def forward_masked(model, atom_fea, nbr_fea, i1, i2, crys, emask):
    atom_fea = model.node_embedding(atom_fea)
    nbr_fea = model.edge_embedding(nbr_fea)
    atom_nbr_fea = torch.cat([atom_fea, torch.zeros_like(atom_fea)], dim=1)
    global_fea = None
    for conv in model.convs:
        nbr_fea, atom_fea, global_fea, atom_nbr_fea = conv_masked(conv, atom_fea, nbr_fea, i1, i2, crys, emask)
    n_crys = int(crys.max()) + 1
    counts = torch.zeros((n_crys, 1), dtype=atom_fea.dtype, device=atom_fea.device).scatter_add(
        0, crys.view(-1, 1), torch.ones((crys.numel(), 1), dtype=atom_fea.dtype, device=atom_fea.device))
    global_fea = global_fea / counts.clamp(min=1.0)
    local = model.local_embedding(atom_nbr_fea)
    glob = model.global_embedding(global_fea[crys])
    return model.classifier(torch.cat([local, glob], dim=-1)).squeeze(-1)


def load_model(path, device):
    ck = torch.load(path, map_location="cpu", weights_only=False)
    model = SolventAtomClassifier(**ck["model_config"]).to(device)
    model.load_state_dict(ck["state_dict"])
    model.eval()
    for p in model.parameters():
        p.requires_grad_(False)
    return model, ck


def build_inputs(cif, ck):
    g = cif2graph(cif, radius=float(ck.get("radius", 8.0)))
    pos = read(cif).get_scaled_positions().astype(np.float32)
    gdf = GaussianDistance(dmin=float(ck.get("dmin", 0.0)), dmax=float(ck.get("radius", 8.0)),
                           step=float(ck.get("step", 0.2)))
    numbers = np.array(g["numbers"], dtype=int)
    atom_fea = np.concatenate([np.eye(119, dtype=np.float32)[numbers], pos], axis=1)
    dij = np.array(g["dij"], dtype=np.float32)
    return (numbers, np.array(g["index1"], dtype=np.int64), np.array(g["index2"], dtype=np.int64),
            dij, torch.tensor(atom_fea), torch.tensor(gdf.expand(dij), dtype=torch.float32))


def batched(atom_fea, nbr_fea, i1, i2, n_copies, device):
    n, e = atom_fea.shape[0], i1.numel()
    off = (torch.arange(n_copies) * n).repeat_interleave(e)
    return (atom_fea.repeat(n_copies, 1).to(device), nbr_fea.repeat(n_copies, 1).to(device),
            (i1.repeat(n_copies) + off).to(device), (i2.repeat(n_copies) + off).to(device),
            torch.arange(n_copies).repeat_interleave(n).to(device))


def explain_batch(model, inputs, n_atom, n_edge, targets, cls, epochs, lr, size_coef, ent_coef):
    atom_fea, nbr_fea, i1, i2, crys = inputs
    B = len(targets)
    dev = atom_fea.device
    std = torch.nn.init.calculate_gain("relu") * np.sqrt(2.0 / (2 * n_atom))
    param = torch.nn.Parameter(torch.randn(B, n_edge, device=dev) * std + 1.0)
    opt = torch.optim.Adam([param], lr=lr)
    tgt = torch.arange(B, device=dev) * n_atom + torch.as_tensor(targets, device=dev)
    want = torch.as_tensor(cls, dtype=torch.float32, device=dev)
    for _ in range(epochs):
        opt.zero_grad()
        m = torch.sigmoid(param)
        p = torch.sigmoid(forward_masked(model, atom_fea, nbr_fea, i1, i2, crys, m.reshape(-1))[tgt])
        pred = -(want * torch.log(p + 1e-12) + (1 - want) * torch.log(1 - p + 1e-12))
        ent = -(m * torch.log(m + 1e-12) + (1 - m) * torch.log(1 - m + 1e-12)).mean(1)
        (pred + size_coef * m.sum(1) + ent_coef * ent).sum().backward()
        opt.step()
    with torch.no_grad():
        m = torch.sigmoid(param)
        p_soft = torch.sigmoid(forward_masked(model, atom_fea, nbr_fea, i1, i2, crys, m.reshape(-1))[tgt])
        p_hard = torch.sigmoid(forward_masked(model, atom_fea, nbr_fea, i1, i2, crys,
                                              (m > 0.5).float().reshape(-1))[tgt])
    return m.cpu().numpy(), p_soft.cpu().numpy(), p_hard.cpu().numpy()


def run_structure(cif, name, model, ck, out_dir, device, batch, epochs, lr, size_coef, ent_coef):
    cid = Path(cif).stem
    out = out_dir / name / cid
    if (out / "atoms.csv").exists():
        print(f"[skip] {name} {cid}: {out / 'atoms.csv'} exists")
        return
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    thr = float(ck.get("threshold", 0.5))
    numbers, i1n, i2n, dij, atom_fea, nbr_fea = build_inputs(cif, ck)
    n_atom, n_edge = len(numbers), len(i1n)
    i1, i2 = torch.tensor(i1n), torch.tensor(i2n)

    one = batched(atom_fea, nbr_fea, i1, i2, 1, device)
    with torch.no_grad():
        logit = model(*one)
        prob = torch.sigmoid(logit).cpu().numpy()
        same = torch.allclose(logit, forward_masked(model, *one, torch.ones(n_edge, device=device)), atol=1e-4)
        p_noedge = torch.sigmoid(forward_masked(model, *one, torch.zeros(n_edge, device=device))).cpu().numpy()
    print(f"{name} {cid}: {n_atom} atoms, {n_edge} edges ({int((dij < 1e-6).sum())} self-edges), "
          f"predicted solvent {int((prob >= thr).sum())} (threshold {thr}), "
          f"masked forward == model: {same}", flush=True)

    cls = (prob >= 0.5).astype(np.float32)
    masks = np.zeros((n_atom, n_edge), dtype=np.float16)
    p_soft = np.zeros(n_atom, dtype=np.float32)
    p_hard = np.zeros(n_atom, dtype=np.float32)
    for start in range(0, n_atom, batch):
        tg = list(range(start, min(start + batch, n_atom)))
        inp = batched(atom_fea, nbr_fea, i1, i2, len(tg), device)
        m, ps, ph = explain_batch(model, inp, n_atom, n_edge, tg, cls[tg], epochs, lr, size_coef, ent_coef)
        masks[tg], p_soft[tg], p_hard[tg] = m, ps, ph
        print(f"  atoms {tg[0]}-{tg[-1]} done ({time.time() - t0:.0f} s)", flush=True)

    mf = masks.astype(np.float32)
    excess = np.clip(mf - np.median(mf, axis=1, keepdims=True), 0, None)
    excess[:, dij < 1e-6] = 0.0
    pred = prob >= thr
    influence = np.zeros(n_atom)
    influence_solv = np.zeros(n_atom)
    np.add.at(influence, i2n, excess.sum(0))
    np.add.at(influence_solv, i2n, excess[pred].sum(0))
    kept = (mf > 0.5).sum(1)
    complete = (p_hard >= 0.5) == (prob >= 0.5)

    cols = {"atom": np.arange(n_atom), "element": [chemical_symbols[int(z)] for z in numbers],
            "probability": prob, "pred": pred.astype(int), "p_no_edges": p_noedge,
            "edge_dependence": prob - p_noedge, "p_soft_mask": p_soft, "p_kept_edges": p_hard,
            "n_edges_kept": kept, "explanation_complete": complete.astype(int),
            "influence": influence, "influence_solvent": influence_solv}
    with open(out / "atoms.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for k in range(n_atom):
            w.writerow([f"{cols[c][k]:.5f}" if isinstance(cols[c][k], (float, np.floating)) else cols[c][k]
                        for c in cols])
    np.savez_compressed(out / "masks.npz", masks=masks, index1=i1n, index2=i2n, dij=dij,
                        numbers=numbers, probability=prob, threshold=thr)
    atoms = read(cif)
    if len(atoms) == n_atom and list(atoms.numbers) == list(numbers):
        for c in cols:
            if c not in ("atom", "element"):
                atoms.new_array(c, np.asarray(cols[c], dtype=float))
        write(out / f"{cid}.extxyz", atoms, format="extxyz")
    else:
        print("  note: ASE atom order differs from the graph; extxyz not written (use atoms.csv)")
    print(f"  -> {out}  ({time.time() - t0:.0f} s; complete explanations "
          f"{int(complete.sum())}/{n_atom})", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif", nargs="+", required=True, help="CIF file(s)")
    ap.add_argument("--model", required=True, choices=["free", "all"],
                    help="free: free-solvent model, all: all-solvent model; also the output subfolder")
    ap.add_argument("--checkpoint", help="model file (default: GTsR/ckpt/<model>_best.pth)")
    ap.add_argument("--out", default=str(HERE / "output"))
    ap.add_argument("--batch", type=int, default=32, help="atoms explained per optimisation")
    ap.add_argument("--epochs", type=int, default=300)
    ap.add_argument("--lr", type=float, default=0.01)
    ap.add_argument("--size-coef", type=float, default=0.005)
    ap.add_argument("--ent-coef", type=float, default=1.0)
    ap.add_argument("--device", default="auto", help="auto, cpu, cuda, cuda:1, ...")
    args = ap.parse_args()

    device = torch.device(("cuda" if torch.cuda.is_available() else "cpu") if args.device == "auto"
                          else args.device)
    checkpoint = Path(args.checkpoint or GTSR / "ckpt" / f"{args.model}_best.pth")
    print(f"device: {device} | model: {args.model} ({checkpoint})", flush=True)
    model, ck = load_model(checkpoint, device)
    for cif in args.cif:
        run_structure(cif, args.model, model, ck, Path(args.out), device, args.batch, args.epochs,
                      args.lr, args.size_coef, args.ent_coef)


if __name__ == "__main__":
    main()
