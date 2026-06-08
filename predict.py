from __future__ import annotations

import argparse
import torch
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from src.data import GaussianDistance
from src.cif_utils import cif2graph, cif2pos, label2cif
from src.GCN import SolventAtomClassifier
from src.utils import load_checkpoint

def main() -> None:
    parser = argparse.ArgumentParser(description="Predict and remove solvent by GTsR.")
    parser.add_argument("--cif", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=0.5)
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    graph = cif2graph(args.cif)
    pos = cif2pos(args.cif)

    gdf = GaussianDistance(dmin=0, dmax=8, step=0.2)

    atom_fea = np.eye(119, dtype=np.float32)[np.array(graph["numbers"], dtype=int)]
    index1 = np.array(graph["index1"], dtype=np.int64)
    index2 = np.array(graph["index2"], dtype=np.int64)
    dij = np.array(graph["dij"], dtype=np.float32)
    nbr_fea = gdf.expand(dij).astype(np.float32)
    atom_fea_pos_tensor = torch.Tensor(np.concatenate([atom_fea, pos], axis=1))

    nbr_fea_tensor = torch.Tensor(nbr_fea)
    index1_tensor = torch.LongTensor(index1)
    index2_tensor = torch.LongTensor(index2)
    atom_fea_tensor = torch.LongTensor([0] * len(atom_fea))
    tensors = tuple(tensor.to(device) for tensor in (atom_fea_pos_tensor, nbr_fea_tensor, index1_tensor, index2_tensor, atom_fea_tensor))

    checkpoint = load_checkpoint(args.checkpoint,device=device)
    model_config = checkpoint["model_config"]
    model = SolventAtomClassifier(**model_config).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    with torch.no_grad():
        prob = torch.sigmoid(model(*tensors)).cpu().numpy()
    
    labels = (prob >= args.threshold).astype(int)
    label2cif(args.cif, labels, args.output)

if __name__ == "__main__":
    main()
