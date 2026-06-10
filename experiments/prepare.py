from __future__ import annotations
import json
import os
from glob import glob
from src.cif_utils import cif2graph, cif2pos, make_label

cif_path = "../benchmark/dataset/clean_dataset/"

graph_path = "../benchmark/dataset/graph"
pos_path = "../benchmark/dataset/pos"

label_path_free = "../benchmark/dataset/label/free"
label_path_all = "../benchmark/dataset/label/all"

os.makedirs(graph_path, exist_ok=True)
os.makedirs(pos_path, exist_ok=True)
os.makedirs(label_path_free, exist_ok=True)
os.makedirs(label_path_all, exist_ok=True)

import pandas as pd
import numpy as np
from tqdm import tqdm


data = pd.read_csv("../benchmark/dataset/dataset_labels.csv")

### label generateing

def solvent_check(label, tol=0.8):
    return np.sum(label == 1) / len(label) > tol

input_list = []
r2f_fail = []
f2a_fail = []

for i, row in tqdm(data.iterrows(), total=len(data)):
    name = row["refcode"]
    
    if row["r2f"] == 1:
        try:
            label_R = make_label(os.path.join(cif_path, f"{name}_R.cif"), os.path.join(cif_path, f"{name}_F.cif"))
            label_F = make_label(os.path.join(cif_path, f"{name}_F.cif"), os.path.join(cif_path, f"{name}_F.cif"))
            if not solvent_check(label_R):
                np.save(os.path.join(label_path_free, f"{name}_R.npy"), label_R)
                input_list.append(f"{name}_R")
            np.save(os.path.join(label_path_free, f"{name}_F.npy"), label_F)
            input_list.append(f"{name}_F")
        except:
            r2f_fail.append(name)
        
        if row["f2a"] == 1:
            try:
                label_R = make_label(os.path.join(cif_path, f"{name}_R.cif"), os.path.join(cif_path, f"{name}_A.cif"))
                label_F = make_label(os.path.join(cif_path, f"{name}_F.cif"), os.path.join(cif_path, f"{name}_A.cif"))
                label_A = make_label(os.path.join(cif_path, f"{name}_A.cif"), os.path.join(cif_path, f"{name}_A.cif"))
                if not solvent_check(label_R):
                    np.save(os.path.join(label_path_all, f"{name}_R.npy"), label_R)
                    input_list.append(f"{name}_R")
                if not solvent_check(label_F):
                    np.save(os.path.join(label_path_all, f"{name}_F.npy"), label_F)
                    input_list.append(f"{name}_F")
                np.save(os.path.join(label_path_all, f"{name}_A.npy"), label_A)
                np.save(os.path.join(label_path_free, f"{name}_A.npy"), label_A)
                input_list.append(f"{name}_A")
            except:
                f2a_fail.append(name)
        else:
            try:
                label_R = make_label(os.path.join(cif_path, f"{name}_R.cif"), os.path.join(cif_path, f"{name}_F.cif"))
                label_F = make_label(os.path.join(cif_path, f"{name}_F.cif"), os.path.join(cif_path, f"{name}_F.cif"))
                if not solvent_check(label_R):
                    np.save(os.path.join(label_path_all, f"{name}_R.npy"), label_R)
                    input_list.append(f"{name}_R")
                np.save(os.path.join(label_path_all, f"{name}_F.npy"), label_F)
                input_list.append(f"{name}_F")
            except:
                f2a_fail.append(name)
    else:
        try:
            label_R = make_label(os.path.join(cif_path, f"{name}_R.cif"), os.path.join(cif_path, f"{name}_R.cif"))
            np.save(os.path.join(label_path_free, f"{name}_R.npy"), label_R)
            input_list.append(f"{name}_R")
        except:
            f2a_fail.append(name)
        
        if row["f2a"] == 1:
            try:
                label_R = make_label(os.path.join(cif_path, f"{name}_R.cif"), os.path.join(cif_path, f"{name}_A.cif"))
                label_A = make_label(os.path.join(cif_path, f"{name}_A.cif"), os.path.join(cif_path, f"{name}_A.cif"))
                if not solvent_check(label_R):
                    np.save(os.path.join(label_path_all, f"{name}_R.npy"), label_R)
                    input_list.append(f"{name}_R")
                np.save(os.path.join(label_path_all, f"{name}_A.npy"), label_A)
                np.save(os.path.join(label_path_free, f"{name}_A.npy"), label_A)
                input_list.append(f"{name}_A")
            except:
                f2a_fail.append(name)
        else:
            try:
                label_R = make_label(os.path.join(cif_path, f"{name}_R.cif"), os.path.join(cif_path, f"{name}_R.cif"))
                np.save(os.path.join(label_path_all, f"{name}_R.npy"), label_R)
                input_list.append(f"{name}_R")
            except:
                f2a_fail.append(name)

with open("../benchmark/dataset/label_fail_list.txt", "w") as f:
    f.write("r2f_fail:\n")
    for name in r2f_fail:
        f.write(f"{name}\n")
    f.write("\nf2a_fail:\n")
    for name in f2a_fail:
        f.write(f"{name}\n")

### graph & pos generateing

input_fail = []
unique_inputs = set(input_list)

for name in tqdm(unique_inputs):
    try:
        graph = cif2graph(os.path.join(cif_path, f"{name}.cif"))
        pos = cif2pos(os.path.join(cif_path, f"{name}.cif"))
        np.save(os.path.join(graph_path, f"{name}.npy"), graph)
        np.save(os.path.join(pos_path, f"{name}.npy"), pos)
    except:
        print(f"r: {name}")
        input_fail.append(name)
        continue

with open("../benchmark/dataset/input_fail_list.txt", "w") as f:
    f.write("input_fail:\n")
    for name in input_fail:
        f.write(f"{name}\n")


def split_data(split_dir="../benchmark/dataset/splits", seed=42, train_ratio=0.8, val_ratio=0.1):
    split_dir = os.path.abspath(split_dir)
    os.makedirs(split_dir, exist_ok=True)

    def extract_refcode(name):
        refcode = os.path.splitext(os.path.basename(name))[0]
        for suffix in ["_IF", "_IA", "_R", "_F", "_A"]:
            if refcode.endswith(suffix):
                refcode = refcode.replace(suffix, "")
                break
        return refcode

    def read_label_ids(label_dir):
        ids = sorted(os.path.splitext(os.path.basename(p))[0] for p in glob(os.path.join(label_dir, "*.npy")))
        return [name for name in ids if os.path.exists(os.path.join(graph_path, f"{name}.npy"))]

    def split_refcodes(refcodes):
        refcodes = sorted(set(refcodes))
        rng = np.random.default_rng(seed)
        shuffled = np.array(refcodes, dtype=object)
        rng.shuffle(shuffled)

        n_total = len(shuffled)
        n_train = int(n_total * train_ratio)
        n_val = int(n_total * val_ratio)

        return {
            "train": set(shuffled[:n_train]),
            "val": set(shuffled[n_train:n_train + n_val]),
            "test": set(shuffled[n_train + n_val:]),
        }

    free_ids = read_label_ids(label_path_free)
    coordinated_ids = read_label_ids(label_path_all)
    all_refcodes = [extract_refcode(name) for name in free_ids + coordinated_ids]
    refcode_splits = split_refcodes(all_refcodes)

    task_ids = {
        "free": free_ids,
        "coordinated": coordinated_ids,
    }
    summary_rows = []

    for task, ids in task_ids.items():
        for split_name, split_refcodes in refcode_splits.items():
            split_ids = [name for name in ids if extract_refcode(name) in split_refcodes]
            output_csv = os.path.join(split_dir, f"{task}_{split_name}.csv")
            pd.DataFrame({"cif_id": sorted(split_ids)}).to_csv(output_csv, index=False)
            summary_rows.append({
                "task": task,
                "split": split_name,
                "structures": len(split_ids),
                "refcodes": len(set(extract_refcode(name) for name in split_ids)),
            })

    pd.DataFrame(summary_rows).to_csv(os.path.join(split_dir, "split_summary.csv"), index=False)
    print(pd.DataFrame(summary_rows))


split_data()