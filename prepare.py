from __future__ import annotations
import json
import os

from src.cif_utils import cif2graph, cif2pos, make_label

cif_path = "./benchmark/dataset/clean_dataset/"

graph_path = "./benchmark/dataset/graph"
pos_path = "./benchmark/dataset/pos"

label_path_free = "./benchmark/dataset/label/free"
label_path_all = "./benchmark/dataset/label/all"

os.makedirs(graph_path, exist_ok=True)
os.makedirs(pos_path, exist_ok=True)
os.makedirs(label_path_free, exist_ok=True)
os.makedirs(label_path_all, exist_ok=True)

import pandas as pd
import numpy as np
from tqdm import tqdm


data = pd.read_csv("./benchmark/dataset/dataset_labels.csv")

### label generateing

def solvent_check(label, tol=0.8):
    return np.sum(label == 1) / len(label) > tol

input_list = []
r2f_fail = []
f2a_fail = []

for i, row in tqdm(data.head(10).iterrows(), total=10):
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
            try:
                label_R = make_label(os.path.join(cif_path, f"{name}_R.cif"), os.path.join(cif_path, f"{name}_IF.cif"))
                label_F = make_label(os.path.join(cif_path, f"{name}_IF.cif"), os.path.join(cif_path, f"{name}_IF.cif"))
                if not solvent_check(label_R):
                    np.save(os.path.join(label_path_free, f"{name}_R.npy"), label_R)
                    input_list.append(f"{name}_R")
                np.save(os.path.join(label_path_free, f"{name}_IF.npy"), label_F)
                input_list.append(f"{name}_IF")
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
                try:
                    label_R = make_label(os.path.join(cif_path, f"{name}_R.cif"), os.path.join(cif_path, f"{name}_IA.cif"))
                    label_F = make_label(os.path.join(cif_path, f"{name}_F.cif"), os.path.join(cif_path, f"{name}_IA.cif"))
                    label_A = make_label(os.path.join(cif_path, f"{name}_IA.cif"), os.path.join(cif_path, f"{name}_IA.cif"))
                    if not solvent_check(label_R):
                        np.save(os.path.join(label_path_all, f"{name}_R.npy"), label_R)
                        input_list.append(f"{name}_R")
                    if not solvent_check(label_F):
                        np.save(os.path.join(label_path_all, f"{name}_IF.npy"), label_F)
                        input_list.append(f"{name}_IF")
                    np.save(os.path.join(label_path_all, f"{name}_IA.npy"), label_A)
                    input_list.append(f"{name}_IA")
                except:
                    f2a_fail.append(name)
        else:
            try:
                
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
                try:
                    label_R = make_label(os.path.join(cif_path, f"{name}_R.cif"), os.path.join(cif_path, f"{name}_IA.cif"))
                    label_A = make_label(os.path.join(cif_path, f"{name}_IA.cif"), os.path.join(cif_path, f"{name}_IA.cif"))
                    if not solvent_check(label_R):
                        np.save(os.path.join(label_path_all, f"{name}_R.npy"), label_R)
                        input_list.append(f"{name}_R")
                    np.save(os.path.join(label_path_all, f"{name}_IA.npy"), label_A)
                    np.save(os.path.join(label_path_free, f"{name}_IA.npy"), label_A)
                    input_list.append(f"{name}_IA")
                except:
                    f2a_fail.append(name)
        else:
            try:
                np.save(os.path.join(label_path_all, f"{name}_R.npy"), label_R)
                input_list.append(f"{name}_R")
            except:
                f2a_fail.append(name)

with open("./benchmark/dataset/label_fail_list.txt", "w") as f:
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

with open("./benchmark/dataset/input_fail_list.txt", "w") as f:
    f.write("input_fail:\n")
    for name in input_fail:
        f.write(f"{name}\n")