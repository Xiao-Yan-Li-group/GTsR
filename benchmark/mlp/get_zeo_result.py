#!/usr/bin/env python3

import os
import glob
import csv
from tqdm import tqdm


path = "/mnt/d/Project/GTSR/benchmark/stability/mlp/zeopp"
OUT = f"/mnt/d/Project/GTSR/benchmark/stability/mlp/results_zeopp.csv"
rows = []

for res_file in tqdm(sorted(glob.glob(os.path.join(path, f"*.res")))):
    base = res_file[:-4]
    name = os.path.basename(base)
    volpo_file = base + ".volpo"

    with open(res_file) as f:
        line = f.readline().split()
        LCD, PLD, LFPD = float(line[1]), float(line[2]), float(line[3])

    with open(volpo_file) as f:
        for i, row in enumerate(f):
            if i == 0:
                volume  = float(row.split('Unitcell_volume:')[1].split()[0])
                density = float(row.split('Density:')[1].split()[0])
                VF      = float(row.split('POAV_Volume_fraction:')[1].split()[0])
                PV      = float(row.split('POAV_cm^3/g:')[1].split()[0])
                break

    rows.append({
        "name":    name,
        "LCD":     LCD,
        "PLD":     PLD,
        "LFPD":    LFPD,
        "volume":  volume,
        "density": density,
        "VF":      VF,
        "PV":      PV,
    })
fieldnames = ["name", "LCD", "PLD", "LFPD", "volume", "density", "VF", "PV"]
with open(OUT, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)