import re
import time
from clean import free_clean, all_clean
from utils import METAL
from ase.io import write
import os, sys, json
from tqdm import tqdm
from glob import glob
from pathlib import Path


metal_list = [element for element, is_metal in METAL.items() if is_metal]

def run_fsr(cif_path, save_folder):
    skin=0.25
    while True:
        printed_formulas, structure = free_clean(cif_path, skin)
        has_metals = False
        for e_s in printed_formulas:
            split_formula = re.findall(r'([A-Z][a-z]?)(\d*)', e_s)
            elements = [match[0] for match in split_formula]
            if any(e in metal_list for e in elements):
                has_metals = True
                skin += 0.05
        if not has_metals:
            new_fn = Path(cif_path).stem + "_core_F.cif"
            write(os.path.join(save_folder, new_fn), structure)
            break
    return printed_formulas

def run_asr(cif_path, save_folder):
    skin=0.25
    while True:
        printed_formulas, structure = all_clean(cif_path, skin)
        has_metals = False
        for e_s in printed_formulas:
            split_formula = re.findall(r'([A-Z][a-z]?)(\d*)', e_s)
            elements = [match[0] for match in split_formula]
            if any(e in metal_list for e in elements):
                has_metals = True
                skin += 0.05
        if not has_metals:
            new_fn = Path(cif_path).stem + "_core_A.cif"
            write(os.path.join(save_folder, new_fn), structure)
            break
    return printed_formulas


cif_folder = sys.argv[1]
save_folder = sys.argv[2]
save_json = sys.argv[3]

results = {}
t_start = time.time()
for cif_path in tqdm(glob(cif_folder+"/*.cif")):
    printed_formulas_fsr=run_fsr(cif_path, save_folder)
    printed_formulas_asr=run_asr(cif_path, save_folder)
    results[Path(cif_path).stem]={}
    results[Path(cif_path).stem]["FSR"] = printed_formulas_fsr
    results[Path(cif_path).stem]["ASR"] = printed_formulas_asr
results["time"] = time.time() - t_start

with open(save_json, "w") as f:
    json.dump(results, f, indent=2)