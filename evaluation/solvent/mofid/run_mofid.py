from mofid.id_constructor import extract_fragments
import os, sys, shutil, time
from tqdm import tqdm
from glob import glob
from pathlib import Path


def get_node_linker_files(structure, prefix):
    os.makedirs(prefix, exist_ok=True)
    extract_fragments(structure, prefix)

cif_folder = sys.argv[1]
save_folder = sys.argv[2]
os.makedirs(save_folder, exist_ok=True)

model = "StandardIsolated" # "AllNode", "MetalOxo", "SingleNode"
mofid_cifs = ["mof_asr.cif", "mof_fsr.cif"]
types = ["A", "F"]

fails = []
t_start = time.time()
for cif_path in tqdm(glob(cif_folder+"/*.cif")):
    new_path=Path(cif_path).stem
    if all(os.path.exists(os.path.join(save_folder, new_path+"_mofid_"+model+"_"+t+".cif")) for t in types):
        continue
    try:
        get_node_linker_files(cif_path, os.path.join(new_path))
        for i, mofid_cif in enumerate(mofid_cifs):
            shutil.move(os.path.join(new_path, model, mofid_cif),
                        os.path.join(save_folder, new_path+"_mofid_"+model+"_"+types[i]+".cif"))
        shutil.rmtree(new_path)
    except:
        shutil.rmtree(new_path)
        fails.append(new_path)

with open("./failed_list.txt", "w") as f:
    f.write(f"time: {time.time() - t_start}\n")
    for fail in fails:
        f.write(fail+"\n")