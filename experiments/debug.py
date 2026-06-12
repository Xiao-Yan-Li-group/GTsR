# from src.cif_utils import label2cif
# import numpy as np
# label = np.load("./webapp/logs/HALHEL_R.npy")
# label2cif("./webapp/logs/HALHEL_R_pacman.cif", label, "./")

# from src.cif_utils import get_sol_smi
# import numpy as np
# smis = get_sol_info("./webapp/logs/predictions/HALHEL_R_pacman_solvent.cif")
# print(smis)

# from gtsr import GTsRunner
# runner = GTsRunner(checkpoint="all")
# result = runner.clean(
#     cif="./webapp/logs/ABEXEM.cif",
#     output="./webapp/logs/",
#     threshold=0.5,
# )
# solvent_smiles = result["solvent_smiles"]
# # probabilities = result["probabilities"]
# print(solvent_smiles)

# from gtsr import GTsRunner
# runner = GTsRunner(checkpoint="stability")
# score = runner.stability(cif="./webapp/logs/ABEXOW_clean.cif")
# print(score)

# from pymatgen.core import Structure
# from pymatgen.io.cif import CifWriter
# struct = Structure.from_file("./webapp/logs/ABEXEM.cif")
# CifWriter(struct).write_file("./webapp/logs/ABEXEM_rw.cif")

import os, sys, time
t_start = time.time()
from tqdm import tqdm
from glob import glob
from gtsr import GTsRunner
cif_folder = sys.argv[1]
save_folder = sys.argv[2]
os.makedirs(save_folder, exist_ok=True)
model = sys.argv[3]
device = sys.argv[4]
fails = []
runner = GTsRunner(checkpoint=model, device=device)
for cif_path in tqdm(glob(cif_folder+"/*.cif")):

        result = runner.clean(
            cif=cif_path,
            output=save_folder,
            threshold=0.5,
        )
        solvent_smiles = result["solvent_smiles"]
        print(solvent_smiles)

with open(save_folder+"/gtsr_"+model+"_"+device+"_gpu.txt", "w") as f:
    f.write(f"time: {time.time() - t_start}\n")
    for fail in fails:
        f.write(fail+"\n")