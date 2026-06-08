from src.cif_utils import make_label, label2cif, cif2pos1, cif2pos2
import numpy as np
from glob import glob

# result =label2cif("./debug/HALHEL_R_pacman.cif", make_label("./debug/HALHEL_R_pacman.cif", "./debug/HALHEL_A_pacman.cif"), "./debug")

# print(result)

pos1 = cif2pos1("./debug/HALHEL_R_pacman.cif")
pos2 = cif2pos2("./debug/HALHEL_R_pacman.cif")

print(pos1, pos2)