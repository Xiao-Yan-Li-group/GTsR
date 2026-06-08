from src.cif_utils import make_label, label2cif
import numpy as np
from glob import glob


# free_labels = glob("./benchmark/dataset/label/free/*.npy")
# for label_path in free_labels:
#     label = np.load(label_path)
#     label2cif("./benchmark/dataset/clean_dataset/" + label_path.split("/")[-1].replace(".npy", ".cif"), label, "./benchmark/dataset/check/free/")

all_labels = glob("./benchmark/dataset/label/all/*.npy")
for label_path in all_labels:
    try:
        label = np.load(label_path)
        label2cif("./benchmark/dataset/clean_dataset/" + label_path.split("/")[-1].replace(".npy", ".cif"), label, "./benchmark/dataset/check/all/")
    except:
        print(f"Failed to process {label_path}")