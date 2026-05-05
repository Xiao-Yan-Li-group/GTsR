from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
from tqdm import tqdm
import os, sys, json
from glob import glob
from pathlib import Path

# atoms =read(struc)
# len(atoms)

def check(cif1: str,
        cif2: str,
        ltol: float = 0.2,
        stol: float = 0.3,
        angle_tol: float = 5.0,
        primitive: bool = True) -> bool:

    s1 = Structure.from_file(cif1)
    s2 = Structure.from_file(cif2)

    matcher = StructureMatcher(
        ltol=ltol,
        stol=stol,
        angle_tol=angle_tol,
        primitive_cell=primitive,
        scale=True,
        attempt_supercell=True,
    )
    return matcher.fit(s1, s2)


clean_folder = sys.argv[1]
ori_folder = sys.argv[2]
save_json = sys.argv[3]

results = {}
for cif_path in tqdm(glob(clean_folder+"/*.cif")):
    result = check(cif_path, os.path.join(ori_folder, Path(cif_path).stem.split("_")[:-1]+".cif"))
    results[Path(cif_path).stem] = result

with open(save_json, "w") as f:
    json.dump(results, f, indent=2)