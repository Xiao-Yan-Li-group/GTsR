from pymatgen.core import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher
import os, sys, json, glob, warnings, emoji, tqdm, pathlib
from ase.io import read
import pandas as pd
from pymatgen.io.ase import AseAtomsAdaptor


warnings.filterwarnings("ignore", category=UserWarning, module="ase")
warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")


def check(cif1: str, cif2: str,
          ltol: float = 0.2, stol: float = 0.3,
          angle_tol: float = 5.0, primitive: bool = True) -> bool:
    try:
        s1, s2 = Structure.from_file(cif1), Structure.from_file(cif2)
    except:
        struc1 = read(cif1)
        s1 = AseAtomsAdaptor.get_structure(struc1)
        struc2 = read(cif2)
        s2 = AseAtomsAdaptor.get_structure(struc2)
    return StructureMatcher(
        ltol=ltol, stol=stol, angle_tol=angle_tol,
        primitive_cell=primitive, scale=False, attempt_supercell=False,
    ).fit(s1, s2)


def check_info(info, refcode, ori_type, clean_type) -> list:
    row = info[info["refcode"] == refcode].iloc[0]
    r2f, f2a = int(row["r2f"]), int(row["f2a"])

    if ori_type in ("A", "IA"):
        return ["_A", "_IA"] if ori_type == "A" else ["_IA", "_A"]

    if clean_type == "F":
        if r2f == 0:
            return ["_R"]
        return ["_IF", "_F"] if ori_type == "IF" else ["_F", "_IF"]

    if ori_type in ("F", "IF"):
        if f2a == 1:
            return ["_A", "_IA"] if ori_type == "F" else ["_IA", "_A"]
        if r2f == 1:
            return ["_F", "_IF"] if ori_type == "F" else ["_IF", "_F"]
        return ["_R"]

    if f2a == 1:
        return ["_A", "_IA"]
    if r2f == 1:
        return ["_F", "_IF"]
    return ["_R"]


def check_one(cif1, cif2) -> bool:
    try:
        return len(read(cif1)) == len(read(cif2)) and check(cif1, cif2)
    except Exception:
        return False


folder_clean   = sys.argv[1]
folder_dataset = sys.argv[2]
info           = pd.read_csv(sys.argv[3])
save_json      = sys.argv[4]

ori_types   = ["R", "F", "IF", "A", "IA"][:]
clean_types = ["A"] # "F", 
results     = {}
W = 54

print(emoji.emojize("\n:magnifying_glass_tilted_right: Structure Matching", language='alias'))
print("=" * W)

for clean_type in clean_types:
    print(f"\n+-- Clean type: {clean_type} " + "-" * (W - 17))
    results[clean_type] = {}
    s = f = 0
    fails = []
    for ori_type in ori_types:
        cifs = glob.glob(folder_dataset + f"/*_{ori_type}.cif")[:]
        print(f"|  > Ori: {ori_type:<4}  [{len(cifs)} CIFs]")
        for cif in tqdm.tqdm(cifs, desc=f"|    {ori_type}", leave=False, ncols=W + 10):
            refcode      = pathlib.Path(cif).stem
            pure_refcode = refcode.rsplit("_", 1)[0]
            candidates   = check_info(info, pure_refcode, ori_type, clean_type)

            matched = any(check_one(os.path.join(folder_dataset, pure_refcode + suff + ".cif"),
                        os.path.join(folder_clean, refcode + "_samosa_" + clean_type + ".cif")) for suff in candidates)
                          # _core_, #_csd_, _mofid_StandardIsolated_, _molsimplify_1_
            
            if matched:
                s += 1
            else:
                tqdm.tqdm.write(emoji.emojize(f'  :x: FAIL  {refcode}', language='alias'))

                f += 1
                fails.append(refcode)
    total = s + f
    ratio = s / total if total else 0
    print("|")
    print("|  Summary")
    print(emoji.emojize(f"|    :check_mark_button: success : {s:>4}  ({ratio:.1%})", language='alias'))
    print(emoji.emojize(f"|    :cross_mark: fail    : {f:>4}", language='alias'))
    print("+" + "-" * (W - 1))
    results[clean_type] = {
        "fails":   fails,
        "summary": {"success": s, "fail": f, "ratio": ratio},
    }

with open(save_json, "w") as f:
    json.dump(results, f, indent=2)