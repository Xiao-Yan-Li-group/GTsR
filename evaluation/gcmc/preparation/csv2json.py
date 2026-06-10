import pandas as pd
import json
from tqdm import tqdm
from mendeleev import element
import numpy as np
import sys
import math
import pandas as pd


component_type = sys.argv[1]
input_csv = sys.argv[2]

data = pd.read_csv(input_csv)

def default(o):
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    raise TypeError(f"{type(o)}")


paras = {}
if component_type=="ff":
    
    paras["cutOffFrameworkVDW"]=data["cutOffFrameworkVDW"][0]
    paras["cutOffMoleculeVDW"]=data["cutOffMoleculeVDW"][0]
    paras["cutOffCoulomb"]=data["cutOffCoulomb"][0]
    paras["shifted"]=data["shifted"][0]
    paras["tailCorrections"]=data["tailCorrections"][0]
    paras["useCharge"]=data["useCharge"][0]

    for i, label in tqdm(enumerate(data["label"])):
        paras[label] = {}
        paras[label]["framework"] = data["framework"][i]
        paras[label]["mass"] = data["mass"][i]
        paras[label]["charge"] = data["charge"][i]
        paras[label]["polarizability"] = data["polarizability"][i]
        if pd.isna(data["atomicNumber"][i]):
            e = element(label)
            paras[label]["atomicNumber"] = e.atomic_number
        else:
            paras[label]["atomicNumber"] = int(data["atomicNumber"][i])
        paras[label]["printToPDB"] = data["printToPDB"][i]
        paras[label]["source"] = data["source"][i]
        paras[label]["epsilon"] = data["epsilon"][i]
        paras[label]["sigma"] = data["sigma"][i]

elif component_type=="adsorbate":
    paras["CriticalTemperature"] = float(data["CriticalTemperature"][0])
    paras["CriticalPressure"] = float(data["CriticalPressure"][0])
    paras["AcentricFactor"] = float(data["AcentricFactor"][0])

    paras["Type"] = data["Type"][0]

    paras["pseudoAtoms"] = []
    for i, atom in enumerate(data["pseudoAtoms"]):
        if pd.isna(atom):
            continue
        pos = [atom]
        pos.append([float(data["x"][i]), float(data["y"][i]), float(data["z"][i])])
        paras["pseudoAtoms"].append(pos)

    if data["Type"][0] == "rigid":
        print("using rigid adsorbate")
    elif data["Type"][0] == "flexible":
        paras["Bonds"] = []
        print("using flexible adsorbate")
        for i, j in zip(data["BondsIndex1"], data["BondsIndex2"]):
        
            if math.isnan(i) or math.isnan(j):
                continue

            paras["Bonds"].append([int(i), int(j)])


with open(input_csv.replace(".csv", ".json"), "w") as f:
    json.dump(paras, f, indent=2, default=default)
