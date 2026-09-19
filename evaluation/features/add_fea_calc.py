from ase.io import read
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor

import sys
import pandas as pd
from tqdm import tqdm

data = pd.read_csv("./zeopp.csv")

def SpaceGroup(structure):   
    atoms = read(structure)
    structure_ = AseAtomsAdaptor.get_structure(atoms)
    result_ = SpacegroupAnalyzer(structure_, symprec=0.01, angle_tolerance=5)
    space_group_number = result_.get_space_group_number()
    return space_group_number

def n_atom(structure):
    atoms = read(structure)
    number_atoms = len(atoms)
    return number_atoms

def Mass(structure):
    atoms = read(structure)
    total_mass = atoms.get_masses().sum()
    return total_mass

path = "../../evaluation/dataset/clean_dataset/"
data_add = []
for name in tqdm(data["refcode"][:], file=sys.stderr):
    try:
        sg = SpaceGroup(path+name+".cif")
        na = n_atom(path+name+".cif")
        ave_mass = Mass(path+name+".cif")/na
        data_add.append([name, sg, na, ave_mass])
        print(name, sg, na, ave_mass)
    except:
        print(name, "fail", "fail", "fail")

pd.DataFrame(data_add, columns=["refcode", "sg", "n_atom", "ave_mass"]).to_csv("add_fea.csv", index=False)