from quippy.descriptors import Descriptor
from pymatgen.analysis.structure_matcher import StructureMatcher,ElementComparator
import numpy as np
from pymatgen.core import Structure
from ase.io import read

def rmsd(cif1,cif2):

    sm  = StructureMatcher(comparator=ElementComparator(),
                           ltol=0.2,
                           stol=0.3,
                           angle_tol=5,
                           primitive_cell=True,
                           scale=True,
                           attempt_supercell=True)

    struct1 = Structure.from_file(cif1,primitive=False)
    struct2 = Structure.from_file(cif2,primitive=False)

    result = sm.get_rms_dist(struct1, struct2)

    return result

def soap(cif1,cif2):
    soap_descriptor = Descriptor("soap cutoff=3 l_max=6 n_max=12 atom_sigma=0.1")

    atoms1 = read(cif1)
    atoms2 = read(cif2)

    soap_1 = soap_descriptor.calc(atoms1)['data']
    soap_2 = soap_descriptor.calc(atoms2)['data']

    d = np.linalg.norm(soap_1 - soap_2)

    return d



import os
import pandas as pd
from tqdm import tqdm


folder1 = "./clean_n"
folder2 = "./results_n/"
n_data = pd.read_csv("./results_n/mlp_results.csv")
data=[]
for cifname in tqdm(n_data["name"][:]):
    print(cifname)
    result = rmsd(os.path.join(folder2, cifname, cifname+"_opt.cif"), os.path.join(folder1, cifname+".cif"))
    d_soap = soap(os.path.join(folder2, cifname, cifname+"_opt.cif"), os.path.join(folder1, cifname+".cif"))
    if result is not None:
        data.append([cifname, result[0], result[1], d_soap])
        print(result[0], result[1], d_soap)
    else:
        data.append([cifname, None, None, d_soap])
        print(None, None, d_soap)
pd.DataFrame(data, columns=["name", "RMSD", "MSD", "SOAP"]).to_csv("./results_n/coor_diff.csv", index=False)


folder1 = "./clean_y"
folder2 = "./results_y/"
n_data = pd.read_csv("./results_y/mlp_results.csv")
data=[]
for cifname in tqdm(n_data["name"][:]):
    print(cifname)
    try:
        result = rmsd(os.path.join(folder2, cifname, cifname+"_opt.cif"), os.path.join(folder1, cifname+".cif"))
        d_soap = soap(os.path.join(folder2, cifname, cifname+"_opt.cif"), os.path.join(folder1, cifname+".cif"))
        if result is not None:
            data.append([cifname, result[0], result[1], d_soap])
            print(result[0], result[1], d_soap)
        else:
            data.append([cifname, None, None, d_soap])
            print(None, None, d_soap)
    except:
        data.append([cifname, None, None, None])
pd.DataFrame(data, columns=["name", "RMSD", "MSD", "SOAP"]).to_csv("./results_y/coor_diff.csv", index=False)