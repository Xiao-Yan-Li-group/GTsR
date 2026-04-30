import json, raspalib
from raspa3_tool.cif import get_label
    
def paras(cif_file, gas_files, paras_file):

    ff_paras = json.load(open(paras_file))
    framework_labels= get_label(cif_file=cif_file)
    
    atomTypes, parameters = [], []
    for label in list(dict.fromkeys(framework_labels)):
        info = ff_paras[label]
        atomTypes.append(
            raspalib.PseudoAtom(
                name=label,
                frameworkType=info["framework"],
                mass=info["mass"],
                charge=info["charge"],
                polarizability=info["polarizability"],
                atomicNumber=info["atomicNumber"],
                printToPDB=info["printToPDB"],
                source=info["source"]
                )
            )
        parameters.append(raspalib.VDWParameters(info["epsilon"], info["sigma"]))

    for gas_file in gas_files:
        adsorbate_labels = json.load(open(gas_file))["pseudoAtoms"]
        
        for label in adsorbate_labels:
            info = ff_paras[label[0]]
            atomTypes.append(
                raspalib.PseudoAtom(
                    name=label[0],
                    frameworkType=info["framework"],
                    mass=info["mass"],
                    charge=info["charge"],
                    polarizability=info["polarizability"],
                    atomicNumber=info["atomicNumber"],
                    printToPDB=info["printToPDB"],
                    source=info["source"]
                    )
                )
            parameters.append(raspalib.VDWParameters(info["epsilon"], info["sigma"]))

    maxcutoff = max(ff_paras["cutOffFrameworkVDW"],
                    ff_paras["cutOffMoleculeVDW"],
                    ff_paras["cutOffCoulomb"])

    return  maxcutoff, raspalib.ForceField(
                        pseudoAtoms=atomTypes,
                        parameters=parameters,
                        mixingRule=raspalib.ForceField.MixingRule.Lorentz_Berthelot,
                        cutOffFrameworkVDW=ff_paras["cutOffFrameworkVDW"],
                        cutOffMoleculeVDW=ff_paras["cutOffMoleculeVDW"],
                        cutOffCoulomb=ff_paras["cutOffCoulomb"],
                        shifted=ff_paras["shifted"],
                        tailCorrections=ff_paras["tailCorrections"],
                        useCharge=ff_paras["useCharge"]
                        )
