# from src.cif_utils import label2cif
# import numpy as np
# label = np.load("./webapp/logs/HALHEL_R.npy")
# label2cif("./webapp/logs/HALHEL_R_pacman.cif", label, "./")


# from src.cif_utils import get_sol_smi
# import numpy as np
# smis = get_sol_info("./webapp/logs/predictions/HALHEL_R_pacman_solvent.cif")
# print(smis)


from gtsr import GTsRunner
runner = GTsRunner(checkpoint="all")
result = runner.predict(
    cif="./webapp/logs/HALHEL_R_pacman.cif",
    output="./webapp/logs/",
    threshold=0.5,
)
solvent_smiles = result["solvent_smiles"]
# probabilities = result["probabilities"]
print(solvent_smiles)
