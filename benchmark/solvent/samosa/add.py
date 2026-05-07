import os,shutil, glob, pathlib, tqdm


fail_list = ["7214458_R.cif", "HOBGAH_R.cif", "ja512973b_si_002_R.cif", "YAPCOL_R.cif "]
for cif in tqdm.tqdm(glob.glob("../../../dataset/clean_dataset/*.cif")[:]):
    clean_cif = os.path.join("../results/", pathlib.Path(cif).stem+"_samosa_A.cif")
    if os.path.exists(clean_cif):
        pass
    else:
        if clean_cif in fail_list:
            pass
        else:
            shutil.copy(cif, clean_cif)
            print(clean_cif)

for cif in tqdm.tqdm(glob.glob("../../../dataset/clean_dataset/*.cif")[:]):
    clean_cif = os.path.join("../results/", pathlib.Path(cif).stem+"_samosa_F.cif")
    if os.path.exists(clean_cif):
        pass
    else:
        if clean_cif in fail_list:
            pass
        else:
            shutil.copy(cif, clean_cif)
            print(clean_cif)