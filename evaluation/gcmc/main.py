from raspa3_tool import ff, system, cif, run


def Widom(cif_file,
          t,
          gas_files,
          paras_file,
          n_cycles,
          print_every):
    
    gases = []
    maxcutoff, forcefield = ff.paras(cif_file=cif_file, gas_files=gas_files, paras_file=paras_file)
    for i, gas_file in enumerate(gas_files):
        uc = cif.unit_cell(cif_file=cif_file, cutoff=maxcutoff)
        adsorbent = system.framework(id=0, cif_file=cif_file, uc=uc, ff=forcefield)
        mcprob = system.probs(widomProbability=1)
        gas = system.adsorbate(id=i, ff=forcefield, file=gas_file, probs=mcprob,
                               n_blocks=5, n_bins=41, f=1.0, thermodynamicIntegration=False)
        gases.append(gas)
    system_job = system.setting(id=0, ff=forcefield, temperature=t, pressure=None,
                                framework=adsorbent, components=gases,
                                n_adsorbates_init=[0] * len(gases), sysMC=mcprob)
    mc = run.setting(system=system_job, n_cycles=n_cycles, n_cycles_init=0,
                     print_every=print_every, save=True)
    run.job(mc)


def Single(cif_file,
           t,
           p,
           gas_file,
           paras_file,
           n_cycles,
           n_cycles_init,
           print_every):
    
    maxcutoff, forcefield = ff.paras(cif_file=cif_file, gas_files=[gas_file], paras_file=paras_file)
    uc = cif.unit_cell(cif_file=cif_file, cutoff=maxcutoff)
    adsorbent = system.framework(id=0, cif_file=cif_file, uc=uc, ff=forcefield)
    mcprob = system.probs(translationProbability=0.5, rotationProbability=0.5,
                          reinsertionCBMCProbability=0.5, swapProbability=0.5)
    gas = system.adsorbate(id=0, ff=forcefield, file=gas_file, probs=mcprob,
                           n_blocks=5, n_bins=41, f=1.0, thermodynamicIntegration=False)
    system_job = system.setting(id=0, ff=forcefield, temperature=t, pressure=p,
                                framework=adsorbent, components=[gas],
                                n_adsorbates_init=[0], sysMC=mcprob)
    mc = run.setting(system=system_job, n_cycles=n_cycles, n_cycles_init=n_cycles_init,
                     print_every=print_every, save=True)
    run.job(mc)


def Mixture(cif_file,
            t,
            p,
            gas_files,
            mol_frac,
            paras_file,
            n_cycles,
            n_cycles_init,
            print_every):
    
    gases = []
    maxcutoff, forcefield = ff.paras(cif_file=cif_file, gas_files=gas_files, paras_file=paras_file)
    for i, gas_file in enumerate(gas_files):
        uc = cif.unit_cell(cif_file=cif_file, cutoff=maxcutoff)
        adsorbent = system.framework(id=0, cif_file=cif_file, uc=uc, ff=forcefield)
        mcprob = system.probs(translationProbability=0.5, rotationProbability=0.5,
                              reinsertionCBMCProbability=0.5, swapProbability=0.5)
        gas = system.adsorbate(id=i, ff=forcefield, file=gas_file, probs=mcprob,
                               n_blocks=5, n_bins=41, f=mol_frac[i], thermodynamicIntegration=False)
        gases.append(gas)
    system_job = system.setting(id=0, ff=forcefield, temperature=t, pressure=p,
                                framework=adsorbent, components=gases,
                                n_adsorbates_init=[0] * len(gases), sysMC=mcprob)
    mc = run.setting(system=system_job, n_cycles=n_cycles, n_cycles_init=n_cycles_init,
                     print_every=print_every, save=True)
    run.job(mc)


RUN_TEMPLATE = '''
import sys, os
sys.path.insert(0, {script_dir!r})
os.chdir({job_folder!r})

from main import {func_name}

{func_name}(
    {kwargs}
)
'''


def make_job(job_folder,
             func_name,
             kwargs: dict,
             script_dir: str):
    
    import os
    os.makedirs(job_folder, exist_ok=True)

    kw_lines = ",\n    ".join(f"{k}={v!r}" for k, v in kwargs.items())

    content = RUN_TEMPLATE.format(
        script_dir=script_dir,
        job_folder=job_folder,
        func_name=func_name,
        kwargs=kw_lines,
    )
    run_py = os.path.join(job_folder, "run.py")
    with open(run_py, "w") as f:
        f.write(content.strip())
    return run_py


if __name__ == '__main__':
    import yaml
    import os
    import sys
    import glob
    import shutil
    from pathlib import Path
    from tqdm import tqdm

    with open(sys.argv[1], "r") as f:
        config = yaml.safe_load(f)

    script_dir = str(Path(__file__).resolve().parent)

    for cif_file in tqdm(glob.glob(config["folder"]["cifs_folder"] + "/*.cif")):
        stem       = Path(cif_file).stem
        job_folder = os.path.join(config["folder"]["jobs_folder"], stem)
        os.makedirs(job_folder, exist_ok=True)

        cif_dest = os.path.join(job_folder, Path(cif_file).name)
        shutil.copy2(cif_file, cif_dest)

        if config["job_type"] == "Widom":
            kwargs = dict(
                cif_file  = cif_dest,
                t         = config["setting"]["temperature"],
                gas_files = config["file"]["molecules_json"],
                paras_file= config["file"]["force_field_json"],
                n_cycles  = config["setting"]["n_cycles"],
                print_every=config["setting"]["print_every"],
            )
            make_job(job_folder, "Widom", kwargs, script_dir)

        elif config["job_type"] == "Single":
            kwargs = dict(
                cif_file    = cif_dest,
                t           = config["setting"]["temperature"],
                p           = config["setting"]["pressure"],
                gas_file    = config["file"]["molecules_json"],
                paras_file  = config["file"]["force_field_json"],
                n_cycles    = config["setting"]["n_cycles"],
                n_cycles_init=config["setting"]["n_cycles_init"],
                print_every = config["setting"]["print_every"],
            )
            make_job(job_folder, "Single", kwargs, script_dir)

        elif config["job_type"] == "Mixture":
            kwargs = dict(
                cif_file    = cif_dest,
                t           = config["setting"]["temperature"],
                p           = config["setting"]["pressure"],
                gas_files   = config["file"]["molecules_json"],
                mol_frac    = config["setting"]["mol_frac"],
                paras_file  = config["file"]["force_field_json"],
                n_cycles    = config["setting"]["n_cycles"],
                n_cycles_init=config["setting"]["n_cycles_init"],
                print_every = config["setting"]["print_every"],
            )
            make_job(job_folder, "Mixture", kwargs, script_dir)
