import sys, os
sys.path.insert(0, '/mnt/d/Project/GTSR/benchmark/gcmc')
os.chdir('/mnt/d/Project/GTSR/benchmark/gcmc/data/jobs/widom/YAPCOL_R_pacman')

from main import Widom

Widom(
    cif_file='/mnt/d/Project/GTSR/benchmark/gcmc/data/jobs/widom/YAPCOL_R_pacman/YAPCOL_R_pacman.cif',
    t=298,
    gas_files=['/mnt/d/Project/GTSR/benchmark/gcmc/preparation/h2o.json'],
    paras_file='/mnt/d/Project/GTSR/benchmark/gcmc/preparation/paras.json',
    n_cycles=5000,
    print_every=1000
)