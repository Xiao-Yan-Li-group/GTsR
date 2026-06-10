from PACMANCharge import pmcharge
import sys
import warnings
warnings.filterwarnings(
    "ignore",
    message="pkg_resources is deprecated",
    category=UserWarning,
)


def run_pacman(cif_file,
                charge_type="DDEC6",
                digits=6,
                atom_type=True,
                neutral=True,
                keep_connect=False):
    pmcharge.predict(
                    cif_file=cif_file,
                    charge_type=charge_type,
                    digits=digits,
                    atom_type=atom_type,
                    neutral=neutral,
                    keep_connect=keep_connect
                    )
    
run_pacman(sys.argv[1])