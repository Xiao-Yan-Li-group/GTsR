import os
import glob
import argparse
import time
import traceback

from ccdc import io
from tqdm import tqdm

arg_handler = argparse.ArgumentParser(description=__doc__)
arg_handler.add_argument(
    'input',
    help='CIF file or directory of CIF files'
)
arg_handler.add_argument(
    '-o', '--output-directory', required=True,
    help='Directory into which to write stripped structures'
)
arg_handler.add_argument(
    '-s', '--solvent-file',
    help='Location of solvent file (directory of .mol2 or single file)'
)
arg_handler.add_argument(
    '-v', '--verbose', default=False, action='store_true',
    help='Print component SMILES for each structure (for debugging)'
)

args = arg_handler.parse_args()

if os.path.isdir(args.input):
    input_files = sorted(glob.glob(os.path.join(args.input, '*.cif')))
    if not input_files:
        raise SystemExit(f'No CIF files found in: {args.input}')
else:
    input_files = [args.input]

os.makedirs(args.output_directory, exist_ok=True)

if not args.solvent_file:
    args.solvent_file = os.path.normpath(os.path.join(
        os.path.dirname(io.__file__), 'resources', 'mercury', 'ccdc_solvents'
    ))

if os.path.isdir(args.solvent_file):
    solvent_smiles = {
        io.MoleculeReader(f)[0].smiles
        for f in glob.glob(os.path.join(args.solvent_file, '*.mol2'))
    }
else:
    solvent_smiles = {m.smiles for m in io.MoleculeReader(args.solvent_file)}


_WATER_SMILES = {'O', '[OH2]', '[O]'}  # CSD may emit any of these for water/isolated O

def is_solvent(c):
    return c.smiles in _WATER_SMILES or (c.smiles is not None and c.smiles in solvent_smiles)


def has_metal(c):
    return any(a.is_metal for a in c.atoms)


def is_free(c, mol):
    for a in c.atoms:
        orig_a = mol.atom(a.label)
        if any(x.is_metal for b in orig_a.bonds for x in b.atoms):
            return False
    return True


def strip_mol(mol, free_only):
    mol_work = mol.copy()
    clone = mol_work.copy()
    clone.remove_bonds(b for b in clone.bonds if any(a.is_metal for a in b.atoms))
    to_remove = [
        c for c in clone.components
        if not has_metal(c) and is_solvent(c) and (not free_only or is_free(c, mol_work))
    ]
    
    mol_work.remove_atoms(
        mol_work.atom(a.label) for c in to_remove for a in c.atoms
    )
    return mol_work


failed = []
processed = 0
start_time = time.time()

for cif_path in tqdm(input_files, desc='Processing', unit='file'):
    stem = os.path.splitext(os.path.basename(cif_path))[0]
    t0 = time.time()
    try:
        for entry in io.EntryReader(cif_path):
            if not entry.has_3d_structure:
                continue
            identifier = stem
            mol_orig = entry.molecule
            mol_orig.normalise_labels()
            mol_orig.assign_bond_types()

            if args.verbose:
                clone_dbg = mol_orig.copy()
                clone_dbg.remove_bonds(
                    b for b in clone_dbg.bonds if any(a.is_metal for a in b.atoms)
                )
                tqdm.write(f'\n[{stem}] components after removing metal bonds:')
                for c in clone_dbg.components:
                    flag = 'SOLVENT' if is_solvent(c) else 'keep'
                    metal = ' (has metal)' if has_metal(c) else ''
                    tqdm.write(f'  [{flag}]{metal} smiles={c.smiles!r}  natoms={len(c.atoms)}')
                tqdm.write(f'  solvent library size: {len(solvent_smiles)}')

            for suffix, free_only in [('A', False), ('F', True)]:
                mol_out = strip_mol(mol_orig, free_only=free_only)
                entry.crystal.molecule = mol_out
                out_path = os.path.join(
                    args.output_directory, f'{identifier}_csd_{suffix}.cif'
                )
                with io.CrystalWriter(out_path) as writer:
                    writer.write(entry.crystal)
            processed += 1
    except Exception:
        elapsed = time.time() - t0
        failed.append((stem, elapsed, traceback.format_exc()))

total_elapsed = time.time() - start_time

failed_log = os.path.join(args.output_directory, 'failed_list.txt')
with open(failed_log, 'w', encoding='utf-8') as f:
    f.write(f'Total time : {total_elapsed:.1f}s\n')
    f.write(f'Processed  : {processed}\n')
    f.write(f'Failed     : {len(failed)}\n')
    if failed:
        f.write('\n--- Failed entries ---\n\n')
        for name, t, tb in failed:
            f.write(f'{name}  ({t:.2f}s)\n{tb}\n')

print(f'\nDone: {processed} processed, {len(failed)} failed — {total_elapsed:.1f}s')
print(f'Log : {failed_log}')
