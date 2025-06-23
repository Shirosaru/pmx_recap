import pmx
from pmx.utils import create_folder
from pmx import gmx, jobscript
import sys
import numpy as np
from IPython.core.display import clear_output
import os,shutil
import re
import subprocess
import glob
import random
from AAtutorial import *

from pmx.model import Model


# === INITIALIZE FREE ENERGY ENVIRONMENT ===
fe = AAtutorial()

# === CONFIGURATION ===
fe.workPath = 'workpath'
# Input PDB (initial, possibly missing atoms)
input_pdb = 'inputs/I74271_clean_clean.pdb'

# Output cleaned PDB (with hydrogens and complete residues)
clean_pdb = 'inputs/I74271_clean_clean_h.pdb'

# Run GROMACS pdb2gmx to generate a complete structure
subprocess.run([
    'gmx', 'pdb2gmx',
    '-f', input_pdb,
    '-o', clean_pdb,
    '-water', 'spce',
    '-ff', 'amber99sb-star-ildn-mut',
    '-ignh'  # <-- Ignore existing hydrogens
], input='1\n', text=True)

# Use the cleaned PDB in PMX
fe.pdbfile = clean_pdb
fe.original_pdbfile = clean_pdb
fe.original_pdbfile = fe.pdbfile
fe.mdpPath = 'inputs/mdppath'
fe.replicas = 1
fe.ff = 'amber99sb-star-ildn-mut'

# === USER MUTATIONS: (chain, pdb_resnum, from_res, to_res) ===
mutations = [
    ('H', 118, 'W', 'S'),  # TRP118 → SER
]

# === Amino Acid Translator ===
aa3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

# === STEP 1: Parse original PDB using Bio.PDB ===
def get_residues_biopdb(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    lookup = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != ' ':  # Skip HETATM, water, etc.
                    continue
                chain_id = chain.id
                resnum = res.id[1]
                resname = res.get_resname().upper()
                lookup[(chain_id, resnum)] = resname
    return lookup

def get_pdb_residue_sequence(pdb_path, target_chain):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("bio", pdb_path)
    pdb_indices = []
    for model in structure:
        for chain in model:
            if chain.id != target_chain:
                continue
            for res in chain:
                if res.id[0] == " ":
                    pdb_indices.append(res.id[1])
    return sorted(set(pdb_indices))

def get_missing_residue_count(pdb_indices, target_resnum):
    return sum(1 for r in range(1, target_resnum) if r not in pdb_indices)

# === Parse PDB ===
pdb_resnames = get_residues_biopdb(fe.pdbfile)

# === Load PMX Model ===
model = Model(fe.pdbfile, rename_atoms=False)

# === Get chain-wise residue sequence from Bio.PDB ===
pdb_chain_residues = {}
for chain_id in set(c for (c, _) in pdb_resnames.keys()):
    pdb_chain_residues[chain_id] = get_pdb_residue_sequence(fe.pdbfile, chain_id)

# === STEP 3: Validate and Register Mutations ===
fe.resids = []
fe.chains = []
fe.edges = []
fe.pdb_labels = []

for chain, pdb_resnum, from_res, to_res in mutations:
    key = (chain, pdb_resnum)

    if key not in pdb_resnames:
        raise ValueError(f"❌ Residue {pdb_resnum} in chain {chain} not found in Bio.PDB-parsed PDB.")

    actual_from = pdb_resnames[key]
    expected_from_3 = aa3.get(from_res.upper(), from_res.upper())
    mutate_to_3 = aa3.get(to_res.upper(), to_res.upper())

    if actual_from != expected_from_3:
        raise ValueError(f"❌ Residue mismatch at {chain}{pdb_resnum}: expected {expected_from_3}, got {actual_from}")

    if expected_from_3 == mutate_to_3:
        raise ValueError(f"❌ Redundant mutation at {chain}{pdb_resnum}: {expected_from_3} → {mutate_to_3}")

    # === Compute corrected PMX index ===
    chain_resnums = pdb_chain_residues[chain]
    missing_before = get_missing_residue_count(chain_resnums, pdb_resnum)
    corrected_pmx_id = pdb_resnum - missing_before

    # === Lookup corresponding PMX residue ===
    pmx_res = next(
        (r for r in model.residues if r.chain.id == chain and r.id == corrected_pmx_id),
        None
    )
    if not pmx_res:
        raise ValueError(f"❌ PMX residue not found for chain {chain}, corrected PMX ID {corrected_pmx_id}")

    # === Register mutation ===
    fe.resids.append(pmx_res.id)
    fe.chains.append(chain)
    fe.edges.append([expected_from_3, mutate_to_3])
    fe.pdb_labels.append((actual_from, pdb_resnum))

    print(f"✅ Registered mutation: {expected_from_3}{pdb_resnum} (chain {chain}, PMX ID {pmx_res.id}) → {mutate_to_3}")

# === STEP 4: Prepare free energy directories ===
fe.prepareFreeEnergyDir()

# === STEP 5: Build hybrid topology ===
fe.hybrid_structure_topology(bVerbose=True)

# === STEP 6: Solvate system ===
fe.boxWaterIons()

# === STEP 7: Prepare & Run EM locally ===
fe.prepare_simulation(simType='em')
fe.run_simulation_locally(simType='em')

# === STEP 8: Setup SLURM job configs ===
fe.JOBsimcpu = 16
fe.JOBqueue = 'SLURM'
fe.JOBsource = ['/etc/profile.d/modules.sh', '/usr/local/gromacs/GMXRC2024']
fe.JOBmodules = ['shared', 'owl/intel-mpi-default', 'cuda91']
fe.JOBbGPU = True
fe.JOBgmx = 'gmx_threads_AVX_256 mdrun -pin on'
fe.JOBpartition = 'p20,p24,p00,p32,p08,p16'

# === STEP 9: Create EM jobscripts ===
fe.prepare_jobscripts(simType='em')

# === STEP 10: Prepare & Submit Equilibration ===
fe.prepare_simulation(simType='eq')
fe.prepare_jobscripts(simType='eq')

# === STEP 11: Prepare transitions ===
fe.workPath = 'workpath_precalculated'
fe.prepare_transitions()

fe.JOBsimcpu = 10
fe.JOBntomp = 2
fe.JOBnotr = 100
fe.JOBqueue = 'SLURM'
fe.JOBsource = ['/etc/profile.d/modules.sh', '/usr/local/gromacs/GMXRC2024']
fe.JOBmodules = ['shared', 'owl/intel-mpi-default', 'cuda91']
fe.JOBbGPU = True
fe.JOBgmx = 'srun gmx_mpi_AVX_256 mdrun -pin on'
fe.JOBpartition = 'p20,p24,p00,p32,p08,p16'

# === STEP 12: Final jobscripts & Analysis ===
fe.prepare_jobscripts(simType='transitions')
fe.run_analysis(bVerbose=True)