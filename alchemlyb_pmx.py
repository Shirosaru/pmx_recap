import os
import subprocess
from pmx.model import Model
from AAtutorial import AAtutorial
from Bio.PDB import PDBParser

# === SETUP WORKING ENVIRONMENT ===
fe = AAtutorial()
fe.workPath = 'workpath'
fe.pdbfile = 'inputs/I74271_clean_clean_h.pdb'
fe.original_pdbfile = fe.pdbfile
fe.mdpPath = 'inputs/mdppath'
fe.replicas = 1
fe.ff = 'amber99sb-star-ildn-mut'

# === HYDROGENATION ===
input_pdb = 'inputs/I74271_clean_clean.pdb'
clean_pdb = fe.pdbfile
subprocess.run([
    'gmx', 'pdb2gmx',
    '-f', input_pdb,
    '-o', clean_pdb,
    '-water', 'spce',
    '-ff', fe.ff,
    '-ignh'
], input='1\n', text=True)

# === MUTATION LIST (CHAIN, RESNUM, FROM, TO) ===
mutations = [('H', 118, 'W', 'S')]

aa3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

def get_residues_biopdb(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    lookup = {}
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != ' ':
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

# === PARSE STRUCTURE ===
pdb_resnames = get_residues_biopdb(fe.pdbfile)
model = Model(fe.pdbfile, rename_atoms=False)
pdb_chain_residues = {
    chain: get_pdb_residue_sequence(fe.pdbfile, chain)
    for chain, _ in pdb_resnames
}

# === REGISTER MUTATIONS ===
fe.resids = []
fe.chains = []
fe.edges = []
fe.pdb_labels = []

for chain, pdb_resnum, from_res, to_res in mutations:
    key = (chain, pdb_resnum)

    if key not in pdb_resnames:
        raise ValueError(f"Residue {pdb_resnum} in chain {chain} not found in PDB.")

    actual_from = pdb_resnames[key]
    expected_from_3 = aa3[from_res.upper()]
    mutate_to_3 = aa3[to_res.upper()]

    if actual_from != expected_from_3:
        raise ValueError(f"Mismatch at {chain}{pdb_resnum}: expected {expected_from_3}, got {actual_from}")

    if expected_from_3 == mutate_to_3:
        raise ValueError(f"Redundant mutation {expected_from_3} → {mutate_to_3}")

    missing_before = get_missing_residue_count(pdb_chain_residues[chain], pdb_resnum)
    corrected_pmx_id = pdb_resnum - missing_before

    pmx_res = next((r for r in model.residues if r.chain.id == chain and r.id == corrected_pmx_id), None)
    if not pmx_res:
        raise ValueError(f"PMX residue not found for chain {chain}, ID {corrected_pmx_id}")

    fe.resids.append(pmx_res.id)
    fe.chains.append(chain)
    fe.edges.append([expected_from_3, mutate_to_3])
    fe.pdb_labels.append((actual_from, pdb_resnum))

    print(f"✅ Mutation registered: {expected_from_3}{pdb_resnum} → {mutate_to_3} (chain {chain})")

# === SETUP WORKDIR ===
fe.prepareFreeEnergyDir()

# === HYBRID TOPOLOGY ===
fe.hybrid_structure_topology(bVerbose=True)

# === SOLVATION + IONS ===
fe.boxWaterIons()

# === ENERGY MINIMIZATION ===
fe.prepare_simulation(simType='em')
fe.run_simulation_locally(simType='em')

# === JOB CONFIG ===
fe.JOBsimcpu = 16
fe.JOBqueue = 'SLURM'
fe.JOBsource = ['/etc/profile.d/modules.sh', '/usr/local/gromacs/GMXRC2024']
fe.JOBmodules = ['shared', 'owl/intel-mpi-default', 'cuda91']
fe.JOBbGPU = True
fe.JOBgmx = 'gmx_threads_AVX_256 mdrun -pin on'
fe.JOBpartition = 'p20,p24,p00,p32,p08,p16'

# === EM JOBS ===
fe.prepare_jobscripts(simType='em')

# === EQUILIBRATION ===
fe.prepare_simulation(simType='eq')
fe.prepare_jobscripts(simType='eq')

# === TRANSITIONS ===
fe.workPath = 'workpath_precalculated'
fe.prepare_transitions()

fe.JOBsimcpu = 10
fe.JOBntomp = 2
fe.JOBnotr = 100
fe.JOBgmx = 'srun gmx_mpi_AVX_256 mdrun -pin on'

fe.prepare_jobscripts(simType='transitions')

# === FINAL ANALYSIS ===
fe.run_analysis(bVerbose=True)
