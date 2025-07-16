import os
import subprocess
import shutil

# Directories
BASE_DIR = "/home2/mmGBSA"
WT_DIR = os.path.join(BASE_DIR, "WT_run")
MT_DIR = os.path.join(BASE_DIR, "MT_run")
MDP_DIR = os.path.join(BASE_DIR, "mdp")

GROMACS = "gmx"         # Or full path to gmx executable
GMXMMPBSA = "gmx_MMPBSA"

# MM/GBSA input file content
MMPBSA_IN = """\
&general
  startframe=1, endframe=-1, interval=1,
  verbose=1,
/

&gb
  igb=5,
/
"""

def run_cmd(cmd, cwd):
    """Run a shell command in specified directory and print output."""
    print(f"[RUN] {cmd}")
    result = subprocess.run(cmd, shell=True, cwd=cwd,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(f"[FAIL] Command failed in {cwd}: {cmd}")
        exit(1)

def fix_pdb_chains(input_pdb, output_pdb, split_resid=220):
    with open(input_pdb) as inp, open(output_pdb, "w") as out:
        for line in inp:
            if line.startswith(("ATOM", "HETATM")):
                resid = int(line[22:26])
                chain = 'A' if resid <= split_resid else 'B'
                # Change chain ID
                newline = line[:21] + chain + line[22:]
                out.write(newline)
            else:
                out.write(line)

def split_pdb_by_chain(input_pdb, workdir):
    """Split a PDB file into Chain_A.pdb and Chain_B.pdb files based on chain ID."""
    chain_a_path = os.path.join(workdir, "Chain_A.pdb")
    chain_b_path = os.path.join(workdir, "Chain_B.pdb")

    with open(input_pdb, "r") as inp, \
         open(chain_a_path, "w") as chain_a, \
         open(chain_b_path, "w") as chain_b:

        for line in inp:
            if line.startswith(("ATOM", "HETATM")):
                chain_id = line[21]
                if chain_id == "A":
                    chain_a.write(line)
                elif chain_id == "B":
                    chain_b.write(line)

    print(f"[INFO] Created Chain_A.pdb and Chain_B.pdb in {workdir}")

# # This is the crucial helper function that was causing issues, now corrected within this context
def create_protein_only_topology(original_top_file, output_top_file):
    """
    Creates a new topology file containing only the protein,
    by filtering out water and ion definitions from the original topol.top,
    while retaining the [ molecules ] section for the protein.
    """
    print(f"[INFO] Creating protein-only topology: {output_top_file}")

    with open(original_top_file, 'r') as infile:
        lines = infile.readlines()

    filtered_lines = []
    in_molecules_section = False

    # Define the exact names of your protein molecules as found in topol.top
    # Based on your previous output:
    protein_molecule_names = ["Protein_chain_A", "Protein_chain_A2"]

    for line in lines:
        stripped_line = line.strip()

        if stripped_line == '[ molecules ]':
            in_molecules_section = True
            filtered_lines.append(line) # Keep the [ molecules ] header
            continue # Go to the next line

        if in_molecules_section:
            # In the molecules section, only include protein molecule lines, comments, or blank lines
            if not stripped_line or stripped_line.startswith(';'): # Keep blank lines and comments
                filtered_lines.append(line)
            else:
                parts = stripped_line.split()
                # Check if the first part of the line (molecule name) is in our list of protein names
                if len(parts) >= 2 and parts[0] in protein_molecule_names:
                    filtered_lines.append(line) # Keep the protein line
        else:
            # Before the molecules section, filter out specific includes for water/ions
            # This list might need to be expanded based on your specific force field and GROMACS version
            if any(inc in line for inc in ['#include "spc.itp"', '#include "tip3p.itp"',
                                           '#include "amber99sb-ildn.ff/spce.itp"',
                                           '#include "amber99sb-ildn.ff/ions.itp"']):
                print(f"[DEBUG] Skipping water/ion include: {line.rstrip()}")
                continue
            filtered_lines.append(line)

    with open(output_top_file, 'w') as outfile:
        outfile.writelines(filtered_lines)

    print(f"[INFO] {output_top_file} created with only protein definitions.")


def prepare_system(workdir):
    """Prepare the molecular system: fix PDB, generate files, run GROMACS steps."""
    print(f"[STAGE] Preparing system in {workdir}")
    os.makedirs(workdir, exist_ok=True)

    model_pdb = os.path.join(workdir, "model.pdb")
    fixed_pdb = os.path.join(workdir, "model_fixed.pdb")

    # Fix chain IDs in PDB
    fix_pdb_chains(model_pdb, fixed_pdb, split_resid=220)

    # Replace original model.pdb with fixed version
    shutil.copy(fixed_pdb, model_pdb)

    # Run GROMACS setup commands (unchanged until creating protein-only files)
    run_cmd(f"{GROMACS} pdb2gmx -f model_fixed.pdb -o processed.gro -water spce -ff amber99sb-ildn", cwd=workdir)
    run_cmd(f"{GROMACS} editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt cubic", cwd=workdir)
    run_cmd(f"{GROMACS} solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top", cwd=workdir)
    run_cmd(f"{GROMACS} grompp -f {MDP_DIR}/ions.mdp -c solvated.gro -p topol.top -o ions.tpr", cwd=workdir)
    run_cmd(f"echo 'SOL' | {GROMACS} genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral", cwd=workdir)
    run_cmd(f"{GROMACS} grompp -f {MDP_DIR}/minim.mdp -c solv_ions.gro -p topol.top -o em.tpr", cwd=workdir)
    run_cmd(f"{GROMACS} mdrun -deffnm em", cwd=workdir)
    run_cmd(f"{GROMACS} grompp -f {MDP_DIR}/nvt.mdp -c em.gro -p topol.top -o nvt.tpr", cwd=workdir)
    run_cmd(f"{GROMACS} mdrun -deffnm nvt", cwd=workdir)
    run_cmd(f"{GROMACS} grompp -f {MDP_DIR}/npt.mdp -c nvt.gro -p topol.top -o npt.tpr -maxwarn 1", cwd=workdir)
    run_cmd(f"{GROMACS} mdrun -deffnm npt", cwd=workdir)
    run_cmd(f"{GROMACS} grompp -f {MDP_DIR}/md.mdp -c npt.gro -p topol.top -o md.tpr", cwd=workdir)
    run_cmd(f"{GROMACS} mdrun -deffnm md", cwd=workdir)

    # Create protein-only trajectory and initial structure for MM/GBSA
    run_cmd(f"echo 'Protein' | {GROMACS} trjconv -s md.tpr -f md.xtc -o md_protein_only.xtc -n index_chain.ndx", cwd=workdir)
    run_cmd(f"echo 'Protein' | {GROMACS} trjconv -s md.tpr -f md.xtc -o protein_only.pdb -n index_chain.ndx -dump 0", cwd=workdir)

    # *** NEW STEP: Create protein-only topology file ***
    original_topol = os.path.join(workdir, "topol.top")
    protein_only_topol = os.path.join(workdir, "protein_only.top")
    # Corrected call: removed the protein_mol_name argument as it's handled internally
    create_protein_only_topology(original_topol, protein_only_topol)

    # Minimal empty.mdp for grompp (if not already existing)
    empty_mdp_path = os.path.join(MDP_DIR, "empty.mdp")
    if not os.path.exists(empty_mdp_path):
        with open(empty_mdp_path, "w") as f:
            f.write("define = -DFLEXIBLE\nintegrator = md\nnsteps = 0\n") # Minimal MDP
        print(f"[INFO] Created minimal empty.mdp for protein-only TPR generation: {empty_mdp_path}")

    # This is the line that needed correction: using protein_only.top for -p
    run_cmd(f"{GROMACS} grompp -f {empty_mdp_path} -c protein_only.pdb -p protein_only.top -o protein_only.tpr -maxwarn 1", cwd=workdir)

    # Other trjconv commands (unchanged)
    run_cmd(f"echo '1\n0' | {GROMACS} trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center", cwd=workdir)
    run_cmd(f"echo 1 | {GROMACS} trjconv -s md.tpr -f md.xtc -o last_frame.pdb -dump -1", cwd=workdir)

    # Write the MM/GBSA input file (unchanged)
    with open(os.path.join(workdir, "mmpbsa.in"), "w") as f:
        f.write(MMPBSA_IN)
    print(f"[INFO] mmpbsa.in created.")

    # Split fixed PDB by chain (unchanged)
    split_pdb_by_chain(model_pdb, workdir)




def create_chain_index_by_residues(tpr_file_path, ndx_file_path, chainA_range, chainB_range, cwd):
    cmd = [GROMACS, 'make_ndx', '-f', os.path.basename(tpr_file_path), '-o', os.path.basename(ndx_file_path)]
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=cwd)

    chain_a_cmd = f"r {chainA_range[0]}-{chainA_range[1]}"
    chain_b_cmd = f"r {chainB_range[0]}-{chainB_range[1]}"

    commands = "\n".join([
        chain_a_cmd,
        "name 17 Chain_A",  # These numbers (17, 18) are for make_ndx internal use
        chain_b_cmd,
        "name 18 Chain_B",
        "q"
    ]) + "\n"

    out, err = p.communicate(commands)
    if p.returncode != 0:
        raise RuntimeError(f"make_ndx failed in {cwd}:\n{err}")
    print("[INFO] Index file created with Chain_A and Chain_B based on residue numbers.")

def run_gmx_mmpbsa(mmpbsa_in, complex_tpr_file, xtc_file, topol_file, ndx_file): # Renamed complex_structure_file to complex_tpr_file
    """Run gmx_MMPBSA using specified inputs and chain groups."""
    cmd = [
        GMXMMPBSA,
        '-O',
        '-i', mmpbsa_in,
        '-cs', complex_tpr_file,      # This will now be protein_only.tpr
        '-ct', xtc_file,              # This is md_protein_only.xtc
        '-cp', topol_file,            # This remains topol.top (full system or protein-only if you generate one)
        '-ci', ndx_file,
        '-cg', 'Chain_A', 'Chain_B',
    ]
    print(f"[RUN] {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[FAIL] gmx_MMPBSA failed:\n{result.stderr}")
    else:
        print("[INFO] gmx_MMPBSA completed successfully.")
        print(result.stdout)




if __name__ == "__main__":
    workdir = WT_DIR

    # ... (prepare_system call and other unchanged code) ...

    prepare_system(workdir) # This should create md.tpr, md.xtc, protein_only.top, and index_chain.ndx

    # --- Files for gmx_MMPBSA ---
    # !!! CRITICAL CHANGE HERE !!!
    # Use the FULL SYSTEM TPR and XTC for gmx_MMPBSA's -cs and -ct flags
    # gmx_MMPBSA will use these to read the trajectory and create intermediate PDBs,
    # then apply your index_chain.ndx to select the protein.
    tpr_for_gmxmmpbsa_cs = os.path.join(workdir, "md.tpr") # This should be md.tpr
    xtc_for_gmxmmpbsa_ct = os.path.join(workdir, "md.xtc") # This should be md.xtc

    # Protein-only topology file for gmx_MMPBSA's -cp flag
    # This is the topology file we created that only describes the protein atoms
    protein_only_topol_for_mmpbsa = os.path.join(workdir, "protein_only.top")

    # Index file for gmx_MMPBSA's -ci and -cg flags
    # This is crucial for gmx_MMPBSA to identify your protein chains in the full system.
    # It must have been created with the FULL SYSTEM TPR as reference.
    ndx = os.path.join(workdir, "index_chain.ndx")
    mmpbsa_input = os.path.join(workdir, "mmpbsa.in")

    # Create the index file (this part is correct - uses full system TPR)
    create_chain_index_by_residues(
        tpr_file_path=tpr_for_gmxmmpbsa_cs, # Use full system TPR for make_ndx
        ndx_file_path=ndx,
        chainA_range=(1, 180), # Adjust these residue ranges to match your protein's actual chains
        chainB_range=(221, 431),
        cwd=workdir
    )

    # Call gmx_MMPBSA with the updated file types
    run_gmx_mmpbsa(mmpbsa_input,
                   tpr_for_gmxmmpbsa_cs,      # Pass full system TPR
                   xtc_for_gmxmmpbsa_ct,      # Pass full system XTC
                   protein_only_topol_for_mmpbsa, # Keep protein-only TOP
                   ndx)
