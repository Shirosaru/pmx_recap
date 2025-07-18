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
  startframe=1,
  interval=1,
  verbose=1,
  keep_files=1,
/

&gb
  igb=5,
  saltcon=0.150,
/

&decomp
  idecomp=1,
  dec_verbose=1,
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
    """Prepare the molecular system for MM/GBSA: clean water/ions, generate protein-only files."""
    print(f"[STAGE] Preparing system in {workdir}")
    os.makedirs(workdir, exist_ok=True)

    model_pdb = os.path.join(workdir, "model.pdb")
    fixed_pdb = os.path.join(workdir, "model_fixed.pdb")

    # Fix chain IDs in PDB based on residue number
    fix_pdb_chains(model_pdb, fixed_pdb, split_resid=220)
    shutil.copy(fixed_pdb, model_pdb)

    # === Full system setup ===
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

    # === Create a full-system index file for later use ===
    # This index file should contain a "Protein" group
    full_system_ndx = os.path.join(workdir, "full_system.ndx")
    cmd_make_ndx_full = [GROMACS, 'make_ndx', '-f', 'md.tpr', '-o', os.path.basename(full_system_ndx)]
    p_full = subprocess.Popen(cmd_make_ndx_full, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=workdir)
    commands_full = "a Protein\nq\n" # Select all protein atoms by default in make_ndx
    out_full, err_full = p_full.communicate(commands_full)
    if p_full.returncode != 0:
        raise RuntimeError(f"make_ndx failed for full system in {workdir}:\n{err_full}")
    print("[INFO] Full system index file created with 'Protein' group.")


    # === Extract initial protein structure from full system trajectory ===
    # Use the full system TPR and the newly created full_system.ndx
    # to extract the protein from the *full system trajectory (md.xtc)*.
    # Group 'Protein' is usually 1 in the default index file after pdb2gmx.
    run_cmd(f"echo 'Protein' | {GROMACS} trjconv -s md.tpr -f md.xtc -o protein_only.pdb -dump 0 -n {os.path.basename(full_system_ndx)}", cwd=workdir)


    # Sanity check: make sure no solvent or ions
    with open(os.path.join(workdir, "protein_only.pdb")) as f:
        pdb_content = f.read()
        if any(res in pdb_content for res in ("SOL", "NA", "CL")):
            raise RuntimeError("protein_only.pdb still contains water or ions!")

    # === Create protein-only topology ===
    original_topol = os.path.join(workdir, "topol.top")
    protein_topol = os.path.join(workdir, "protein_only.top")
    create_protein_only_topology(original_topol, protein_topol)

    # === Create protein-only TPR (This is still correct for gmx_MMPBSA's -cs) ===
    # This TPR describes ONLY the protein.
    empty_mdp = os.path.join(MDP_DIR, "empty.mdp")
    if not os.path.exists(empty_mdp):
        with open(empty_mdp, "w") as f:
            f.write("define = -DFLEXIBLE\nintegrator = md\nnsteps = 0\n")
    run_cmd(f"{GROMACS} grompp -f {empty_mdp} -c protein_only.pdb -p protein_only.top -o protein_only.tpr -maxwarn 1", cwd=workdir)


    # === Extract protein-only trajectory from the full trajectory using full system TPR/NDX ===
    # This is the corrected step for protein_only.xtc
    # We use 'Protein' from the 'full_system.ndx' and the 'md.tpr' to correctly select
    # protein atoms from the 'md.xtc' and write them to 'protein_only.xtc'.
    run_cmd(f"echo 'Protein' | {GROMACS} trjconv -s md.tpr -f md.xtc -n {os.path.basename(full_system_ndx)} -o protein_only.xtc", cwd=workdir)


    # === Optional outputs (last frame, centered noPBC) ===
    # Select 'Protein' for centering AND 'Protein' for output.
    run_cmd(f"echo 'Protein\nProtein' | {GROMACS} trjconv -s md.tpr -f md.xtc -o md_noPBC_protein.xtc -pbc mol -center -n {os.path.basename(full_system_ndx)}", cwd=workdir)

    # Select 'Protein' for output for the last frame pdb.
    run_cmd(f"echo 'Protein' | {GROMACS} trjconv -s md.tpr -f md.xtc -o last_frame_protein.pdb -dump -1 -n {os.path.basename(full_system_ndx)}", cwd=workdir)


    # === Write MM/GBSA input file ===
    with open(os.path.join(workdir, "mmpbsa.in"), "w") as f:
        f.write(MMPBSA_IN)

    # === Split initial PDB by chains ===
    split_pdb_by_chain(model_pdb, workdir)

    print(f"[INFO] System prepared for MM/GBSA in {workdir}")


def create_chain_index_by_residues(tpr_file_path, ndx_file_path, chainA_range, chainB_range, cwd):
    print(f"[INFO] Creating index file for chains: {ndx_file_path}")

    # Use a direct input string for gmx make_ndx
    # We now know (from previous debugging) that these will likely be groups 10 and 11.
    # This assumes no other custom groups are added before these.
    # If your default group list changes, these numbers might need adjustment.
    commands = f"""
r {chainA_range[0]}-{chainA_range[1]}
r {chainB_range[0]}-{chainB_range[1]}
name 10 Chain_A  # Assign name to the group you expect to be number 10
name 11 Chain_B  # Assign name to the group you expect to be number 11
q
"""
    # Verify these numbers (10 and 11) by running `gmx make_ndx -f protein_only.tpr`
    # interactively, adding `r 1-180` and `r 181-431`, then typing `ls` to see assigned group numbers.
    # Based on your previous output, default groups go up to 9. So 10 and 11 are next.


    print(f"[DEBUG] Commands for make_ndx:\n{commands}")

    # The -f flag takes the structure file (tpr_file_path), not an existing ndx.
    # make_ndx will create the index file from scratch or add to it based on interactive commands.
    cmd = [GROMACS, 'make_ndx', '-f', os.path.basename(tpr_file_path), '-o', os.path.basename(ndx_file_path)]
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=cwd)

    out, err = p.communicate(commands)
    if p.returncode != 0:
        print(f"[ERROR] make_ndx failed in {cwd}:")
        print(f"STDOUT:\n{out}")
        print(f"STDERR:\n{err}")
        raise RuntimeError(f"make_ndx failed for {ndx_file_path}")
    print("[INFO] Index file created with Chain_A and Chain_B based on residue numbers.")



def create_chain_index_by_residues_robust(tpr_file_path, ndx_file_path, chainA_range, chainB_range, cwd):
    print(f"[INFO] Creating index file for chains (robust method): {ndx_file_path}")

    temp_ndx_file = os.path.join(cwd, "temp_make_ndx.ndx")
    cmd_init = [GROMACS, 'make_ndx', '-f', os.path.basename(tpr_file_path), '-o', os.path.basename(temp_ndx_file)]
    p_init = subprocess.Popen(cmd_init, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=cwd)

    initial_commands = f"""
r {chainA_range[0]}-{chainA_range[1]}
r {chainB_range[0]}-{chainB_range[1]}
q
"""
    out_init, err_init = p_init.communicate(initial_commands)

    if p_init.returncode != 0:
        print(f"[ERROR] Initial make_ndx failed in {cwd}:")
        print(f"STDOUT:\n{out_init}")
        print(f"STDERR:\n{err_init}")
        raise RuntimeError(f"Initial make_ndx failed for {ndx_file_path}")

    # --- THIS IS THE CRUCIAL CHANGE ---
    # Based on the stdout, the last default group was 9. So the next two are 10 and 11.
    group_numbers = [10, 11]
    # --- END CRUCIAL CHANGE ---

    # Step 2: Rename the groups using a new make_ndx call on the temp file
    cmd_rename = [GROMACS, 'make_ndx', '-f', os.path.basename(temp_ndx_file), '-o', os.path.basename(ndx_file_path)]
    p_rename = subprocess.Popen(cmd_rename, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=cwd)

    rename_commands = f"""
name {group_numbers[0]} Chain_A
name {group_numbers[1]} Chain_B
q
"""
    print(f"[DEBUG] Rename commands for make_ndx:\n{rename_commands}")

    out_rename, err_rename = p_rename.communicate(rename_commands)
    if p_rename.returncode != 0:
        print(f"[ERROR] Renaming make_ndx groups failed in {cwd}:")
        print(f"STDOUT:\n{out_rename}")
        print(f"STDERR:\n{err_rename}")
        raise RuntimeError(f"Renaming groups in make_ndx failed for {ndx_file_path}")

    os.remove(temp_ndx_file) # Clean up temporary file
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

    prepare_system(workdir) # This should create md.tpr, md.xtc, protein_only.top, protein_only.tpr, protein_only.xtc, and full_system.ndx

    # --- Files for gmx_MMPBSA ---
    # !!! CRITICAL CHANGE HERE !!!
    # gmx_MMPBSA requires the -cs and -ct to contain ONLY the atoms that will be used for MM/GBSA.
    # Therefore, we use the protein-only TPR and XTC.
    tpr_for_gmxmmpbsa_cs = os.path.join(workdir, "protein_only.tpr") # Changed to protein_only.tpr
    xtc_for_gmxmmpbsa_ct = os.path.join(workdir, "protein_only.xtc") # Changed to protein_only.xtc

    # Protein-only topology file for gmx_MMPBSA's -cp flag (already correct)
    protein_only_topol_for_mmpbsa = os.path.join(workdir, "protein_only.top")

    # Index file for gmx_MMPBSA's -ci and -cg flags
    # This index file MUST be created with the protein_only.tpr as reference
    # so its atom numbers match those in protein_only.xtc and protein_only.top.
    ndx = os.path.join(workdir, "index_chain.ndx")
    mmpbsa_input = os.path.join(workdir, "mmpbsa.in")

    # Create the index file using the PROTEIN-ONLY TPR
    # This is crucial for the index file atom numbers to align with protein_only.tpr/xtc
    create_chain_index_by_residues(
        tpr_file_path=tpr_for_gmxmmpbsa_cs, # Use protein_only.tpr for make_ndx
        ndx_file_path=ndx,
        chainA_range=(1, 180), # These ranges are RESIDUE numbers, which are absolute in the protein.
        chainB_range=(181, 431), # **Ensure this is (181, 431) if chain A has 180 residues and chain B starts after it.**
                                 # (Previously you had 221, 431, assuming 220 was split point. Re-check your actual residue numbering)
        cwd=workdir
    )

    # Call gmx_MMPBSA with the corrected protein-only files
    run_gmx_mmpbsa(mmpbsa_input,
                   tpr_for_gmxmmpbsa_cs,      # Pass protein_only TPR
                   xtc_for_gmxmmpbsa_ct,      # Pass protein_only XTC
                   protein_only_topol_for_mmpbsa, # Keep protein-only TOP
                   ndx)
