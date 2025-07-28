import os
import subprocess
import shutil
from collections import defaultdict
import re # Added for more robust group number parsing if needed, but simplified below

# --- Configuration (same as before) ---
BASE_DIR = "/home2/mmGBSA/test2"
WT_SOURCE_DIR = os.path.join(BASE_DIR, "WT")
MT_SOURCE_DIR = os.path.join(BASE_DIR, "Mut")

WT_RUN_DIR = os.path.join(BASE_DIR, "WT_run_auto")
MT_RUN_DIR = os.path.join(BASE_DIR, "MT_run_auto")

MDP_DIR = os.path.join(BASE_DIR, "mdp")

GROMACS = "gmx"
GMXMMPBSA = "gmx_MMPBSA"

# MM/GBSA input file content (same as before)
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

# Helper function to run commands (same as before)
def run_cmd(cmd, cwd):
    print(f"[RUN] {cmd}")
    result = subprocess.run(cmd, shell=True, cwd=cwd,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(f"[FAIL] Command failed in {cwd}: {cmd}")
        raise RuntimeError(f"Command failed with exit code {result.returncode}")
    print(f"[INFO] Command finished successfully in {cwd}: {cmd}")

# relabel_and_renumber_pdb_chains function (same as before)
def relabel_and_renumber_pdb_chains(input_pdb, output_pdb):
    """
    Reads an input PDB, reassigns chain IDs 'A' and 'B', and renumbers
    residues of Chain B to follow sequentially after Chain A.
    Inserts a TER card between chains to help GROMACS recognize separate molecules.

    This function attempts to detect the split point between chains based on
    a sudden jump in residue number or a new chain identifier starting from 1.
    It returns the determined lengths of Chain A and Chain B.

    Args:
        input_pdb (str): Path to the original PDB file.
        output_pdb (str): Path for the processed output PDB file (e.g., model.pdb).

    Returns:
        tuple: (chain_a_length, chain_b_length) in terms of number of residues.
    """
    print(f"[INFO] Relabeling and renumbering chains in {input_pdb} -> {output_pdb}")

    output_lines = []
    
    chain_a_max_resid = 0
    
    all_resids_A = set()
    all_resids_B = set()

    is_processing_chain_a = True
    chain_a_assigned_id = 'A'
    chain_b_assigned_id = 'B'

    last_original_chain_id = ''
    last_original_resid = 0

    with open(input_pdb, 'r') as infile:
        for line in infile:
            if line.startswith(("ATOM", "HETATM")):
                original_chain_id = line[21].strip()
                original_resid = int(line[22:26])

                if is_processing_chain_a:
                    # Heuristic for chain break:
                    # If we've processed some residues for chain A, AND
                    # current residue is 1 (or significantly lower than previous) AND
                    # the original chain ID has changed (e.g., from A to B) OR
                    # the residue number has reset (1) while previous was high.
                    if (len(all_resids_A) > 0 and original_resid == 1 and
                        (original_chain_id != last_original_chain_id or original_resid < last_original_resid)):
                        
                        output_lines.append("TER\n") # Insert TER before the new chain
                        print(f"[DEBUG] Detected chain transition. Appending TER. Old Chain: {last_original_chain_id} Resid: {last_original_resid}. New Chain: {original_chain_id} Resid: {original_resid}")
                        is_processing_chain_a = False
                        
                        # Calculate new_resid for the first residue of chain B
                        new_resid = (chain_a_max_resid + 1) # Chain B starts after Chain A's last original residue
                        all_resids_B.add(new_resid)
                        newline = line[:21] + chain_b_assigned_id + f"{new_resid:4d}" + line[26:]
                        output_lines.append(newline)

                    else: # Still processing Chain A
                        all_resids_A.add(original_resid)
                        chain_a_max_resid = max(chain_a_max_resid, original_resid) # Max original resid for chain A
                        newline = line[:21] + chain_a_assigned_id + line[22:]
                        output_lines.append(newline)

                else: # Processing Chain B
                    # Renumber Chain B residues sequentially after Chain A's last original residue
                    new_resid = chain_a_max_resid + (original_resid - 1) + 1 # (original_resid - 1) makes it 0-indexed offset
                    all_resids_B.add(new_resid)
                    newline = line[:21] + chain_b_assigned_id + f"{new_resid:4d}" + line[26:]
                    output_lines.append(newline)
                
                last_original_chain_id = original_chain_id
                last_original_resid = original_resid
            else:
                output_lines.append(line) # Keep non-ATOM/HETATM lines

    # Add a final TER card at the end of the file if needed
    if output_lines and output_lines[-1].strip().startswith(("ATOM", "HETATM")):
        output_lines.append("TER\n")
        print("[DEBUG] Appended final TER card.")

    with open(output_pdb, "w") as outfile:
        outfile.writelines(output_lines)
    
    actual_chain_a_length = len(all_resids_A)
    actual_chain_b_length = len(all_resids_B)

    # Fallback/sanity check if only one chain was detected
    if actual_chain_a_length > 0 and actual_chain_b_length == 0 and not is_processing_chain_a:
        # If we switched to B but B ended up empty, means it was likely a single chain.
        # Recalculate based on 'A' and if 'B' was truly not found, set its length to 0.
        print("[WARNING] Chain B appeared to be intended but was empty. Re-evaluating chain lengths.")
        all_resids_in_output = set()
        for line in output_lines:
            if line.startswith(("ATOM", "HETATM")):
                current_chain = line[21].strip()
                current_resid = int(line[22:26])
                if current_chain == chain_a_assigned_id:
                    all_resids_A.add(current_resid)
                elif current_chain == chain_b_assigned_id:
                    all_resids_B.add(current_resid)
        actual_chain_a_length = len(all_resids_A)
        actual_chain_b_length = len(all_resids_B) # This might still be 0, which is correct
        
    print(f"[INFO] Chains relabeled and renumbered to {output_pdb}. Chain A length: {actual_chain_a_length}, Chain B length: {actual_chain_b_length}.")
    
    return actual_chain_a_length, actual_chain_b_length


# create_protein_only_topology function (same as before)
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
    
    for line in lines:
        stripped_line = line.strip()

        if stripped_line == '[ molecules ]':
            in_molecules_section = True
            filtered_lines.append(line)
            continue

        if in_molecules_section:
            if not stripped_line or stripped_line.startswith(';'):
                filtered_lines.append(line)
            else:
                parts = stripped_line.split()
                if len(parts) >= 2:
                    if not any(mol_name in parts[0].upper() for mol_name in ["SOL", "HOH", "WAT", "NA", "CL", "ION"]):
                        filtered_lines.append(line)
        else:
            if any(inc in line for inc in ['#include "spc.itp"', '#include "tip3p.itp"',
                                           '#include "amber99sb-ildn.ff/spce.itp"',
                                           '#include "amber99sb-ildn.ff/ions.itp"']):
                print(f"[DEBUG] Skipping water/ion include: {line.rstrip()}")
                continue
            filtered_lines.append(line)

    with open(output_top_file, 'w') as outfile:
        outfile.writelines(filtered_lines)

    print(f"[INFO] {output_top_file} created with only protein definitions.")

# prepare_system function (same as before)
def prepare_system(initial_pdb_source_path, workdir):
    """
    Prepare the molecular system for MM/GBSA.
    Corrects chain IDs and renumbers residues in the initial PDB, then runs full GROMACS setup,
    and prepares protein-only files and index for gmx_MMPBSA.
    """
    print(f"\n--- [STAGE] Preparing system in {workdir} ---")
    
    if os.path.exists(workdir):
        print(f"[INFO] Cleaning existing working directory: {workdir}")
        shutil.rmtree(workdir)
    os.makedirs(workdir, exist_ok=True)
    
    if not os.path.exists(initial_pdb_source_path):
        raise FileNotFoundError(f"Source PDB not found at: {initial_pdb_source_path}")

    model_pdb = os.path.join(workdir, "model.pdb")
    
    temp_pdb_for_relabel = os.path.join(workdir, "input_for_relabel.pdb")
    shutil.copy(initial_pdb_source_path, temp_pdb_for_relabel)

    actual_chain_a_length, actual_chain_b_length = relabel_and_renumber_pdb_chains(
        temp_pdb_for_relabel, model_pdb
    )
        
    run_cmd(f"{GROMACS} pdb2gmx -f {os.path.basename(model_pdb)} -o processed.gro -water spce -ff amber99sb-ildn", cwd=workdir)
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

    full_system_ndx = os.path.join(workdir, "full_system.ndx")
    cmd_make_ndx_full = [GROMACS, 'make_ndx', '-f', 'md.tpr', '-o', os.path.basename(full_system_ndx)]
    p_full = subprocess.Popen(cmd_make_ndx_full, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=workdir)
    commands_full = "a Protein\nq\n"
    out_full, err_full = p_full.communicate(commands_full)
    if p_full.returncode != 0:
        raise RuntimeError(f"make_ndx failed for full system in {workdir}:\n{err_full}")
    print("[INFO] Full system index file created with 'Protein' group.")

    run_cmd(f"echo 'Protein' | {GROMACS} trjconv -s md.tpr -f md.xtc -o protein_only.pdb -dump 0 -n {os.path.basename(full_system_ndx)}", cwd=workdir)

    with open(os.path.join(workdir, "protein_only.pdb")) as f:
        pdb_content = f.read()
        if any(res in pdb_content for res in ("SOL", "NA", "CL", "WAT", "HOH")):
            raise RuntimeError("protein_only.pdb still contains water or ions after trjconv!")

    original_topol = os.path.join(workdir, "topol.top")
    protein_topol = os.path.join(workdir, "protein_only.top")
    create_protein_only_topology(original_topol, protein_topol)

    empty_mdp = os.path.join(MDP_DIR, "empty.mdp")
    if not os.path.exists(empty_mdp):
        with open(empty_mdp, "w") as f:
            f.write("define = -DFLEXIBLE\nintegrator = md\nnsteps = 0\n")
    run_cmd(f"{GROMACS} grompp -f {empty_mdp} -c protein_only.pdb -p protein_only.top -o protein_only.tpr -maxwarn 1", cwd=workdir)

    run_cmd(f"echo 'Protein' | {GROMACS} trjconv -s md.tpr -f md.xtc -n {os.path.basename(full_system_ndx)} -o protein_only.xtc", cwd=workdir)

    run_cmd(f"echo 'Protein\nProtein' | {GROMACS} trjconv -s md.tpr -f md.xtc -o md_noPBC_protein.xtc -pbc mol -center -n {os.path.basename(full_system_ndx)}", cwd=workdir)
    run_cmd(f"echo 'Protein' | {GROMACS} trjconv -s md.tpr -f md.xtc -o last_frame_protein.pdb -dump -1 -n {os.path.basename(full_system_ndx)}", cwd=workdir)

    with open(os.path.join(workdir, "mmpbsa.in"), "w") as f:
        f.write(MMPBSA_IN)

    print(f"[INFO] System prepared for MM/GBSA in {workdir}")
    
    return actual_chain_a_length, actual_chain_b_length

# get_chain_residue_ranges function (same as before)
def get_chain_residue_ranges(pdb_path):
    """
    Parses a PDB file and returns a dict of residue ID ranges for each chain.
    """
    chain_residues = defaultdict(set)

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                chain = line[21].strip()
                resid = int(line[22:26])
                chain_residues[chain].add(resid)

    residue_ranges = {}
    for chain_id, resids in chain_residues.items():
        residue_ranges[chain_id] = (min(resids), max(resids))

    return residue_ranges

# --- REVISED create_chain_index_by_residues (Simplified and Corrected) ---
def create_chain_index_by_residues(tpr_file_path, ndx_file_path, chainA_range, chainB_range, cwd):
    """
    Creates an index file (.ndx) with groups for Chain A and Chain B
    by selecting atoms based on their residue numbers using gmx make_ndx.
    This version creates the groups and names them in a single interactive session.
    """
    print(f"[INFO] Creating index file {os.path.basename(ndx_file_path)} for Chain A (resid {chainA_range[0]}-{chainA_range[1]}) and Chain B (resid {chainB_range[0]}-{chainB_range[1]})")

    cmd = [GROMACS, 'make_ndx', '-f', os.path.basename(tpr_file_path), '-o', os.path.basename(ndx_file_path)]
    
    # We need to determine the default group numbers GROMACS assigns after `r` command
    # and then use `name` on those. This is slightly heuristic but usually reliable for 2 new groups.
    # We will simulate the interactive session to get the group numbers.
    
    # First, run make_ndx and just print groups to find existing highest number
    temp_check_ndx = os.path.join(cwd, "temp_check.ndx")
    p_check = subprocess.run([GROMACS, 'make_ndx', '-f', os.path.basename(tpr_file_path), '-o', os.path.basename(temp_check_ndx)], 
                             input="q\n", capture_output=True, text=True, cwd=cwd)
    
    highest_default_group_num = -1
    for line in p_check.stdout.splitlines():
        match = re.match(r'^\s*(\d+)\s+.*', line)
        if match:
            try:
                group_num = int(match.group(1))
                highest_default_group_num = max(highest_default_group_num, group_num)
            except ValueError:
                pass
    
    os.remove(temp_check_ndx) # Clean up temp file

    # Calculate the expected new group numbers
    chain_a_group_num = highest_default_group_num + 1
    chain_b_group_num = highest_default_group_num + 2

    print(f"[DEBUG] Expected Chain_A group number: {chain_a_group_num}, Chain_B group number: {chain_b_group_num}")

    commands = [
        f"r {chainA_range[0]}-{chainA_range[1]}", # Creates a new group
        f"r {chainB_range[0]}-{chainB_range[1]}", # Creates another new group
        f"name {chain_a_group_num} Chain_A", # Renames the first new group
        f"name {chain_b_group_num} Chain_B", # Renames the second new group
        "q" # Save and quit
    ]
    input_str = "\n".join(commands) + "\n"

    print(f"[DEBUG] make_ndx input:\n{input_str}")
    process = subprocess.run(cmd, input=input_str, capture_output=True, text=True, cwd=cwd)

    if process.returncode != 0:
        print(f"[FAIL] make_ndx failed in {cwd}:\nSTDOUT:\n{process.stdout}\nSTDERR:\n{process.stderr}")
        raise RuntimeError(f"make_ndx failed with exit code {process.returncode}")
    else:
        print("[INFO] Index file created successfully with Chain_A and Chain_B groups.")
        print(f"STDOUT:\n{process.stdout}")


# run_gmx_mmpbsa function (same as before)
def run_gmx_mmpbsa(mmpbsa_in, complex_tpr_file, xtc_file, topol_file, ndx_file):
    """Run gmx_MMPBSA using specified inputs and chain groups."""
    cmd = [
        GMXMMPBSA,
        '-O',
        '-i', mmpbsa_in,
        '-cs', complex_tpr_file,
        '-ct', xtc_file,
        '-cp', topol_file,
        '-ci', ndx_file,
        '-cg', 'Chain_A', 'Chain_B',
    ]
    print(f"[RUN] {' '.join(cmd)}")
    run_cwd = os.path.dirname(mmpbsa_in) 
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=run_cwd)
    if result.returncode != 0:
        print(f"[FAIL] gmx_MMPBSA failed in {run_cwd}:\n{result.stderr}")
        raise RuntimeError(f"gmx_MMPBSA failed with exit code {result.returncode}")
    else:
        print("[INFO] gmx_MMPBSA completed successfully.")
        print(result.stdout)


# Main Execution Block (UPDATED with your configuration paths)
if __name__ == "__main__":
    
    # Ensure MDP_DIR exists
    os.makedirs(MDP_DIR, exist_ok=True)
    # Create dummy .mdp files if they don't exist
    for mdp_file in ["ions.mdp", "minim.mdp", "nvt.mdp", "npt.mdp", "md.mdp", "empty.mdp"]: # Added empty.mdp here
        path = os.path.join(MDP_DIR, mdp_file)
        if not os.path.exists(path):
            print(f"[INFO] Creating dummy {mdp_file} in {MDP_DIR}. Please replace with actual MDPs.")
            with open(path, "w") as f:
                if mdp_file == "empty.mdp":
                    f.write("define = -DFLEXIBLE\nintegrator = md\nnsteps = 0\n")
                else:
                    f.write(f"; Dummy {mdp_file}\nintegrator = md\nnsteps = 1000 ; Adjust as needed for actual simulation\n")

    # --- WT System Processing ---
    wt_initial_pdb_source = os.path.join(WT_SOURCE_DIR, "model.pdb") # Assuming model.pdb is directly in WT_SOURCE_DIR
    wt_workdir = WT_RUN_DIR

    wt_chain_a_length, wt_chain_b_length = prepare_system(wt_initial_pdb_source, wt_workdir)

    wt_protein_only_pdb_path = os.path.join(wt_workdir, "protein_only.pdb")
    wt_chain_residue_ranges = get_chain_residue_ranges(wt_protein_only_pdb_path)
    
    print(f"[INFO] WT Detected residue ranges from protein_only.pdb: {wt_chain_residue_ranges}")

    wt_tpr_for_gmxmmpbsa_cs = os.path.join(wt_workdir, "protein_only.tpr")
    wt_xtc_for_gmxmmpbsa_ct = os.path.join(wt_workdir, "protein_only.xtc")
    wt_protein_only_topol_for_mmpbsa = os.path.join(wt_workdir, "protein_only.top")
    wt_ndx = os.path.join(wt_workdir, "index_chain.ndx")
    wt_mmpbsa_input = os.path.join(wt_workdir, "mmpbsa.in")

    create_chain_index_by_residues(
        tpr_file_path=wt_tpr_for_gmxmmpbsa_cs,
        ndx_file_path=wt_ndx,
        chainA_range=wt_chain_residue_ranges.get('A', (1,1)),
        chainB_range=wt_chain_residue_ranges.get('B', (1,1)),
        cwd=wt_workdir
    )

    run_gmx_mmpbsa(wt_mmpbsa_input,
                   wt_tpr_for_gmxmmpbsa_cs,
                   wt_xtc_for_gmxmmpbsa_ct,
                   wt_protein_only_topol_for_mmpbsa,
                   wt_ndx)

    print(f"\n--- [STAGE] WT MM/GBSA calculation complete ---")


    # --- MT System Processing (similar structure) ---
    mt_initial_pdb_source = os.path.join(MT_SOURCE_DIR, "model.pdb")
    mt_workdir = MT_RUN_DIR

    mt_chain_a_length, mt_chain_b_length = prepare_system(mt_initial_pdb_source, mt_workdir)

    mt_protein_only_pdb_path = os.path.join(mt_workdir, "protein_only.pdb")
    mt_chain_residue_ranges = get_chain_residue_ranges(mt_protein_only_pdb_path)
    
    print(f"[INFO] MT Detected residue ranges from protein_only.pdb: {mt_chain_residue_ranges}")

    mt_tpr_for_gmxmmpbsa_cs = os.path.join(mt_workdir, "protein_only.tpr")
    mt_xtc_for_gmxmmpbsa_ct = os.path.join(mt_workdir, "protein_only.xtc")
    mt_protein_only_topol_for_mmpbsa = os.path.join(mt_workdir, "protein_only.top")
    mt_ndx = os.path.join(mt_workdir, "index_chain.ndx")
    mt_mmpbsa_input = os.path.join(mt_workdir, "mmpbsa.in")

    create_chain_index_by_residues(
        tpr_file_path=mt_tpr_for_gmxmmpbsa_cs,
        ndx_file_path=mt_ndx,
        chainA_range=mt_chain_residue_ranges.get('A', (1,1)),
        chainB_range=mt_chain_residue_ranges.get('B', (1,1)),
        cwd=mt_workdir
    )

    run_gmx_mmpbsa(mt_mmpbsa_input,
                   mt_tpr_for_gmxmmpbsa_cs,
                   mt_xtc_for_gmxmmpbsa_ct,
                   mt_protein_only_topol_for_mmpbsa,
                   mt_ndx)

    print(f"\n--- [STAGE] MT MM/GBSA calculation complete ---")

    print("\n--- All calculations finished. Review gmx_MMPBSA.log files for results. ---")