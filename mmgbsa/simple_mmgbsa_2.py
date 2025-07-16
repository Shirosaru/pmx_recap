import os
import subprocess
import shutil
import glob # Import glob for pattern matching

# Directories
BASE_DIR = "/home2/mmGBSA"
WT_DIR = os.path.join(BASE_DIR, "WT_run")
MT_DIR = os.path.join(BASE_DIR, "MT_run")
MDP_DIR = os.path.join(BASE_DIR, "mdp")

GROMACS = "gmx"          # Or full path to gmx executable
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
    print(f"[INFO] Fixing PDB chains for {input_pdb}...")
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
    print(f"[INFO] Chains fixed and saved to {output_pdb}")


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

    # CRITICAL FIX HERE:
    # Based on the debug output, your protein molecules are named 'Protein_chain_A' and 'Protein_chain_A2'
    protein_molecule_names = ["Protein_chain_A", "Protein_chain_A2"]

    # Also, keep includes for force fields and atom types
    # This list should be comprehensive for your force field
    keep_includes = ['#include "amber99sb-ildn.ff/forcefield.itp"',
                     '#include "amber99sb-ildn.ff/ffnonbonded.itp"',
                     '#include "amber99sb-ildn.ff/ffbonded.itp"',
                     '#include "amber99sb-ildn.ff/aminoacids.rtp"',
                     '#include "amber99sb-ildn.ff/dna.rtp"',
                     '#include "amber99sb-ildn.ff/rna.rtp"',
                     '#include "amber99sb-ildn.ff/atomtypes.atp"',
                     '#include "amber99sb-ildn.ff/residuetypes.dat"']

    for line in lines:
        stripped_line = line.strip()

        # Keep includes that are part of the force field itself
        if any(inc in line for inc in keep_includes):
            filtered_lines.append(line)
            continue

        if stripped_line == '[ molecules ]':
            in_molecules_section = True
            filtered_lines.append(line) # Keep the [ molecules ] header
            continue

        if in_molecules_section:
            if not stripped_line or stripped_line.startswith(';'): # Keep blank lines and comments
                filtered_lines.append(line)
            else:
                parts = stripped_line.split()
                # Only keep lines where the first part (molecule name) is in our protein_molecule_names list
                if len(parts) >= 2 and parts[0] in protein_molecule_names:
                    filtered_lines.append(line) # Keep the protein molecule line
                    print(f"[DEBUG] KEEPING molecule line: {line.rstrip()}") # Added for explicit debug
                else:
                    print(f"[DEBUG] Skipping molecule line in [ molecules ]: {line.rstrip()}")
        else:
            # Filter out specific includes for water/ions BEFORE [ molecules ] section
            if any(inc in line for inc in ['#include "spc.itp"', '#include "tip3p.itp"',
                                           '#include "amber99sb-ildn.ff/spce.itp"',
                                           '#include "amber99sb-ildn.ff/ions.itp"']):
                print(f"[DEBUG] Skipping water/ion include: {line.rstrip()}")
                continue
            filtered_lines.append(line)

    with open(output_top_file, 'w') as outfile:
        outfile.writelines(filtered_lines)

    print(f"[INFO] {output_top_file} created with only protein definitions.")


def clean_workdir(workdir):
    """Removes GROMACS backup files and other intermediate files."""
    print(f"[INFO] Cleaning up working directory: {workdir}")
    # List of file patterns to remove
    patterns_to_remove = [
        '#*#',          # GROMACS backup files (e.g., #topol.top.1#)
        '*.log',        # GROMACS log files (em.log, nvt.log, md.log etc.)
        '*.trr',        # Trajectories
        '*.xtc',        # Trajectories (md.xtc is kept if it's the final output you want to preserve for MMPBSA)
        '*.edr',        # Energy files
        '*.cpt',        # Checkpoint files
        '*.tpr',        # TPR files (em.tpr, nvt.tpr, npt.tpr, ions.tpr - protein_only.tpr is re-generated)
        '*.gro',        # GRO files (processed.gro, boxed.gro, solvated.gro, solv_ions.gro, em.gro, nvt.gro, npt.gro)
        '*.pdb',        # PDB files (except initial model.pdb and fixed one, which are copied/re-generated)
                        # Be careful if you have other PDBs you want to keep.
        '_GMXMMPBSA_*', # Intermediate files created by gmx_MMPBSA
        'gmx_MMPBSA.log', # Main gmx_MMPBSA log
        'protein_only.top', # Will be recreated
        'index_chain.ndx', # Will be recreated
        'mmpbsa.in' # Will be recreated
    ]

    # Specific files to explicitly keep or handle differently
    files_to_keep = [
        'model.pdb', # Original input PDB
        'model_fixed.pdb', # Initial fixed PDB (copied over model.pdb anyway)
        'md.tpr',    # The final full system TPR for MMPBSA input
        'md.xtc',    # The final full system XTC for MMPBSA input
        'last_frame.pdb' # If you want to keep the final frame PDB
    ]

    for pattern in patterns_to_remove:
        for f in glob.glob(os.path.join(workdir, pattern)):
            if os.path.basename(f) not in files_to_keep:
                try:
                    os.remove(f)
                    print(f"  Removed: {f}")
                except OSError as e:
                    print(f"  Error removing {f}: {e}")
    print(f"[INFO] Cleanup complete for {workdir}.")




def create_chain_index_by_residues(tpr_file_path, ndx_file_path, chainA_range, chainB_range, cwd):
    # CRITICAL CHANGE: Use protein_only.tpr for make_ndx to ensure atom numbering matches
    # the protein-only context. This is the TPR generated with protein_only.pdb and protein_only.top.
    protein_only_tpr_path = os.path.join(cwd, "protein_only.tpr")

    cmd = [GROMACS, 'make_ndx', '-f', os.path.basename(protein_only_tpr_path), '-o', os.path.basename(ndx_file_path)]
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=cwd)

    chain_a_cmd = f"r {chainA_range[0]}-{chainA_range[1]}"
    chain_b_cmd = f"r {chainB_range[0]}-{chainB_range[1]}"

    commands = "\n".join([
        chain_a_cmd,
        "name 17 Chain_A", # Ensure these group numbers and names are consistent with gmx_MMPBSA expected
        chain_b_cmd,
        "name 18 Chain_B",
        "q"
    ]) + "\n"

    out, err = p.communicate(commands)
    if p.returncode != 0:
        raise RuntimeError(f"make_ndx failed in {cwd}:\n{err}")
    print("[INFO] Index file created with Chain_A and Chain_B based on residue numbers for protein_only.tpr.")



def prepare_system(workdir):
    """Prepare the molecular system: fix PDB, generate files, run GROMACS steps."""
    print(f"[STAGE] Preparing system in {workdir}")
    os.makedirs(workdir, exist_ok=True)

    clean_workdir(workdir) # Clean up first

    model_pdb = os.path.join(workdir, "model.pdb")
    fixed_pdb = os.path.join(workdir, "model_fixed.pdb")

    if not os.path.exists(model_pdb):
        print(f"[ERROR] Required file {model_pdb} not found. Please place your initial PDB there.")
        exit(1)

    fix_pdb_chains(model_pdb, fixed_pdb, split_resid=220)
    shutil.copy(fixed_pdb, model_pdb)

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

    # --- NEW / MODIFIED STEPS FOR GMX_MMPBSA INPUTS ---

    # 1. Create a protein-only PDB from the last frame of the full MD
    # This 'protein_only.pdb' will be used to generate the protein_only.tpr
    run_cmd(f"echo 'Protein' | {GROMACS} trjconv -s md.tpr -f md.xtc -o protein_only.pdb -dump 0", cwd=workdir)

    # 2. Create the protein-only topology file (from full topol.top)
    original_topol = os.path.join(workdir, "topol.top")
    protein_only_topol = os.path.join(workdir, "protein_only.top")
    create_protein_only_topology(original_topol, protein_only_topol)

    # 3. Create a minimal empty.mdp for grompp (if not already existing)
    empty_mdp_path = os.path.join(MDP_DIR, "empty.mdp")
    if not os.path.exists(empty_mdp_path):
        with open(empty_mdp_path, "w") as f:
            f.write("define = -DFLEXIBLE\nintegrator = md\nnsteps = 0\n")
        print(f"[INFO] Created minimal empty.mdp for protein-only TPR generation: {empty_mdp_path}")

    # 4. Generate the protein-only TPR. This TPR defines the atoms and groups for the protein-only system.
    run_cmd(f"{GROMACS} grompp -f {empty_mdp_path} -c protein_only.pdb -p protein_only.top -o protein_only.tpr -maxwarn 1", cwd=workdir)

    # 5. Generate the protein-only XTC trajectory.
    #    It's CRITICAL that this trajectory corresponds to the protein_only.tpr.
    #    We use the 'Protein' group to extract only protein atoms from the full trajectory.
    run_cmd(f"echo 'Protein' | {GROMACS} trjconv -s md.tpr -f md.xtc -o protein_only.xtc", cwd=workdir)


    # --- Other trjconv commands (unchanged, not directly used by gmx_MMPBSA inputs) ---
    run_cmd(f"echo '1\n0' | {GROMACS} trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center", cwd=workdir)
    run_cmd(f"echo 1 | {GROMACS} trjconv -s md.tpr -f md.xtc -o last_frame.pdb -dump -1", cwd=workdir)

    # Write the MM/GBSA input file
    with open(os.path.join(workdir, "mmpbsa.in"), "w") as f:
        f.write(MMPBSA_IN)
    print(f"[INFO] mmpbsa.in created.")

    # Split fixed PDB by chain
    split_pdb_by_chain(model_pdb, workdir)


def run_gmx_mmpbsa(mmpbsa_in, complex_tpr_file, xtc_file, topol_file, ndx_file):
    """Run gmx_MMPBSA using specified inputs and chain groups.
    CRITICAL: Now use the protein-only TPR and XTC for gmx_MMPBSA.
    """
    cmd = [
        GMXMMPBSA,
        '-O',
        '-i', mmpbsa_in,
        '-cs', complex_tpr_file,      # This should now be protein_only.tpr
        '-ct', xtc_file,              # This should now be protein_only.xtc
        '-cp', topol_file,            # This is already protein_only.top
        '-ci', ndx_file,
        '-cg', 'Chain_A', 'Chain_B',
    ]
    print(f"[RUN] {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(mmpbsa_in))
    if result.returncode != 0:
        print(f"[FAIL] gmx_MMPBSA failed. Return code: {result.returncode}")
        print(f"STDOUT:\n{result.stdout}")
        print(f"STDERR:\n{result.stderr}")
        print(f"Please also check the gmx_MMPBSA.log file in {os.path.dirname(mmpbsa_in)} for more details.")
        exit(1)
    else:
        print("[INFO] gmx_MMPBSA completed successfully.")
        print(result.stdout)


if __name__ == "__main__":
    workdir = WT_DIR

    prepare_system(workdir)

    # --- Files for gmx_MMPBSA ---
    # Now explicitly use the protein-only TPR and XTC here
    tpr_protein_only = os.path.join(workdir, "protein_only.tpr")
    xtc_protein_only = os.path.join(workdir, "protein_only.xtc")

    protein_only_topol_for_mmpbsa = os.path.join(workdir, "protein_only.top") # This is already correct

    ndx = os.path.join(workdir, "index_chain.ndx")
    mmpbsa_input = os.path.join(workdir, "mmpbsa.in")

    # Create the index file. This step is still correct as it uses protein_only.tpr
    create_chain_index_by_residues(
        tpr_file_path=tpr_protein_only, # Still good to pass this for clarity, although function uses explicit path
        ndx_file_path=ndx,
        chainA_range=(1, 180),
        chainB_range=(181, 431),
        cwd=workdir
    )

    # Call gmx_MMPBSA with the newly created protein-only files
    run_gmx_mmpbsa(mmpbsa_input, tpr_protein_only, xtc_protein_only, protein_only_topol_for_mmpbsa, ndx)