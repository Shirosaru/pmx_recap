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
    """
    Fix PDB chain IDs by splitting residues at split_resid:
    residues <= split_resid get chain A, else chain B.
    """
    with open(input_pdb) as inp, open(output_pdb, "w") as out:
        for line in inp:
            if line.startswith(("ATOM", "HETATM")):
                resid = int(line[22:26])
                chain = 'A' if resid <= split_resid else 'B'
                new_line = line[:21] + chain + line[22:]
                out.write(new_line)
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

    # Run GROMACS setup commands
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
    run_cmd(f"echo -e '1\n0' | {GROMACS} trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center", cwd=workdir)
    run_cmd(f"echo 1 | {GROMACS} trjconv -s md.tpr -f md.xtc -o last_frame.pdb -dump -1", cwd=workdir)

    # Write the MM/GBSA input file
    with open(os.path.join(workdir, "mmpbsa.in"), "w") as f:
        f.write(MMPBSA_IN)
    print(f"[INFO] mmpbsa.in created.")

    # Split fixed PDB by chain
    split_pdb_by_chain(model_pdb, workdir)

def create_chain_index(tpr_file, ndx_file):
    """Create index file with chain A and B groups (groups 17 and 18)."""
    cmd = ['gmx', 'make_ndx', '-f', tpr_file, '-o', ndx_file]
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    commands = "\n".join([
        "r Chain_A",
        "name 17 Chain_A",
        "r Chain_B",
        "name 18 Chain_B",
        "q"
    ]) + "\n"

    out, err = p.communicate(commands)
    if p.returncode != 0:
        raise RuntimeError(f"make_ndx failed:\n{err}")
    print("[INFO] Index file created with Chain_A and Chain_B groups.")

def run_gmx_mmpbsa(mmpbsa_in, tpr_file, xtc_file, topol_file, ndx_file, chainA_pdb, chainB_pdb):
    """Run gmx_MMPBSA using specified inputs and chain groups."""
    cmd = [
        GMXMMPBSA,
        '-O',
        '-i', mmpbsa_in,
        '-cs', tpr_file,
        '-ct', xtc_file,
        '-cp', topol_file,
        '-ci', ndx_file,
        '-cg', '17', '18',   # ❌ old, unsafe
        '-rs', chainA_pdb,
        '-ls', chainB_pdb
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[FAIL] gmx_MMPBSA failed:\n{result.stderr}")
    else:
        print("[INFO] gmx_MMPBSA completed successfully.")
        print(result.stdout)

if __name__ == "__main__":
    # Adjust these paths as needed
    workdir = WT_DIR  # or MT_DIR depending on your run
    prepare_system(workdir)

    tpr = os.path.join(workdir, "md.tpr")
    ndx = os.path.join(workdir, "index_chain.ndx")
    mmpbsa_input = os.path.join(workdir, "mmpbsa.in")
    xtc = os.path.join(workdir, "md_noPBC.xtc")
    topol = os.path.join(workdir, "topol.top")
    chainA_pdb = os.path.join(workdir, "Chain_A.pdb")
    chainB_pdb = os.path.join(workdir, "Chain_B.pdb")

    create_chain_index(tpr, ndx)
    run_gmx_mmpbsa(mmpbsa_input, tpr, xtc, topol, ndx, chainA_pdb, chainB_pdb)
