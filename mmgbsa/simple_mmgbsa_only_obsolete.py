#!/usr/bin/env python3
import os
import subprocess
import re
from collections import defaultdict

# ────────────────────────────────────────────────
# PATHS
# ────────────────────────────────────────────────
GROMACS_2025 = "/opt/gromacs/bin/gmx"
GMXMMPBSA_BIN = os.path.expanduser("~/miniconda3/envs/gmxMMPBSA/bin/gmx_MMPBSA")

# ────────────────────────────────────────────────
# UTILS
# ────────────────────────────────────────────────
def run_cmd(cmd_list, cwd):
    """Run a command as list, show output, and stop if fails."""
    print(f"[RUN] {' '.join(cmd_list)}")
    result = subprocess.run(cmd_list, cwd=cwd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"Command failed: {' '.join(cmd_list)}")

def get_chain_residue_ranges(pdb_path):
    """Detect residue ID ranges per chain from a PDB file."""
    chain_residues = defaultdict(set)
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                chain = line[21].strip()
                resid = int(line[22:26])
                chain_residues[chain].add(resid)
    return {c: (min(r), max(r)) for c, r in chain_residues.items()}

# ────────────────────────────────────────────────
# INDEX CREATION (with GROMACS 2025)
# ────────────────────────────────────────────────
def create_chain_index_by_residues(tpr_file, ndx_file, chainA_range, chainB_range, cwd):
    """Use gmx 2025.1 to create index groups for MMPBSA."""
    print(f"[INFO] Creating index file {os.path.basename(ndx_file)} using GROMACS 2025.1")

    # First check existing groups
    temp_check = os.path.join(cwd, "temp_check.ndx")
    check = subprocess.run(
        [GROMACS_2025, "make_ndx", "-f", os.path.basename(tpr_file), "-o", os.path.basename(temp_check)],
        input="q\n",
        capture_output=True,
        text=True,
        cwd=cwd
    )
    ndx_output = check.stdout + "\n" + check.stderr

    highest_group = -1
    for line in ndx_output.splitlines():
        match = re.match(r"^\s*(\d+)\s+.*", line)
        if match:
            highest_group = max(highest_group, int(match.group(1)))

    if os.path.exists(temp_check):
        os.remove(temp_check)
    if highest_group < 0:
        print(ndx_output)
        raise RuntimeError("[ERROR] No groups found in TPR. Check your file.")

    # Build commands
    commands = []
    if chainA_range:
        commands += [f"r {chainA_range[0]}-{chainA_range[1]}", f"name {highest_group+1} Chain_A"]
    if chainB_range:
        commands += [f"r {chainB_range[0]}-{chainB_range[1]}", f"name {highest_group+2} Chain_B"]
    commands.append("q")

    process = subprocess.run(
        [GROMACS_2025, "make_ndx", "-f", os.path.basename(tpr_file), "-o", os.path.basename(ndx_file)],
        input="\n".join(commands) + "\n",
        capture_output=True,
        text=True,
        cwd=cwd
    )
    if process.returncode != 0:
        print(process.stdout)
        print(process.stderr)
        raise RuntimeError("[FAIL] make_ndx failed.")
    print("[INFO] Index file created successfully.")

# ────────────────────────────────────────────────
# RUN GMXMMPBSA
# ────────────────────────────────────────────────
def run_gmx_mmpbsa(mmpbsa_in, tpr, xtc, top, ndx, cwd):
    """Run gmx_MMPBSA using its own binary, with TPR from GROMACS 2025.1."""
    cmd = [
        GMXMMPBSA_BIN,
        "-O",
        "-i", mmpbsa_in,
        "-cs", tpr,
        "-ct", xtc,
        "-cp", top,
        "-ci", ndx,
        "-cg", "Chain_A", "Chain_B"
    ]
    print(f"[RUN] {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError("gmx_MMPBSA failed.")
    print(result.stdout)
    print("[INFO] gmx_MMPBSA completed successfully.")

# ────────────────────────────────────────────────
# MAIN WORKFLOW
# ────────────────────────────────────────────────
if __name__ == "__main__":
    BASEDIR = "/home2/mmGBSA/test2"
    RUN_DIRS = ["WT_run_auto", "MT_run_auto"]

    for run_name in RUN_DIRS:
        WORKDIR = os.path.join(BASEDIR, run_name)
        print(f"\n==============================")
        print(f"[INFO] Starting MMGBSA run for: {run_name}")
        print(f"==============================")

        PDB = os.path.join(WORKDIR, "protein_only.pdb")
        TPR = os.path.join(WORKDIR, "protein_only.tpr")  # using 2025.1 TPR
        XTC = os.path.join(WORKDIR, "protein_only.xtc")
        TOP = os.path.join(WORKDIR, "protein_only.top")
        NDX = os.path.join(WORKDIR, "index_chain.ndx")
        MMPBSA_IN = os.path.join(WORKDIR, "mmpbsa.in")

        # Step 1: get chain ranges from pdb
        chain_ranges = get_chain_residue_ranges(PDB)
        print(f"[INFO] Chain residue ranges detected: {chain_ranges}")

        # Step 2: create index using gmx 2025.1
        create_chain_index_by_residues(
            tpr_file=TPR,
            ndx_file=NDX,
            chainA_range=chain_ranges.get('A', (1, 1)),
            chainB_range=chain_ranges.get('B', (1, 1)),
            cwd=WORKDIR
        )

        # Step 3: run MMPBSA
        try:
            run_gmx_mmpbsa(
                mmpbsa_in=MMPBSA_IN,
                tpr=TPR,
                xtc=XTC,
                top=TOP,
                ndx=NDX,
                cwd=WORKDIR
            )
            print(f"[SUCCESS] MMPBSA completed for {run_name}")
        except RuntimeError as e:
            print(f"[ERROR] MMPBSA failed for {run_name}: {e}")
