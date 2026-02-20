# ==========================
# md_only.py  (GPU PART)
# ==========================
import os
from simple_mmgbsa_8 import (
    prepare_system,
    get_chain_residue_ranges,
    create_chain_index_by_residues,
    get_zinc_resids_from_pdb,
    MDP_DIR,
    WT_SOURCE_DIR,
    WT_RUN_DIR,
    MT_SOURCE_DIR,
    MT_RUN_DIR,
)

if __name__ == "__main__":

    # Ensure MDP_DIR exists
    os.makedirs(MDP_DIR, exist_ok=True)

    # Create dummy .mdp files if not present
    for mdp_file in ["ions.mdp", "minim.mdp", "nvt.mdp", "npt.mdp", "md.mdp", "empty.mdp"]:
        path = os.path.join(MDP_DIR, mdp_file)
        if not os.path.exists(path):
            print(f"[INFO] Creating dummy {mdp_file} in {MDP_DIR}")
            with open(path, "w") as f:
                if mdp_file == "empty.mdp":
                    f.write("define = -DFLEXIBLE\nintegrator = md\nnsteps = 0\n")
                else:
                    f.write(f"; Dummy {mdp_file}\nintegrator = md\nnsteps = 1000\n")

    # -------------------------------
    # WT MD WORKFLOW
    # -------------------------------
    wt_initial_pdb_source = os.path.join(WT_SOURCE_DIR, "WTModel_1.pdb")
    wt_workdir = WT_RUN_DIR

    try:
        prepare_system(wt_initial_pdb_source, wt_workdir)

        wt_protein_only_pdb = os.path.join(wt_workdir, "protein_only.pdb")
        wt_chain_residue_ranges = get_chain_residue_ranges(wt_protein_only_pdb)
        print(f"[INFO] WT residue ranges: {wt_chain_residue_ranges}")

        zinc_resids = get_zinc_resids_from_pdb(wt_protein_only_pdb)

        # Create chain index (needed later for MMGBSA)
        create_chain_index_by_residues(
            tpr_file_path=os.path.join(wt_workdir, "md.tpr"),
            ndx_file_path=os.path.join(wt_workdir, "index_chain.ndx"),
            chainA_range=wt_chain_residue_ranges.get("A", (1, 1)),
            chainB_range=wt_chain_residue_ranges.get("B", (1, 1)),
            zinc_resids=zinc_resids,
            cwd=wt_workdir
        )

        print("\n--- [STAGE] WT MD complete ---\n")

    except Exception as e:
        print(f"[ERROR] WT MD failed: {e}")

    # -------------------------------
    # MT MD WORKFLOW
    # -------------------------------
    mt_initial_pdb_source = os.path.join(MT_SOURCE_DIR, "MutModel_1.pdb")
    mt_workdir = MT_RUN_DIR

    try:
        prepare_system(mt_initial_pdb_source, mt_workdir)

        mt_protein_only_pdb = os.path.join(mt_workdir, "protein_only.pdb")
        mt_chain_residue_ranges = get_chain_residue_ranges(mt_protein_only_pdb)
        print(f"[INFO] MT residue ranges: {mt_chain_residue_ranges}")

        zinc_resids = get_zinc_resids_from_pdb(mt_protein_only_pdb)

        create_chain_index_by_residues(
            tpr_file_path=os.path.join(mt_workdir, "md.tpr"),
            ndx_file_path=os.path.join(mt_workdir, "index_chain.ndx"),
            chainA_range=mt_chain_residue_ranges.get("A", (1, 1)),
            chainB_range=mt_chain_residue_ranges.get("B", (1, 1)),
            zinc_resids=zinc_resids,
            cwd=mt_workdir
        )

        print("\n--- [STAGE] MT MD complete ---\n")

    except Exception as e:
        print(f"[ERROR] MT MD failed: {e}")

    print("\n--- All MD finished. Now run mmgbsa_only.py on CPU. ---")
