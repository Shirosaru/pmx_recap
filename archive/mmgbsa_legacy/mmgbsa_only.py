# ==========================
# mmgbsa_only.py (CPU PART)
# ==========================
import os
from simple_mmgbsa_8 import (
    get_chain_residue_ranges,
    create_chain_index_by_residues,
    get_zinc_resids_from_pdb,
    run_gmx_mmpbsa,
    WT_RUN_DIR,
    MT_RUN_DIR,
)

if __name__ == "__main__":

    # ------------------------
    # WT MMGBSA
    # ------------------------
    wt_workdir = WT_RUN_DIR
    try:
        wt_pdb = os.path.join(wt_workdir, "protein_only.pdb")
        wt_chain_ranges = get_chain_residue_ranges(wt_pdb)

        zinc_resids = get_zinc_resids_from_pdb(wt_pdb)

        chain_a_group_num, chain_b_group_num = create_chain_index_by_residues(
            tpr_file_path=os.path.join(wt_workdir, "md.tpr"),
            ndx_file_path=os.path.join(wt_workdir, "index_chain.ndx"),
            chainA_range=wt_chain_ranges.get("A", (1, 1)),
            chainB_range=wt_chain_ranges.get("B", (1, 1)),
            zinc_resids=zinc_resids,
            cwd=wt_workdir
        )

        run_gmx_mmpbsa(
            mmpbsa_in=os.path.join(wt_workdir, "mmpbsa.in"),
            complex_tpr_file=os.path.join(wt_workdir, "md.tpr"),
            xtc_file=os.path.join(wt_workdir, "md.xtc"),
            topol_file=os.path.join(wt_workdir, "protein_only.top"),
            ndx_file=os.path.join(wt_workdir, "index_chain.ndx"),
            chain_a_group_num=chain_a_group_num,
            chain_b_group_num=chain_b_group_num
        )

        print("\n--- [STAGE] WT MMGBSA complete ---\n")

    except Exception as e:
        print(f"[ERROR] WT MMGBSA failed: {e}")

    # ------------------------
    # MT MMGBSA
    # ------------------------
    mt_workdir = MT_RUN_DIR
    try:
        mt_pdb = os.path.join(mt_workdir, "protein_only.pdb")
        mt_chain_ranges = get_chain_residue_ranges(mt_pdb)

        zinc_resids = get_zinc_resids_from_pdb(mt_pdb)

        chain_a_group_num, chain_b_group_num = create_chain_index_by_residues(
            tpr_file_path=os.path.join(mt_workdir, "md.tpr"),
            ndx_file_path=os.path.join(mt_workdir, "index_chain.ndx"),
            chainA_range=mt_chain_ranges.get("A", (1, 1)),
            chainB_range=mt_chain_ranges.get("B", (1, 1)),
            zinc_resids=zinc_resids,
            cwd=mt_workdir
        )

        run_gmx_mmpbsa(
            mmpbsa_in=os.path.join(mt_workdir, "mmpbsa.in"),
            complex_tpr_file=os.path.join(mt_workdir, "md.tpr"),
            xtc_file=os.path.join(mt_workdir, "md.xtc"),
            topol_file=os.path.join(mt_workdir, "protein_only.top"),
            ndx_file=os.path.join(mt_workdir, "index_chain.ndx"),
            chain_a_group_num=chain_a_group_num,
            chain_b_group_num=chain_b_group_num
        )

        print("\n--- [STAGE] MT MMGBSA complete ---\n")

    except Exception as e:
        print(f"[ERROR] MT MMGBSA failed: {e}")

    print("\n--- All MMGBSA calculations finished. ---")
