import pandas as pd
import os
import re

# --- Configuration ---
BASE_DIR = "/home2/mmGBSA/test2"
WT_RUN_DIR = os.path.join(BASE_DIR, "WT_run_auto")
MT_RUN_DIR = os.path.join(BASE_DIR, "MT_run_auto")


def read_delta_total_from_mmpbsa_file(filepath):
    """
    Parses a FINAL_RESULTS_MMPBSA.dat file to find the ΔTOTAL binding energy.
    This function is tailored to the specific file format provided in the debugging output.
    
    Args:
        filepath (str): Path to the FINAL_RESULTS_MMPBSA.dat file.
        
    Returns:
        tuple: A tuple containing (delta_total_avg, delta_total_sd) in kcal/mol,
               or (None, None) if the file could not be parsed.
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            
        start_parsing = False
        
        for i, line in enumerate(lines):
            # Step 1: Find the anchor line for the Delta results
            if "Delta (Complex - Receptor - Ligand):" in line:
                start_parsing = True
                # The data table starts a few lines after this.
                # Let's start looking from the next line.
                continue
            
            # Step 2: Once we've found the anchor, look for the ΔTOTAL line
            if start_parsing:
                # The ΔTOTAL line is also a consistent marker
                if line.strip().startswith("ΔTOTAL"):
                    # Use a regular expression to split the line by one or more spaces
                    parts = re.split(r'\s+', line.strip())
                    
                    # The parts are expected to be: ['ΔTOTAL', 'Average', 'SD(Prop.)', 'SD', ...]
                    # The value we want is the 'Average' which is the second item.
                    # The standard deviation is the 'SD' value which is the fourth item.
                    
                    # Ensure we have enough parts to avoid an IndexError
                    if len(parts) >= 4:
                        try:
                            delta_total_avg = float(parts[1])
                            # SD is the 4th value from the start (index 3)
                            delta_total_sd = float(parts[3])
                            return delta_total_avg, delta_total_sd
                        except (ValueError, IndexError):
                            print(f"WARNING: Found 'ΔTOTAL' line but could not parse values in {filepath}")
                            return None, None
                            
        # If we get here, the section or the ΔTOTAL line was not found
        print(f"WARNING: Could not find 'ΔTOTAL' binding energy in {filepath}")
        return None, None
        
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return None, None
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None, None


if __name__ == "__main__":
    
    result_file_name = "FINAL_RESULTS_MMPBSA.dat"
    
    # --- Process WT Results ---
    print("\n--- WT Results ---")
    wt_results_file_path = os.path.join(WT_RUN_DIR, result_file_name)
    delta_g_wt, sd_wt = read_delta_total_from_mmpbsa_file(wt_results_file_path)
    if delta_g_wt is not None:
        print(f"ΔG_WT: {delta_g_wt:.2f} kcal/mol (SD: {sd_wt:.2f} kcal/mol)")
    else:
        print("Failed to get WT results.")

    # --- Process MT Results ---
    print("\n--- MT Results ---")
    mt_results_file_path = os.path.join(MT_RUN_DIR, result_file_name)
    delta_g_mt, sd_mt = read_delta_total_from_mmpbsa_file(mt_results_file_path)
    if delta_g_mt is not None:
        print(f"ΔG_MT: {delta_g_mt:.2f} kcal/mol (SD: {sd_mt:.2f} kcal/mol)")
    else:
        print("Failed to get MT results.")

    # --- Calculate Delta Delta G ---
    print("\n--- Delta Delta G Calculation ---")
    if delta_g_wt is not None and delta_g_mt is not None:
        delta_delta_g = delta_g_mt - delta_g_wt
        
        # Calculate the propagated standard deviation
        delta_delta_g_sd = (sd_wt**2 + sd_mt**2)**0.5
        
        print(f"ΔG_WT: {delta_g_wt:.2f} ± {sd_wt:.2f} kcal/mol")
        print(f"ΔG_MT: {delta_g_mt:.2f} ± {sd_mt:.2f} kcal/mol")
        print(f"ΔΔG (ΔG_MT - ΔG_WT): {delta_delta_g:.2f} ± {delta_delta_g_sd:.2f} kcal/mol")
        print("\nA positive ΔΔG typically indicates that the mutation destabilizes the binding.")
        print("A negative ΔΔG typically indicates that the mutation stabilizes the binding.")
    else:
        print("Could not calculate Delta Delta G. Ensure both WT and MT results were parsed correctly.")