import os
import re

# --- Configuration ---
BASE_DIR = "/home2/mmGBSA/test2"
WT_RUN_DIR = os.path.join(BASE_DIR, "WT_run_auto")
MT_RUN_DIR = os.path.join(BASE_DIR, "MT_run_auto")


def debug_mmpbsa_file_to_file(filepath, output_file):
    """
    Writes the contents of the file with line numbers to a specified output file.
    
    Args:
        filepath (str): Path to the FINAL_RESULTS_MMPBSA.dat file.
        output_file (str): Path to the file where debugging output will be written.
    """
    
    with open(output_file, 'a') as out_f: # 'a' for append mode
        out_f.write(f"\n--- DEBUGGING FILE: {filepath} ---\n")
        
        try:
            with open(filepath, 'r') as f:
                for i, line in enumerate(f):
                    out_f.write(f"Line {i+1:3d}: {line}") # Write the line with number
        except FileNotFoundError:
            out_f.write(f"Error: File not found at {filepath}\n")
            out_f.write("--------------------------------------------------\n")
            return
        except Exception as e:
            out_f.write(f"Error reading file at {filepath}: {e}\n")
        
        out_f.write("--------------------------------------------------\n")


if __name__ == "__main__":
    
    result_file_name = "FINAL_RESULTS_MMPBSA.dat"
    output_debug_file = os.path.join(BASE_DIR, "mmpbsa_debug_output.txt")
    
    # Clear the previous debug output file if it exists
    if os.path.exists(output_debug_file):
        os.remove(output_debug_file)

    # Debug the WT results file and write to the output file
    wt_results_file_path = os.path.join(WT_RUN_DIR, result_file_name)
    debug_mmpbsa_file_to_file(wt_results_file_path, output_debug_file)
    
    # Debug the MT results file and append to the same output file
    mt_results_file_path = os.path.join(MT_RUN_DIR, result_file_name)
    debug_mmpbsa_file_to_file(mt_results_file_path, output_debug_file)

    print(f"\n--- Debugging complete. The output has been saved to: {output_debug_file} ---")