import os
import argparse
from antibody_pipeline import AntibodyPipeline

def main():
    parser = argparse.ArgumentParser(description="Compare Binding Energy of WT vs Mut Antibody-Antigen")
    parser.add_argument("--steps", type=int, default=10000, help="Simulation steps (default 10k)")
    args = parser.parse_args()

    files = ["WTModel_0.pdb", "MutModel_0.pdb"]
    results = {}

    for pdb in files:
        if not os.path.exists(pdb):
            print(f"Error: {pdb} not found!")
            continue
            
        pipeline = AntibodyPipeline(pdb, output_dir=f"output_{pdb.replace('.pdb','')}")
        
        # 1. Run Simulation
        if not os.path.exists(pipeline.trajectory_dcd):
            pipeline.run_simulation(steps=args.steps)
        else:
            print(f"Trajectory for {pdb} exists, skipping simulation.")
            
        # 2. Analyze Binding
        avg, std = pipeline.calculate_binding_energy_gbsa(start_frame=10) # Skip first 10 frames as equilibration
        results[pdb] = (avg, std)

    print("\n" + "="*40)
    print("       BINDING ENERGY COMPARISON")
    print("="*40)
    print(f"{'Structure':<15} | {'Binding Energy (kcal/mol)':<25}")
    print("-" * 42)
    
    for pdb, (u, s) in results.items():
        print(f"{pdb:<15} | {u:6.2f} +/- {s:5.2f}")
        
    if len(results) == 2:
        wt_e = results["WTModel_0.pdb"][0]
        mut_e = results["MutModel_0.pdb"][0]
        diff = mut_e - wt_e
        print("-" * 42)
        print(f"Delta (Mut - WT): {diff:.2f} kcal/mol")
        if diff < 0:
            print("Conclusion: MUTANT binds STRONGER (more negative energy).")
        else:
            print("Conclusion: WT binds STRONGER (or Mutant binds weaker).")
    print("="*40)

if __name__ == "__main__":
    main()
