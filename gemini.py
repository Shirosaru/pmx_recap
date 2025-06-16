import os
import subprocess
from pmx.model import Model
from pmx.alchemy import mutate
from pmx import gmx

# === SETTINGS ===
input_pdb = "/home2/pmx/test/I74271_clean_clean.pdb"
output_dir = "/home2/pmx/test/fep_mutation_output"
forcefield = "amber99sb-ildn"
water_model = "tip3p"
mutation_chain = "H"
mutation_resid = 50
wt_resname = "ARG"
mut_resname = "ALA"

# === CREATE OUTPUT DIRECTORY ===
os.makedirs(output_dir, exist_ok=True)

# === 1. LOAD AND MUTATE THE PDB FILE USING PMX ===
print("Step 1: Mutating residue with PMX...")

# Load the PDB file
m = Model(input_pdb, rename_atoms=True)

# Perform mutation
m2 = mutate(
    m=m,
    mut_resid=mutation_resid,
    mut_resname=mut_resname,
    ff=forcefield,
    mut_chain=mutation_chain,
    verbose=True
)

# Save the mutated PDB file
mutated_pdb = os.path.join(output_dir, "mutated.pdb")
m2.write(mutated_pdb)
print(f"✅ Mutated PDB saved to: {mutated_pdb}")

# === 2. PREPROCESS THE MUTATED PDB WITH GROMACS ===
print("Step 2: Preprocessing mutated PDB with GROMACS...")

# Run pdb2gmx to generate topology and processed PDB
gmx.pdb2gmx(
    f=mutated_pdb,
    o=os.path.join(output_dir, "conf.pdb"),
    p=os.path.join(output_dir, "topol.top"),
    ff=forcefield,
    water=water_model
)

# === 3. GENERATE HYBRID TOPOLOGY WITH PMX ===
print("Step 3: Generating hybrid topology with PMX...")

# Load the topology file
topol = Topology(os.path.join(output_dir, "topol.top"))

# Generate hybrid topology
pmxtop, _ = gen_hybrid_top(topol)

# Save the new topology
newtop_top = os.path.join(output_dir, "newtop.top")
pmxtop.write(newtop_top)
print(f"✅ Hybrid topology saved to: {newtop_top}")

# === 4. SYSTEM SETUP FOR SIMULATION ===
print("Step 4: Setting up the system for simulation...")

# Generate the system using GROMACS
gmx.editconf(
    f=os.path.join(output_dir, "conf.pdb"),
    o=os.path.join(output_dir, "conf.gro"),
    box="10.0"
)

gmx.solvate(
    f=os.path.join(output_dir, "conf.gro"),
    o=os.path.join(output_dir, "solv.gro"),
    p=os.path.join(output_dir, "topol.top")
)

gmx.grompp(
    f="mdp/minim.mdp",
    c=os.path.join(output_dir, "solv.gro"),
    p=os.path.join(output_dir, "topol.top"),
    o=os.path.join(output_dir, "em.tpr")
)

gmx.mdrun(
    s=os.path.join(output_dir, "em.tpr"),
    o=os.path.join(output_dir, "em.trr"),
    x=os.path.join(output_dir, "em.xtc"),
    g=os.path.join(output_dir, "em.log"),
    e=os.path.join(output_dir, "em.edr")
)

print("\n🎉 Mutation and system setup completed successfully.")