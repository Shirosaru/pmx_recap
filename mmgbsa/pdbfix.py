from pdbfixer import PDBFixer
from openmm.app import PDBFile, Modeller
from openmm.app.element import hydrogen
from openmm.app import Topology

input_pdb = "WTModel_1.pdb"
output_pdb = "no_hydrogens.pdb"

# Load and fix the structure
fixer = PDBFixer(filename=input_pdb)
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.removeHeterogens(keepWater=False)

# DON'T add hydrogens
# fixer.addMissingHydrogens(pH=7.0)  ← commented out

# Remove existing hydrogens using Modeller
modeller = Modeller(fixer.topology, fixer.positions)
modeller.delete([atom for atom in modeller.topology.atoms() if atom.element == hydrogen])

# Save the fixed PDB without hydrogens
with open(output_pdb, 'w') as out_file:
    PDBFile.writeFile(modeller.topology, modeller.positions, out_file)

print(f"✅ Hydrogen-free structure saved as: {output_pdb}")
