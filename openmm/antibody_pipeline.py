import os
import sys
import shutil
import numpy as np
try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from pdbfixer import PDBFixer
    import parmed as pmd
    import mdtraj as md
except ImportError as e:
    print(f"Missing dependency: {e}")
    sys.exit(1)

class AntibodyPipeline:
    def __init__(self, pdb_input, output_dir="output", chain_receptor="A", chain_ligand="B"):
        self.pdb_input = pdb_input
        self.output_dir = output_dir
        self.chain_receptor = chain_receptor
        self.chain_ligand = chain_ligand
        
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.filename_base = os.path.splitext(os.path.basename(pdb_input))[0]
        self.trajectory_dcd = os.path.join(self.output_dir, f"{self.filename_base}_traj.dcd")
        self.topology_pdb = os.path.join(self.output_dir, f"{self.filename_base}_solvated.pdb")
        self.dry_pdb = os.path.join(self.output_dir, f"{self.filename_base}_dry.pdb") # Complex for analysis

    def prepare_system(self):
        print(f"[{self.filename_base}] Loading and fixing PDB...")
        fixer = PDBFixer(filename=self.pdb_input)
        
        # Identify chains to keep
        desired_chains = [self.chain_receptor, self.chain_ligand]
        chains_to_remove = [c.id for c in fixer.topology.chains() if c.id not in desired_chains]
        if chains_to_remove:
            print(f"Removing chains: {chains_to_remove}")
            fixer.removeChains(chains_to_remove)

        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        return fixer

    def run_simulation(self, steps=10000):
        fixer = self.prepare_system()
        
        # Save Dry Structure for Analysis later (Reference Topology)
        with open(self.dry_pdb, 'w') as f:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

        print(f"[{self.filename_base}] Solvating...")
        modeller = app.Modeller(fixer.topology, fixer.positions)
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        modeller.addSolvent(forcefield, padding=1.0*unit.nanometers, ionicStrength=0.15*unit.molar)

        # Create System
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, 
                                         nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds)
        
        # Platform
        try:
            platform = mm.Platform.getPlatformByName('CUDA')
        except:
            try:
                platform = mm.Platform.getPlatformByName('OpenCL')
            except:
                platform = mm.Platform.getPlatformByName('CPU')
        
        print(f"[{self.filename_base}] Using platform: {platform.getName()}")
        
        integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
        simulation = app.Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(modeller.positions)
        
        print(f"[{self.filename_base}] Minimizing...")
        simulation.minimizeEnergy()
        
        # Output PDB of solvated system
        with open(self.topology_pdb, 'w') as f:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

        # Run
        simulation.reporters.append(app.DCDReporter(self.trajectory_dcd, 100)) # Save every 100 steps
        simulation.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, speed=True))
        
        print(f"[{self.filename_base}] Simulating {steps} steps...")
        simulation.step(steps)
        print(f"[{self.filename_base}] Simulation Done.")

    def calculate_binding_energy_gbsa(self, start_frame=0, end_frame=None, stride=1):
        """
        Calculates binding energy using OBC2 GBSA model on the dry trajectory.
        E_binding = E_complex - E_receptor - E_ligand
        """
        print(f"[{self.filename_base}] Calculating Binding Energy (GBSA)...")
        
        # Load Trajectory (Solvated DCD)
        traj = md.load(self.trajectory_dcd, top=self.topology_pdb)
        
        # Remove solvent for GBSA calculation
        # Identify Receptor and Ligand atoms
        # topology is mdtraj topology
        # We need to map this back to OpenMM topology for Energy calculation
        
        # 1. Load Dry Topology (Complex)
        pdb_complex = app.PDBFile(self.dry_pdb)
        forcefield = app.ForceField('amber14-all.xml', 'implicit/obc2.xml')
        
        # We need to extract coordinates from the solvated trajectory corresponding to the dry atoms
        # MDTraj is good for this slicing
        dry_indices = traj.topology.select(f"chainid 0 or chainid 1") # Assuming 0 is A, 1 is B if we cleaned it right. 
        # Better: Select by Protein
        dry_indices = traj.topology.select("protein")
        
        dry_traj = traj.atom_slice(dry_indices)
        
        # Sliced Trajectory for analysis
        if end_frame is None: end_frame = len(dry_traj)
        analysis_frames = dry_traj[start_frame:end_frame:stride]
        
        # Create Systems
        # System Complex
        system_complex = forcefield.createSystem(pdb_complex.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
        
        # System Receptor (Chain A / First Chain) - create by Modeller
        modeller_rec = app.Modeller(pdb_complex.topology, pdb_complex.positions)
        # Delete Chain B
        chains = list(modeller_rec.topology.chains())
        # Assuming Chain A is index 0, Chain B is index 1
        to_delete = [chains[1]] 
        modeller_rec.delete(to_delete)
        system_receptor = forcefield.createSystem(modeller_rec.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
        
        # System Ligand (Chain B / Second Chain)
        modeller_lig = app.Modeller(pdb_complex.topology, pdb_complex.positions)
        chains = list(modeller_lig.topology.chains())
        to_delete = [chains[0]]
        modeller_lig.delete(to_delete)
        system_ligand = forcefield.createSystem(modeller_lig.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
        
        # Helper to get energy
        def get_potential_energy(syst, positions):
            integrator = mm.VerletIntegrator(0.001)
            # Use CPU for analysis to avoid context thrashing or small system overhead on GPU
            platform = mm.Platform.getPlatformByName('CPU') 
            context = mm.Context(syst, integrator, platform)
            context.setPositions(positions)
            state = context.getState(getEnergy=True)
            return state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)

        energies = []
        
        # Map indices for Receptor and Ligand from the Complex Dry topology
        # Re-using mdtraj topology of the dry structure to map atoms might be complex if atom order changes
        # But Modeller usually preserves order of remaining atoms.
        # We will assume:
        # Complex Pos = Frame Pos
        # Receptor Pos = First N atoms of Frame Pos (if Receptor is first chain)
        # Ligand Pos = Last M atoms of Frame Pos
        
        n_atoms_rec = list(modeller_rec.topology.atoms())[-1].index + 1
        # Quick verify: Reconstruct indices
        # Actually safer way: mapping. But if we stick to Chain A then B order:
        
        for k in range(len(analysis_frames)):
            frame_pos = analysis_frames.xyz[k] * 10.0 # nm to Angstrom? No, OpenMM uses nm? 
            # MDTraj xyz is in nanometers. OpenMM expects nm.
            # Convert numpy array to Quantity or list of Vec3?
            # context.setPositions accepts (N,3) array
            
            # Complex Energy
            e_complex = get_potential_energy(system_complex, frame_pos)
            
            # Receptor Energy (First N atoms)
            pos_rec = frame_pos[:n_atoms_rec]
            e_rec = get_potential_energy(system_receptor, pos_rec)
            
            # Ligand Energy (Rest)
            pos_lig = frame_pos[n_atoms_rec:]
            e_lig = get_potential_energy(system_ligand, pos_lig)
            
            ddG = e_complex - (e_rec + e_lig)
            energies.append(ddG)
            
        avg_binding = np.mean(energies)
        std_binding = np.std(energies)
        print(f"[{self.filename_base}] Average Binding Energy: {avg_binding:.2f} +/- {std_binding:.2f} kcal/mol")
        return avg_binding, std_binding
