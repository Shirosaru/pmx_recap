import mdtraj as md
from openmm.app import *
from openmm import *
from openmm.unit import *

# 1. Setup - Load your AMBER files
# Replace 'complex.prmtop' with your antibody topology
prmtop = AmberPrmtopFile('antibody_complex.prmtop')
inpcrd = AmberInpcrdFile('antibody_complex.inpcrd')

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                             constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picosecond)
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)

# 2. Minimize and Equilibrate
print("Minimizing...")
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(300*kelvin)

# 3. Production Run (e.g., 1ns)
# Saving to .dcd for conversion later
simulation.reporters.append(DCDReporter('trajectory.dcd', 5000))
simulation.reporters.append(StateDataReporter(sys.stdout, 5000, step=True, potentialEnergy=True, temperature=True))

print("Running Production...")
simulation.step(500000) # 1ns simulation

# 4. Convert DCD to NetCDF (Required for MMPBSA.py)
print("Converting trajectory for MMPBSA...")
traj = md.load('trajectory.dcd', top='antibody_complex.prmtop')
traj.save_netcdf('trajectory.nc')