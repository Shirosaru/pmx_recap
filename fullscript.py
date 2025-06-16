import pmx
from pmx.utils import create_folder
from pmx import gmx, jobscript
import sys
import numpy as np
from IPython.core.display import clear_output
import os,shutil
import re
import subprocess
import glob
import random
from AAtutorial import *

# initialize the free energy environment object: it will store the main parameters for the calculations
fe = AAtutorial( )

# set the workpath
fe.workPath = 'workpath'
# set the path to the pdb file
fe.pdbfile = 'inputs/protein.pdb'
# set the path to the molecular dynamics parameter files
fe.mdpPath = 'inputs/mdppath'
# set the number of replicas (several repetitions of calculation are useful to obtain reliable statistics)
fe.replicas = 1
# choose the forcefield of your choice
fe.ff = 'amber99sb-star-ildn-mut'
# provide amino acids as single-letter code and the resid
fe.edges = [ ['W', 'F'] ]
fe.resids = [6]

# finally, let's prepare the overall free energy calculation directory structure
fe.prepareFreeEnergyDir( )

fe.hybrid_structure_topology(bVerbose=False)

fe.boxWaterIons( )

fe.prepare_simulation( simType='em' )
fe.run_simulation_locally( simType='em')

# set several parameters
fe.JOBsimcpu = 16
fe.JOBqueue = 'SLURM'
fe.JOBsource = ['/etc/profile.d/modules.sh','/usr/local/gromacs/GMXRC2024']
fe.JOBmodules = ['shared','owl/intel-mpi-default','cuda91']
fe.JOBbGPU = True
fe.JOBgmx = 'gmx_threads_AVX_256 mdrun -pin on'
fe.JOBpartition= 'p20,p24,p00,p32,p08,p16' #'short'

# create the jobscripts
fe.prepare_jobscripts(simType='em')

fe.prepare_simulation( simType='eq')
fe.prepare_jobscripts(simType='eq')

fe.workPath = 'workpath_precalculated'
fe.prepare_transitions()

fe.JOBsimcpu = 10
fe.JOBntomp = 2
fe.JOBnotr = 100 # the number of transitions
fe.JOBqueue = 'SLURM'
fe.JOBsource = ['/etc/profile.d/modules.sh','/usr/local/gromacs/GMXRC2024']
fe.JOBmodules = ['shared','owl/intel-mpi-default','cuda91']
fe.JOBbGPU = True
fe.JOBgmx = 'srun gmx_mpi_AVX_256 mdrun -pin on'
fe.JOBpartition= 'p20,p24,p00,p32,p08,p16' #'short'

fe.prepare_jobscripts(simType='transitions')

fe.run_analysis( bVerbose=True)