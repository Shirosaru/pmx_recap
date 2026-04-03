# PMX Recap

Validation framework for chain A and chain B interaction metrics using MMGBSA and OpenMM approaches.

## Project Structure

### Core Directories

- **`mmgbsa/`** - MMGBSA (Molecular Mechanics Poisson Boltzmann Surface Area) validation
  - `run_mmgbsa_compare.py` - Main script to run and compare MMGBSA calculations
  - `mdp/` - GROMACS MDP parameter files for simulations
  - `Mut/` & `WT/` - Directories for mutant and wild-type structures

- **`openmm/`** - OpenMM validation approach
  - `OpenMM.py` - Main OpenMM simulation and analysis script
  - `compare_binding.py` - Script to compare binding metrics
  - `antibody_pipeline.py` - Antibody-specific pipeline
  - `setup_env.sh` - Environment setup script
  - Input PDB files for WT and mutant models

- **`input/`** - Input data and structures

### Archive

- **`archive/dev/`** - Archived development and experimental scripts
- **`archive/mmgbsa_legacy/`** - Archived legacy MMGBSA variants and analysis scripts

## Usage

Run MMGBSA validation:
```bash
cd mmgbsa/
python run_mmgbsa_compare.py
```

Run OpenMM validation:
```bash
cd openmm/
python OpenMM.py
```

Compare job security cultures (US at-will vs Japan lifetime employment):
```bash
python compare_job_security.py
```

## License

See LICENSE file for details.
