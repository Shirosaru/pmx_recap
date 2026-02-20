import os
import subprocess
import glob
import re
import shutil
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser
from pmx import gmx, Model, mutate, Topology, gen_hybrid_top
from pmx.utils import create_folder
import pmx.jobscript



class AAtutorial:
    """
    Setup and run free-energy PMX/GROMACS workflows.
    """

    def __init__(self, **kwargs):
        # Initialize GROMACS API
        gmx.set_gmxlib()

        # Results containers
        self.resultsAll = pd.DataFrame()
        self.resultsSummary = pd.DataFrame()

        # Paths and simulation parameters
        self.workPath = './'
        self.pdbfile = 'inputs/protein.pdb'
        self.original_pdbfile = self.pdbfile
        self.mdpPath = os.path.join(self.workPath, 'mdp')

        self.edges = {}
        self.resids = {}
        self.chains = []

        self.replicas = 5
        self.simTypes = ['em', 'eq', 'transitions']
        self.states = ['stateA', 'stateB']

        self.ff = 'amber99sb-star-ildn-mut'
        self.boxshape = 'dodecahedron'
        self.boxd = 1.0
        self.water = 'tip3p'
        self.conc = 0.15
        self.pname = 'NaJ'
        self.nname = 'ClJ'

        self.JOBqueue = 'SGE'
        self.JOBsimtime = 24
        self.JOBnotr = 100
        self.JOBsimcpu = 8
        self.JOBntomp = 1
        self.JOBbGPU = True
        self.JOBmodules = []
        self.JOBsource = []
        self.JOBexport = []
        self.JOBgmx = 'gmx mdrun'
        self.JOBpartition = ''

        for k, v in kwargs.items():
            setattr(self, k, v)

    def _read_path(self, path):
        return os.path.abspath(path)

    def _clean_backup_files(self, path):
        for f in glob.glob(os.path.join(path, '*#')):
            os.remove(f)

    def _get_specific_path(self, edge=None, state=None, r=None, sim=None):
        path = self.workPath
        if edge:
            path = os.path.join(path, edge)
        if state:
            path = os.path.join(path, state)
        if r:
            path = os.path.join(path, f'run{r}')
        if sim:
            path = os.path.join(path, sim)
        return path

    def _be_verbose(self, process, bVerbose=False):
        out, err = process.communicate()
        if bVerbose and out:
            print(out.decode())
        if err:
            print(err.decode())

    def grompp(self, f, c, p, o, r=None, po=None, maxwarn=None, other_flags=None):
        cmd = ['gmx', 'grompp', '-f', f, '-c', c, '-p', p, '-o', o]
        if r:
            cmd += ['-r', r]
        if po:
            cmd += ['-po', po]
        if maxwarn is not None:
            cmd += ['-maxwarn', str(maxwarn)]
        if other_flags:
            extra = other_flags.split() if isinstance(other_flags, str) else other_flags
            cmd += extra
        print("🌐 Running:", ' '.join(cmd))
        subprocess.run(cmd, check=True)

    def prepareFreeEnergyDir(self):
        # Load edges and setup folder hierarchy
        self._read_edges()
        self.mdpPath = self._read_path(self.mdpPath)
        self.workPath = self._read_path(self.workPath)
        create_folder(self.workPath)
        self._create_folder_structure()
        self._print_summary()
        self._print_folder_structure()
        print("✅ Directory setup complete\n")

    def _read_edges(self):
        if isinstance(self.edges, str) and os.path.isfile(self.edges):
            self._read_edges_from_file()
        else:
            self.edges = {
                f"{e[0]}{self.resids[i]}{e[1]}": e
                for i, e in enumerate(self.edges)
            }

    def _read_edges_from_file(self):
        # Implement file parsing to populate self.edges
        pass

    def _create_folder_structure(self):
        for edge in self.edges:
            for state in self.states:
                for r in range(1, self.replicas + 1):
                    for sim in self.simTypes:
                        p = self._get_specific_path(edge, state, r, sim)
                        create_folder(p)

    def _print_summary(self):
        print("📋 Setup Summary")
        print(f"   Work path: {self.workPath}")
        print(f"   PDB file : {self.pdbfile}")
        print(f"   MDP folder: {self.mdpPath}")
        print(f"   Replicas : {self.replicas}")
        print("   Mutations:")
        for e in self.edges:
            print(f"     - {e}")

    def _print_folder_structure(self):
        print("\n📂 Directory Layout:")
        print(f"{self.workPath}/")
        print("└── <edge>/")
        print("    └── stateA/stateB/")
        print("        └── run1..N/")
        print("            └── em/eq/transitions/\n")

    def hybrid_structure_topology(self, edges=None, bVerbose=False):
        print("🔧 Building hybrid topology for mutations...")

        if edges is None:
            edges = self.edges

        # Parse original PDB
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("orig", self.original_pdbfile)
        pdb_map = {
            (chain.id, resid.id[1]): resid.get_resname().upper()
            for model in struct
            for chain in model
            for resid in chain
            if resid.id[0] == ' '
        }

        # PMX model
        pmx_model = Model(self.original_pdbfile, rename_atoms=True)
        pmx_map = {
            (r.chain.id, r.id): r.atoms[0].resnr
            for r in pmx_model.residues
            if r.atoms
        }

        # Preprocess original PDB with pdb2gmx to add missing hydrogens
        print("🔄 Running pdb2gmx on original PDB to add hydrogens...")
        preprocessed_pdb = os.path.join(self.workPath, "processed_input.pdb")
        gmx.pdb2gmx(
            f=self.original_pdbfile,
            o=preprocessed_pdb,
            p=os.path.join(self.workPath, "dummy.top"),
            ff=self.ff,
            water=self.water,
            other_flags='-ignh -missing'
        )

        for i, edge in enumerate(edges):
            if isinstance(edge, str):
                m = re.match(r'^([A-Z]{3})(\d+)([A-Z]{3})$', edge.upper())
                if not m:
                    raise ValueError(f"Unrecognized mutation format: {edge}")
                old3, resnum_str, new3 = m.groups()
                resnum = int(resnum_str)
            else:
                old3, new3 = edge
                resnum = pmx_map[(self.chains[i], self.resids[i])]

            chain_id = self.chains[i]
            pmx_id = self.resids[i]
            pdb_name = pdb_map.get((chain_id, resnum), 'UNK')

            print(f"\n➡️  {old3}{resnum}{new3} → PMX chain {chain_id} ID {pmx_id} ({pdb_name}{resnum})")

            outdir = self._get_specific_path(edge=edge)

            # Mutate using the preprocessed structure (with hydrogens)
            mutated = mutate(
                m=Model(preprocessed_pdb, rename_atoms=True),
                mut_resid=pmx_id,
                mut_resname=new3,
                ff=self.ff
            )
            mutant_pdb = os.path.join(outdir, "mutant.pdb")
            mutated.write(mutant_pdb)

            # pdb2gmx on the mutant to generate topology
            gmx.pdb2gmx(
                f=mutant_pdb,
                o=os.path.join(outdir, "conf.pdb"),
                p=os.path.join(outdir, "topol_prev.top"),
                ff=self.ff,
                water=self.water,
                other_flags='-ignh -missing'
            )

            # Hybrid topology generation
            top = Topology(os.path.join(outdir, "topol_prev.top"), ff=self.ff)
            hyb, _ = gen_hybrid_top(top)
            hyb.write(os.path.join(outdir, "topol.top"), scale_mass=0.33)

        print("🎉 Hybrid structures ready\n")

    def boxWaterIons(self, edges=None, bBoxProt=True, bWatProt=True, bIonProt=True):
        import shutil

        print('=============================')
        print('🔧 Step: Box, Solvate, Add Ions')
        print('=============================')

        if edges is None:
            edges = self.edges

        for edge in edges:
            print(f"\n➡️ Processing edge: {edge}")
            out = self._get_specific_path(edge=edge)
            conf = os.path.join(out, 'conf.pdb')
            topol = os.path.join(out, 'topol.top')

            # === Ensure hybrid .itp files are present ===
            for required_itp in ['stateA.itp', 'stateB.itp']:
                itp_path = os.path.join(out, required_itp)
                if not os.path.exists(itp_path):
                    # Try copying from hybrid source path
                    src = os.path.join(out, required_itp)  # Same as destination folder
                    if not os.path.exists(itp_path):
                        src = os.path.join(out, required_itp)
                        if os.path.exists(src) and src != itp_path:
                            print(f"📄 Copying {required_itp} from {src}")
                            shutil.copy(src, itp_path)
                        elif not os.path.exists(src):
                            raise FileNotFoundError(f"❌ {required_itp} not found in {out}")
                    if os.path.exists(src):
                        print(f"📄 Copying {required_itp} from hybrid path.")
                        shutil.copy(src, itp_path)
                    else:
                        raise FileNotFoundError(f"❌ {required_itp} not found in {out} or {self.hybridPath}")

            # === Patch topol.top if needed ===
            with open(topol, 'r') as f:
                topol_lines = f.readlines()

            if not any('stateA.itp' in line for line in topol_lines):
                insert_idx = next((i for i, line in enumerate(topol_lines)
                                if line.strip().startswith('[ moleculetype ]')), len(topol_lines))

                topol_lines.insert(insert_idx, '#include "stateA.itp"\n')
                topol_lines.insert(insert_idx + 1, '#include "stateB.itp"\n')

                with open(topol, 'w') as f:
                    f.writelines(topol_lines)

                print("🔧 Patched topol.top to include stateA.itp and stateB.itp.")

            # === Box ===
            if bBoxProt:
                box = os.path.join(out, 'box.pdb')
                print("📦 Running editconf to create box")
                try:
                    subprocess.run(['gmx', 'editconf', '-f', conf, '-o', box, '-bt', self.boxshape, '-d', str(self.boxd)],
                                check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as e:
                    print("⚠️ editconf failed, retrying with explicit -box:", e.stderr)
                    subprocess.run(['gmx', 'editconf', '-f', conf, '-o', box, '-bt', self.boxshape,
                                    '-box', str(self.boxd), str(self.boxd), str(self.boxd)],
                                check=True)

            # === Solvate ===
            if bWatProt:
                box = os.path.join(out, 'box.pdb')
                water = os.path.join(out, 'water.pdb')
                print("💧 Running solvate")
                subprocess.run(['gmx', 'solvate', '-cp', box, '-cs', 'spc216.gro',
                                '-p', topol, '-o', water],
                            check=True, capture_output=True)

            # === GROMPP for Ion ===
            if bIonProt:
                water = os.path.join(out, 'water.pdb')
                tpr = os.path.join(out, 'tpr.tpr')
                mdout = os.path.join(out, 'mdout.mdp')

                print("⚙️ Running grompp for ion setup")
                result = subprocess.run(['gmx', 'grompp',
                                        '-f', os.path.join(self.mdpPath, 'em_l0.mdp'),
                                        '-c', water,
                                        '-p', topol,
                                        '-o', tpr,
                                        '-r', water,
                                        '-po', mdout,
                                        '-maxwarn', '4'],
                                        capture_output=True, text=True)

                if "No default" in result.stderr:
                    print("❌ Topology has missing bonded parameters:\n", result.stderr)
                    raise RuntimeError("Hybrid topology has missing bonded parameters. Fix topology.")

                print("✅ grompp passed. Proceeding to genion")
                subprocess.run(['gmx', 'genion', '-s', tpr, '-p', topol,
                                '-o', os.path.join(out, 'ions.pdb'),
                                '-conc', str(self.conc), '-neutral',
                                '-pname', self.pname, '-nname', self.nname],
                            input=b'WE\n', check=True)

            # === Clean backup files ===
            self._clean_backup_files(out)

        print('\n✅ Finished: boxWaterIons\n')


    def _prepare_single_tpr(self, simpath, toppath, state, simType, empath=None, frameNum=0):
        mdp_map = {'em': 'em', 'eq': 'eq', 'transitions': 'ti'}
        prefix = mdp_map.get(simType, '')
        mdp_file = os.path.join(
            self.mdpPath,
            f"{prefix}_l0.mdp" if state == 'stateA' else f"{prefix}_l1.mdp"
        )

        in_str = {
            'em': os.path.join(toppath, 'ions.pdb'),
            'eq': os.path.join(empath, 'confout.gro') if empath else None,
            'transitions': os.path.join(simpath, f"frame{frameNum}.gro")
        }[simType]

        out_tpr = os.path.join(simpath, "tpr.tpr")
        out_mdout = os.path.join(simpath, "mdout.mdp")
        print(f"📂 Preparing TPR: {simpath}")
        self.grompp(
            f=mdp_file,
            c=in_str,
            p=os.path.join(toppath, "topol.top"),
            o=out_tpr,
            r=in_str,
            po=out_mdout,
            maxwarn=4
        )
        self._clean_backup_files(simpath)

    def prepare_simulation(self, edges=None, simType='em'):
        if edges is None:
            edges = self.edges
        print(f"🧪 Preparing '{simType}' runs")
        for edge in edges:
            topo = self._get_specific_path(edge=edge)
            for state in self.states:
                for r in range(1, self.replicas + 1):
                    sim_dir = self._get_specific_path(edge, state, r, simType)
                    em_dir = self._get_specific_path(edge, state, r, 'em')
                    self._prepare_single_tpr(sim_dir, topo, state, simType, empathem=em_dir)
        print(f"✅ {simType} preparation done\n")

    def _run_mdrun(self, tpr, ener, confout, mdlog, trr, dhdl=None, xtc=None, bVerbose=False):
        cmd = ['gmx', 'mdrun', '-s', tpr, '-e', ener, '-c', confout, '-g', mdlog]
        if xtc:
            cmd += ['-dhdl', dhdl, '-x', xtc, '-o', trr]
        else:
            cmd += ['-o', trr]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self._be_verbose(proc, bVerbose=bVerbose)
        proc.wait()

    def run_simulation_locally(self, edges=None, simType='em', bVerbose=False):
        if edges is None:
            edges = self.edges
        print(f"⚙️ Running '{simType}' simulations")
        for edge in edges:
            for state in self.states:
                for r in range(1, self.replicas + 1):
                    sim_dir = self._get_specific_path(edge, state, r, simType)
                    tpr = os.path.join(sim_dir, "tpr.tpr")
                    ener = os.path.join(sim_dir, "ener.edr")
                    confout = os.path.join(sim_dir, "confout.gro")
                    mdlog = os.path.join(sim_dir, "md.log")
                    trr = os.path.join(sim_dir, "traj.trr")
                    xtc = os.path.join(sim_dir, "traj.xtc") if simType != 'em' else None
                    dhdl = os.path.join(sim_dir, "dhdl.xvg") if simType != 'em' else None

                    self._run_mdrun(tpr, ener, confout, mdlog, trr, dhdl=dhdl, xtc=xtc, bVerbose=bVerbose)
                    self._clean_backup_files(sim_dir)

        print(f"✅ Local runs for '{simType}' completed\n")

    def prepare_transitions(self, edges=None, bGenTpr=True):
        if edges is None:
            edges = self.edges
        print("🌀 Generating transition snapshots")
        for edge in edges:
            for state in self.states:
                for r in range(1, self.replicas + 1):
                    eq_dir = self._get_specific_path(edge, state, r, 'eq')
                    ti_dir = self._get_specific_path(edge, state, r, 'transitions')
                    create_folder(ti_dir)

                    # fetch frame0 and rename
                    subprocess.run([
                        'gmx', 'trjconv',
                        '-s', os.path.join(eq_dir, "tpr.tpr"),
                        '-f', os.path.join(eq_dir, "traj.trr"),
                        '-o', os.path.join(ti_dir, "frame.gro"),
                        '-sep', '-ur', 'compact', '-pbc', 'mol', '-b', '2000'
                    ], input=b"0\n", check=True)
                    os.rename(os.path.join(ti_dir, "frame0.gro"), os.path.join(ti_dir, "frame100.gro"))
                    self._clean_backup_files(ti_dir)

                    if bGenTpr:
                        for i in range(1, 101):
                            self._prepare_single_tpr(
                                simpath=ti_dir,
                                toppath=self._get_specific_path(edge),
                                state=state,
                                simType='transitions',
                                empath=eq_dir,
                                frameNum=i
                            )
        print("✅ Transitions prepared\n")

    def run_analysis(self, edges=None, bProt=True, bParseOnly=False, bVerbose=False):
        if edges is None:
            edges = self.edges
        print("📊 Running analysis across replicas")
        for edge in edges:
            if bProt:
                for r in range(1, self.replicas + 1):
                    analysis_dir = os.path.join(self._get_specific_path(edge), f"analyse{r}")
                    create_folder(analysis_dir)
                    aA = os.path.join(self._get_specific_path(edge, 'stateA', r, 'transitions'))
                    aB = os.path.join(self._get_specific_path(edge, 'stateB', r, 'transitions'))
                    cmd = (
                        f"pmx analyse -fA {aA}/*xvg -fB {aB}/*xvg "
                        f"-o {analysis_dir}/results.txt "
                        f"-oA {analysis_dir}/integ0.dat -oB {analysis_dir}/integ1.dat "
                        f"-w {analysis_dir}/wplot.png -t 298 -b 100"
                    )
                    subprocess.run(cmd, shell=True, check=True)
                    if bVerbose:
                        with open(os.path.join(analysis_dir, "results.txt")) as f:
                            print("".join([l for l in f if "ANALYSIS" in l or "INTEG" in l]))

        print("✅ Analysis completed\n")

    def prepare_jobscripts(self, edges=None, simType='em'):
        if edges is None:
            edges = self.edges
        print(f"📁 Generating job scripts for: {simType}")
        jobdir = os.path.join(self.workPath, f"{simType}_jobscripts")
        create_folder(jobdir)
        counter = 0

        for edge in edges:
            for state in self.states:
                for r in range(1, self.replicas + 1):
                    sim_dir = self._get_specific_path(edge, state, r, simType)
                    base_script = os.path.join(jobdir, f"jobscript{counter}")
                    job = pmx.jobscript.Jobscript(
                        fname=base_script,
                        queue=self.JOBqueue,
                        simcpu=self.JOBsimcpu,
                        jobname=f"j{edge}_{state}_{r}_{simType}",
                        modules=self.JOBmodules,
                        source=self.JOBsource,
                        gmx=self.JOBgmx,
                        partition=self.JOBpartition,
                        bGPU=self.JOBbGPU
                    )

                    job.cmds = [f"cd {sim_dir}", "$GMXRUN -s tpr.tpr"]
                    if simType == 'transitions':
                        total = self.JOBnotr
                        step = self.JOBsimcpu
                        for idx_start in range(0, total, step):
                            idx_end = min(idx_start+step, total)
                            self._create_jobscript(
                                fname=os.path.join(jobdir, f"jobscript{counter}"),
                                jobname=f"{edge}_run{r}_chunk{idx_start}",
                                modules=self.JOBmodules,
                                gmx=self.JOBgmx,
                                partition=self.JOBpartition,
                                simpath=sim_dir,
                                indx=[idx_start, idx_end]
                            )
                            counter += 1
                    else:
                        job.create_jobscript()
                        counter += 1

        # Write submission wrapper
        submit_py = os.path.join(jobdir, "submit.py")
        with open(submit_py, "w") as fp:
            fp.write("import os\n")
            for i in range(counter):
                if self.JOBqueue == 'SGE':
                    cmd = f"qsub jobscript{i}"
                else:
                    cmd = f"sbatch jobscript{i}"
                fp.write(f"os.system('{cmd}')\n")

        print(f"✅ {counter} job scripts generated in {jobdir}\n")

    def _create_jobscript(self, fname, jobname, modules, gmx, partition, simpath, indx):
        with open(fname, "w") as fp:
            fp.write("#!/bin/bash\n")
            fp.write(f"#SBATCH --job-name={jobname}\n")
            fp.write("#SBATCH --nodes=1\n")
            fp.write(f"#SBATCH --ntasks-per-node={indx[1]-indx[0]}\n")
            fp.write(f"#SBATCH --cpus-per-task={self.JOBntomp}\n")
            fp.write("#SBATCH --mail-type=none\n")
            fp.write(f"#SBATCH --time={self.JOBsimtime}:00:00\n")
            fp.write(f"#SBATCH --partition={partition}\n")
            if self.JOBbGPU:
                fp.write("#SBATCH --gres=gpu:1\n")
            fp.write("\n")
            for src in self.JOBsource:
                fp.write(f"source {src}\n")
            for mod in modules:
                fp.write(f"module load {mod}\n")
            fp.write(f'export GMXRUN="{gmx}"\n')
            fp.write("cd $TMPDIR\n")
            fp.write(f"for i in {{{indx[0]+1}..{indx[1]}}}; do\n")
            fp.write("  mkdir ti$i; cp {}/ti$i.tpr ti$i/tpr.tpr\n".format(simpath))
            fp.write("done\n")
            fp.write("$GMXRUN -s tpr.tpr -multidir ti{{$indx[0]+1..$indx[1]}}\n")
            fp.write("for i in {{{0}+1..{1}}}; do cp ti$i/dhdl.xvg {2}/dhdl$i.xvg; done\n".format(indx[0], indx[1], simpath))

