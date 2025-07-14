import os, subprocess, shutil, datetime

def run_cmd(cmd, cwd=None, input_text=None, label=""):
    logdir = os.path.join(cwd or ".", "logs")
    os.makedirs(logdir, exist_ok=True)
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    logfile = os.path.join(logdir, f"{label.replace(' ','_')}_{ts}.log")
    print(f"[RUN] {label} > {' '.join(cmd)}")
    proc = subprocess.run(cmd, cwd=cwd, input=input_text, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with open(logfile, "w") as f:
        f.write(f"CMD: {' '.join(cmd)}\n\nSTDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}")
    if proc.returncode:
        print(f"[ERROR] {label} failed. See {logfile}")
        raise RuntimeError(label)
    print(f"[INFO] {label} done.")
    return proc.stdout

def clean_env(wd, remove_top=False):
    for fn in os.listdir(wd):
        if fn.startswith("#") and fn.endswith("#"):
            os.remove(os.path.join(wd, fn))
        if remove_top and fn in ("topol.top", "solvated.gro"):
            os.remove(os.path.join(wd, fn))
    logs = os.path.join(wd, "logs")
    if os.path.isdir(logs): shutil.rmtree(logs)

def prepare_system(name, pdb, wd, mdp_dir):
    os.makedirs(wd, exist_ok=True)
    clean_env(wd, remove_top=True)
    for mdp in ("ions","minim","nvt","npt","md"):
        shutil.copy(os.path.join(mdp_dir, f"{mdp}.mdp"), os.path.join(wd, f"{mdp}.mdp"))

    run_cmd(["gmx","pdb2gmx","-f",pdb,"-o","processed.gro","-ff","amber99sb-ildn","-water","tip3p"], cwd=wd, input_text="6\n1\n", label="pdb2gmx")
    run_cmd(["gmx","editconf","-f","processed.gro","-o","boxed.gro","-c","-d","1.0","-bt","cubic"], cwd=wd, label="editconf")
    run_cmd(["gmx","solvate","-cp","boxed.gro","-cs","spc216.gro","-o","solvated.gro","-p","topol.top"], cwd=wd, label="solvate")
    run_cmd(["gmx","grompp","-f","ions.mdp","-c","solvated.gro","-p","topol.top","-o","ions.tpr"], cwd=wd, label="grompp ions")
    run_cmd(["gmx","genion","-s","ions.tpr","-o","ionized.gro","-p","topol.top","-neutral","-pname","NA","-nname","CL"], cwd=wd, input_text="SOL\n", label="genion")

    steps = [("minim","ionized.gro"), ("nvt","em.gro"), ("npt","nvt.gro"), ("md","npt.gro")]
    for step, inp in steps:
        run_cmd(["gmx","grompp","-f",f"{step}.mdp","-c",inp,"-p","topol.top","-o",f"{step}.tpr","-maxwarn","1"], cwd=wd, label=f"grompp {step}")
        run_cmd(["gmx","mdrun","-deffnm",step], cwd=wd, label=f"mdrun {step}")

    run_cmd(["gmx","trjconv","-s","md.tpr","-f","md.xtc","-o","md_noPBC.xtc","-pbc","mol","-center"], cwd=wd,
            input_text="1\n1\n", label="trjconv")

def parse_index_groups_detailed(index_file):
    group_data = {}
    current_group = None
    with open(index_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('[') and line.endswith(']'):
                current_group = line[1:-1].strip()
                group_data[current_group] = []
            elif current_group and line:
                group_data[current_group].extend(map(int, line.split()))
    return group_data

def fix_pdb_chains(infile, outfile, split_resid):
    with open(infile) as inp, open(outfile, 'w') as out:
        for line in inp:
            if line.startswith(('ATOM', 'HETATM')):
                resid = int(line[22:26])
                chain = 'A' if resid <= split_resid else 'B'
                newline = line[:21] + chain + line[22:]
                out.write(newline)
            else:
                out.write(line)

def log_dir_contents(wd, label):
    log_file = os.path.join(wd, f"debug_ls_{label}.txt")
    with open(log_file, "w") as f:
        for root, dirs, files in os.walk(wd):
            for file in sorted(files):
                path = os.path.join(root, file)
                size = os.path.getsize(path)
                f.write(f"{path} ({size} bytes)\n")

def run_mmgbsa(name, wd, pdb, chain1='A', chain2='B'):
    print(f"[STAGE] MM/GBSA for {name}")

    fixed_pdb = os.path.join(wd, "output_fixed.pdb")
    fix_pdb_chains(pdb, fixed_pdb, split_resid=300)

    input_text = (
        f"chain {chain1}\n"
        f"name 1 Chain1\n"
        f"chain {chain2}\n"
        f"name 2 Chain2\n"
        "q\n"
    )

    run_cmd(
        ["gmx", "make_ndx", "-f", "md.tpr", "-o", "index.ndx"],
        cwd=wd, input_text=input_text, label="make_ndx"
    )

    # Parse the index and create new renamed version
    groups = parse_index_groups_detailed(os.path.join(wd, "index.ndx"))
    if "Chain1" not in groups or "Chain2" not in groups:
        raise RuntimeError("Index groups missing Chain1/Chain2")

    complex_atoms = sorted(set(groups["Chain1"]) | set(groups["Chain2"]))
    output_ndx = os.path.join(wd, "index_renamed.ndx")
    with open(output_ndx, "w") as f:
        for g, atoms in groups.items():
            f.write(f"[ {g} ]\n")
            for i in range(0, len(atoms), 15):
                f.write(" ".join(str(a) for a in atoms[i:i+15]) + "\n")
            f.write("\n")
        f.write("[ Complex ]\n")
        for i in range(0, len(complex_atoms), 15):
            f.write(" ".join(str(a) for a in complex_atoms[i:i+15]) + "\n")

    # Dump it for inspection
    with open(output_ndx) as f:
        with open(os.path.join(wd, "logs", "index_renamed_dump.log"), "w") as logf:
            logf.write(f.read())

    idx = list(groups.keys())
    cg1, cg2 = idx.index("Chain1") + 1, idx.index("Chain2") + 1

    # Log contents before gmx_MMPBSA
    log_dir_contents(wd, "before_gmxMMPBSA")

    try:
        run_cmd([
            "gmx_MMPBSA", "-O",
            "-i", "../mmpbsa.in", "-cs", "md.tpr",
            "-ci", "index_renamed.ndx",
            "-cg", str(cg1), str(cg2),
            "-ct", "md_noPBC.xtc", "-cp", "topol.top"
        ], cwd=wd, label="gmx_MMPBSA")
    except Exception as e:
        for f in ["_GMXMMPBSA_COM.pdb", "_GMXMMPBSA_REC.pdb", "_GMXMMPBSA_LIG.pdb"]:
            fpath = os.path.join(wd, f)
            if os.path.exists(fpath):
                shutil.copy(fpath, os.path.join(wd, f"DEBUG_{f}"))
        log_dir_contents(wd, "after_gmxMMPBSA")
        raise

    out = os.path.join(wd, "FINAL_RESULTS_MMPBSA.dat")
    with open(out) as fh:
        for line in fh:
            if "Complex" in line:
                val = float(next(fh).split()[3])
                print(f"[SUCCESS] ΔG = {val:.3f} kcal/mol")
                return val
    raise RuntimeError("ΔG result not found")



def log_dir_contents(wd, label):
    log_file = os.path.join(wd, f"debug_ls_{label}.txt")
    with open(log_file, "w") as f:
        for root, dirs, files in os.walk(wd):
            for file in sorted(files):
                path = os.path.join(root, file)
                size = os.path.getsize(path)
                f.write(f"{path} ({size} bytes)\n")

def main():
    mdp_dir = "./mdp"
    systems = {
        "WT": ("/home2/mmGBSA/WT/WTModel_0.pdb", "WT_run"),
        "MT": ("/home2/mmGBSA/Mut/MutModel_0.pdb", "MT_run")
    }
    results = {}
    for name, (pdb, wd) in systems.items():
        print(f"\n=== {name} ===")
        try:
            prepare_system(name, pdb, wd, mdp_dir)
            results[name] = run_mmgbsa(name, wd, pdb, chain1='A', chain2='B')
        except Exception as e:
            results[name] = None
            print(f"[FAIL] {name}: {e}")

    with open("ddg_results.txt", "w") as out:
        out.write("System\tΔG\n")
        for k, v in results.items():
            out.write(f"{k}\t{v if isinstance(v, float) else 'ERROR'}\n")
        if all(isinstance(v, float) for v in results.values()):
            ddg = results["MT"] - results["WT"]
            out.write(f"\nΔΔG (MT-WT) = {ddg:.3f}\n")
            print(f"[FINAL] ΔΔG = {ddg:.3f}")
        else:
            print("[FINAL] Calculation failed")

if __name__ == "__main__":
    main()