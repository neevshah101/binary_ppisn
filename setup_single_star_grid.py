#!/usr/bin/env python3
"""
Create a grid of MESA single-star PPISN runs from the single_star_ppisn_progenitor
template (one run per MASSES x METALLICITIES pair), write a SLURM script for each,
and commit + push the new runs to GitHub.

    runs/single_star/z_{Z}/m_{M}[_Pextra{X}][_w_{W}]/
        (Pextra suffix only if PEXTRA != 1; w suffix only for rotating stars)
        inlist_extra    <- initial_mass, new_Z, Zbase, Pextra_factor
        inlist1         <- new_rotation_flag, new_omega_div_omega_crit
        submit.slurm    <- (gitignored)

Existing run folders are skipped, not overwritten.

Usage: edit USER CONFIG below, then
    python setup_single_star_grid.py            # create, commit, push
    python setup_single_star_grid.py --submit   # also sbatch every new run
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from setup_run import (IGNORE, ACCOUNT, PARTITION, CPUS, TIME,
                       fmt, fmt_z, fortran, patch, set_pextra, set_rotation,
                       git, ask_commit_message)

# ── USER CONFIG ────────────────────────────────────────────────────────────────

MASSES = [100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 275.0, 300.0]            # initial masses (Msun)
METALLICITIES = [1e-4, 1e-3]      # Z (new_Z and Zbase)
PEXTRA = 1                  # Pextra_factor (MESA default 1)
ROT = True                 # rotation on/off (all runs)
OMEGA = 0.2                 # initial omega/omega_crit (used only if ROT)

# ── PATHS ──────────────────────────────────────────────────────────────────────

REPO = Path(__file__).parent.resolve()
TEMPLATE = REPO / "single_star_ppisn_progenitor"
RUNS_DIR = REPO / "runs" / "single_star"

# ── HELPERS ────────────────────────────────────────────────────────────────────

def write_slurm(run_dir, job_name):
    script = run_dir / "submit.slurm"
    script.write_text(f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={ACCOUNT}
#SBATCH --partition={PARTITION}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={CPUS}
#SBATCH --time={TIME}
#SBATCH --output={run_dir}/slurm-%j.out

echo "=== Job $SLURM_JOB_ID on $(hostname), start: $(date) ==="
export OMP_NUM_THREADS=16
cd "{run_dir}" || exit 1
make clean && ./rn
echo "=== End: $(date) ==="
""")
    script.chmod(0o755)
    return script

def make_run(mass, z):
    name = f"m_{fmt(mass)}"
    if PEXTRA != 1:
        name += f"_Pextra{fmt(PEXTRA).replace('.', 'p')}"
    if ROT:
        name += f"_w_{fmt(OMEGA)}"
    run_dir = RUNS_DIR / f"z_{fmt_z(z)}" / name
    if run_dir.exists():
        print(f"SKIP {run_dir.relative_to(REPO)} (exists)")
        return None

    shutil.copytree(TEMPLATE, run_dir, ignore=IGNORE)
    patch(run_dir / "inlist_extra", {
        "initial_mass": fortran(mass),
        "new_Z": fmt_z(z),
        "Zbase": fmt_z(z),
    })
    set_pextra(run_dir / "inlist_extra", fortran(PEXTRA))
    set_rotation(run_dir / "inlist1", ROT, OMEGA)
    slurm = write_slurm(run_dir, f"ss_z{fmt_z(z)}_{name}")
    print(f"Created {run_dir.relative_to(REPO)}  (M={mass}, Z={fmt_z(z)}, "
          f"Pextra_factor={PEXTRA}, w={OMEGA if ROT else 'off'})")
    return run_dir, slurm

# ── MAIN ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--submit", action="store_true", help="sbatch every new run after pushing")
    args = ap.parse_args()

    created = [r for z in METALLICITIES for m in MASSES if (r := make_run(m, z))]
    if not created:
        sys.exit("Nothing new to create.")

    git("add", *(str(d) for d, _ in created), str(Path(__file__).resolve()))
    git("status", "--short")
    git("commit", "-m", ask_commit_message())
    git("push")

    for _, slurm in created:
        if args.submit:
            subprocess.run(["sbatch", str(slurm)], check=True)
        else:
            print(f"Submit with: sbatch {slurm}")

if __name__ == "__main__":
    main()
