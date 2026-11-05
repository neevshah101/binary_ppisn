#!/usr/bin/env python3
"""
Create one MESA binary PPISN run from the binary_ppisn_progenitor template,
write a SLURM script for it, and commit + push the new run to GitHub.

    runs/z_{Z}/m1_{M1}_q_{Q}_p_{P}[_Pextra{X}]/   (suffix only if PEXTRA != 1)
        inlist_binary   <- m1, m2, initial_period_in_days
        inlist_extra    <- new_Z, Zbase, Pextra_factor (both stars)
        submit.slurm    <- (gitignored)

Usage: edit USER CONFIG below, then
    python setup_run.py            # create, commit, push
    python setup_run.py --submit   # also sbatch the job
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

# ── USER CONFIG ────────────────────────────────────────────────────────────────

M1 = 150.0          # donor mass (Msun)
Q = 0.7             # mass ratio q = m2/m1
PERIOD = 20.0       # initial orbital period (days)
Z = 1e-4            # metallicity (new_Z and Zbase)
PEXTRA = 1.5        # Pextra_factor (MESA default 1; only used with atm_option='T_tau', the default)

# SLURM
ACCOUNT = "binboom"
PARTITION = "standard"
CPUS = 8
TIME = "8:00:00"

# ── PATHS ──────────────────────────────────────────────────────────────────────

REPO = Path(__file__).parent.resolve()
TEMPLATE = REPO / "binary_ppisn_progenitor"
RUNS_DIR = REPO / "runs"
IGNORE = shutil.ignore_patterns(
    "build", ".mesa_temp_cache", ".ipynb_checkpoints",
    "LOGS*", "photos", "png", "*.mod", "*.data", "star", "binary", ".restart",
)

# ── HELPERS ────────────────────────────────────────────────────────────────────

def fmt(v):
    """150.0 -> '150', 20.5 -> '20.5'."""
    return str(int(v)) if float(v).is_integer() else str(v)

def fmt_z(v):
    """1e-4 -> '1d-4', 2.5e-4 -> '2.5d-4'."""
    mant, exp = f"{v:e}".split("e")
    return f"{mant.rstrip('0').rstrip('.')}d{int(exp)}"

def fortran(v):
    """150 -> '150d0'."""
    return f"{fmt(v)}d0"

def patch(path, values):
    """Replace `key = <number>` at the start of a line; fail if a key is missing."""
    text = path.read_text()
    for key, val in values.items():
        text, n = re.subn(
            rf"^(\s*{re.escape(key)}\s*=\s*)[-+\d.eEdD]+",
            lambda m: m.group(1) + val, text, flags=re.M,
        )
        if n != 1:
            sys.exit(f"ERROR: expected 1 '{key}' in {path}, found {n}")
    path.write_text(text)

def set_pextra(path, value):
    """Insert Pextra_factor at the top of the &controls namelist."""
    text = path.read_text()
    if re.search(r"^\s*Pextra_factor\s*=", text, flags=re.M):
        sys.exit(f"ERROR: Pextra_factor already set in {path}")
    text, n = re.subn(r"^(&controls[^\n]*\n)", rf"\g<1>Pextra_factor = {value}\n",
                      text, flags=re.M)
    if n != 1:
        sys.exit(f"ERROR: expected 1 '&controls' in {path}, found {n}")
    path.write_text(text)

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

def git(*args, **kw):
    return subprocess.run(["git", "-C", str(REPO), *args], check=True, **kw)

def ask_commit_message():
    msg = ""
    while not msg:
        try:
            msg = input("Commit message (required): ").strip()
        except (EOFError, KeyboardInterrupt):
            sys.exit("\nAborted: no commit message. Run folder was created but not committed.")
    return msg

# ── MAIN ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--submit", action="store_true", help="sbatch the job after pushing")
    args = ap.parse_args()

    m2 = round(M1 * Q, 6)
    name = f"m1_{fmt(M1)}_q_{fmt(Q)}_p_{fmt(PERIOD)}"
    if PEXTRA != 1:
        name += f"_Pextra{fmt(PEXTRA).replace('.', 'p')}"
    run_dir = RUNS_DIR / f"z_{fmt_z(Z)}" / name
    if run_dir.exists():
        sys.exit(f"ERROR: {run_dir} already exists")

    shutil.copytree(TEMPLATE, run_dir, ignore=IGNORE)
    patch(run_dir / "inlist_binary", {
        "m1": fortran(M1),
        "m2": fortran(m2),
        "initial_period_in_days": fortran(PERIOD),
    })
    patch(run_dir / "inlist_extra", {"new_Z": fmt_z(Z), "Zbase": fmt_z(Z)})
    set_pextra(run_dir / "inlist_extra", fortran(PEXTRA))
    slurm = write_slurm(run_dir, f"z{fmt_z(Z)}_{name}")
    print(f"Created {run_dir.relative_to(REPO)}  (m1={M1}, q={Q}, m2={m2}, P={PERIOD} d, Z={fmt_z(Z)}, Pextra_factor={PEXTRA})")

    git("add", str(run_dir), str(Path(__file__).resolve()))
    git("status", "--short")
    git("commit", "-m", ask_commit_message())
    git("push")

    if args.submit:
        subprocess.run(["sbatch", str(slurm)], check=True)
    else:
        print(f"Submit with: sbatch {slurm}")

if __name__ == "__main__":
    main()
