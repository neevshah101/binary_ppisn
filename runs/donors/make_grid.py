#!/usr/bin/env python3
"""
Create the donor grid from binary_ppisn/binary_ppisn_progenitor (non-rotating, both stars):

    runs/donors/z_{Z}/m1_{M1}_q_{Q}_p_{P}/
        inlist_binary <- evolve_both_stars = .false., m1, m2 = q*m1, initial_period_in_days
        inlist_extra  <- new_Z, Zbase
        inlist1/2     <- new_rotation_flag = .false., new_omega_div_omega_crit = 0

Also writes runs/donors/grid_runs.txt (one run dir per line, the order used by
submit_grid.sh) and a RUN_LOG.md in each run. Existing run dirs are skipped.

    python make_grid.py
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent.resolve()
REPO = HERE.parents[1]                       # binary_ppisn/
sys.path.insert(0, str(REPO))
from setup_run import IGNORE, TEMPLATE, fmt, fmt_z, fortran, patch, set_rotation  # noqa: E402

Z_VALUES = [1e-3, 1e-4]
M1_VALUES = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0]
Q_VALUES = [0.7]
P_VALUES = [20.0, 100.0, 500.0]

RUNLOG = Path.home() / ".claude/skills/mesa-runlog/runlog.py"


def write_mesa_version(run_dir):
    """Record the branch/commit of the MESA in $MESA_DIR (the template's MESA_VERSION may be stale)."""
    mesa_dir = os.environ["MESA_DIR"]
    git = lambda *a: subprocess.run(["git", "-C", mesa_dir, *a], check=True,
                                    capture_output=True, text=True).stdout.strip()
    (run_dir / "MESA_VERSION").write_text(
        f"MESA branch: {git('branch', '--show-current')}\n"
        f"MESA commit: {git('rev-parse', 'HEAD')}\n"
        f"MESA_DIR:    {mesa_dir}\n")


def main():
    runs, created = [], 0
    for z in Z_VALUES:
        for m1 in M1_VALUES:
            for q in Q_VALUES:
                m2 = round(m1 * q, 6)
                for p in P_VALUES:
                    run_dir = HERE / f"z_{fmt_z(z)}" / f"m1_{fmt(m1)}_q_{fmt(q)}_p_{fmt(p)}"
                    runs.append(run_dir)
                    if run_dir.exists():
                        print(f"[SKIP]    {run_dir.relative_to(HERE)}")
                        continue
                    shutil.copytree(TEMPLATE, run_dir, ignore=IGNORE)
                    patch(run_dir / "inlist_binary", {
                        "evolve_both_stars": ".false.",   # donor-only grid: star 2 is a point mass
                        "m1": fortran(m1),
                        "m2": fortran(m2),
                        "initial_period_in_days": fortran(p),
                    })
                    patch(run_dir / "inlist_extra", {"new_Z": fmt_z(z), "Zbase": fmt_z(z)})
                    set_rotation(run_dir / "inlist1", False, 0)
                    set_rotation(run_dir / "inlist2", False, 0)
                    write_mesa_version(run_dir)
                    created += 1
                    print(f"[CREATED] {run_dir.relative_to(HERE)}  (m2={m2})")
                    if RUNLOG.exists():
                        subprocess.run([sys.executable, str(RUNLOG), "init", str(run_dir),
                                        "--parent", str(TEMPLATE), "-m",
                                        f"runs/donors/make_grid.py: m1={fmt(m1)}, q={fmt(q)} (m2={m2}), "
                                        f"P={fmt(p)} d, Z={fmt_z(z)}, non-rotating (w1=w2=0), "
                                        f"evolve_both_stars=.false. (donor only)."],
                                       check=True, stdout=subprocess.DEVNULL)

    (HERE / "grid_runs.txt").write_text("".join(f"{r}\n" for r in runs))
    print(f"\n{created} created, {len(runs) - created} skipped; {len(runs)} runs in grid_runs.txt "
          f"-> array indices 0-{len(runs) - 1}")


if __name__ == "__main__":
    main()
