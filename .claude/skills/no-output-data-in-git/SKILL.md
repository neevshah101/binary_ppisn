---
name: no-output-data-in-git
description: Use before any git add, commit, or push in the binary_ppisn repo. Keeps MESA/SLURM output data (LOGS, photos, png, .data, .mod, slurm logs, movies, PDFs, build artifacts) out of GitHub; only inputs and source code are committed.
---

# No output data in git

Only commit run *inputs* and *code*. Never stage, commit, or push output data.

## Never commit (mirrors `.gitignore`)

| Category | Patterns |
|---|---|
| MESA outputs / models | `*.mod`, `*.data`, `LOGS*/`, `LOGS*`, `photos/`, `png/`, `*.png` |
| SLURM | `*.slurm`, `*.sbatch`, `slurm*` (incl. `slurm-*.out`), `*.sh` (except `make_movies.sh`) |
| Text dumps | `*.txt` |
| Build artifacts | `build/`, `.mesa_temp_cache/`, `star`, `binary`, `.restart` |
| Misc | `.ipynb_checkpoints/`, `__pycache__/`, `*.mp4`, `*.pdf` |

Exception kept by `.gitignore`: `!make_movies.sh` (may be committed).

If `.gitignore` changes, it is the source of truth; re-read it before committing.

## What is fine to commit

`inlist*`, `src/*.f90`, `Makefile`, `re`, `rn`, `*columns.list`, `README.rst`, `MESA_VERSION`, `make_movies.sh`, Python scripts (`*.py`), notebooks only if asked.

## Procedure

1. `git status --short` and read `.gitignore`.
2. Stage explicit paths or directories and let `.gitignore` filter them. Never `git add -A` / `git add .` at the repo root blindly, and **never `git add -f`** on an ignored file unless the user names that file explicitly.
3. Before committing, run `git diff --cached --name-only` and confirm none of the staged paths match the patterns above (in particular no `.pdf`, `.png`, `.mp4`, `.data`, `.mod`, `LOGS*`, `photos/`, `slurm-*.out`). Unstage any that do with `git restore --staged <path>`.
4. Flag large files: anything staged > ~1 MB is almost certainly output; stop and ask.
5. Commit with a short message describing the run/code change, then push.
