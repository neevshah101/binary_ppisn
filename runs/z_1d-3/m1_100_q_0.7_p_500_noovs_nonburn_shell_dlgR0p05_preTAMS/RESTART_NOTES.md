# Restart notes

- 2026-11-09: rerun in place (no separate copy) from photo `x500` (`./re x500` in `submit.slurm`), SLURM job 24591806.
- Change from the original run: `delta_lgR_limit = 0.01d0` added to `inlist_both` (default 0.5). `delta_lgR_hard_limit` left at default (-1, off).
- Original run ended at model 522 with `overflow from L2 (R_L2) surface` (SLURM job 24591105, `slurm-24591105.out`): R1 jumped 702 -> 1001 Rsun (dlogR = 0.154) in one step, past r_L2 ~ 960 Rsun.
- MESA appends to LOGS*/history.data and binary_history.data on restart: rows with model numbers >= 500 after the restart point supersede the original ones.
- Job 24591806 reproduced the original run exactly (same jump at model 522, L2 termination): `delta_lgR_limit` is ignored unless
  `delta_lgR_limit_min_lgR` (default 1d99) is lowered (`$MESA_DIR/star/private/timestep.f90:958`).
- 2026-11-09 (2nd restart): added `delta_lgR_limit_min_lgR = 0d0` to `inlist_both`, restarted again from `x500`, SLURM job 24591848. Currently active settings: `delta_lgR_limit = 0.01d0`, `delta_lgR_limit_min_lgR = 0d0`.
