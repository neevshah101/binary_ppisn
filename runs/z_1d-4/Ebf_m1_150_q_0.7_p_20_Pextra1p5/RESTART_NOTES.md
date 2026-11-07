# Restart notes

- 2026-11-06 (1st restart, 5%): rerun in place (no separate copy) from photo `4000` (`./re 4000`).
- Change from the original run: wind mass-loss rate enhanced by 5% when surface T < 1e4 K
  (`x_ctrl(3) = 0.05d0`, `x_ctrl(4) = 1d4` in `inlist_both`; applied in `VMS_framework` in
  `src/run_star_extras.f90`). Same routine as `../m1_150_q_0.7_p_20_coolwind1p02` (which uses 2%).
- Original run (no enhancement) ended at `min_timestep_limit` (SLURM job 24545717, log `slurm-24545717.out`).
- 5% run: SLURM job 24559742; reached model 4206, then `min_timestep_limit`.
- 2026-11-06 (2nd restart, 10%): `x_ctrl(3)` raised to `0.10d0`, restarted again from photo `4000` (`./re 4000`),
  SLURM job 24559942. Currently active setting: 10%.
- MESA appends to LOGS*/history.data and binary_history.data on restart: rows after the restart point with
  repeated model numbers (>= 4000) come from the latest restart and supersede the original ones.
