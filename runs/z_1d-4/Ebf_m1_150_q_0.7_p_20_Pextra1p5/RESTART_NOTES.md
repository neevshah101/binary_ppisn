# Restart notes

- 2026-11-06: rerun in place (no separate copy) from photo `4000` (`./re 4000`).
- Change from the original run: wind mass-loss rate enhanced by 5% when surface T < 1e4 K
  (`x_ctrl(3) = 0.05d0`, `x_ctrl(4) = 1d4` in `inlist_both`; applied in `VMS_framework` in
  `src/run_star_extras.f90`). Same routine as `../m1_150_q_0.7_p_20_coolwind1p02` (which uses 2%).
- Original run (no enhancement) ended at `min_timestep_limit` (SLURM job 24545717). Its LOGS1, LOGS2,
  binary_history.data and slurm log are saved in `backup_pre_coolwind1p05/` (not in git).
- Submitted as SLURM job 24559742.
- Everything in LOGS*/ and binary_history.data after model 4000 is from the 5%-enhanced run.
