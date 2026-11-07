# Restart notes

- Copied from `../m1_150_q_0.7_p_20` (including its photos and LOGS).
- Restarted from photo `x037` (`./re x037`), not run from scratch.
- Only change from the parent run: wind mass-loss rate enhanced by 2% when surface T < 1e4 K
  (`x_ctrl(3) = 0.02d0`, `x_ctrl(4) = 1d4` in `inlist_both`; applied in `VMS_framework` in `src/run_star_extras.f90`).
- First submitted 2026-11-06 as SLURM job 24559029.
