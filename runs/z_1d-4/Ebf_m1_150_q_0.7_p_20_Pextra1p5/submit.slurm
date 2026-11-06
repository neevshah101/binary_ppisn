#!/bin/bash
#SBATCH --job-name=z1d-4_Ebf_m1_150_q_0.7_p_20_Pextra1p5
#SBATCH --account=binboom
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --output=/groups/mrenzo/Neev/Projects/Pair-Instability/binary_ppisn/runs/z_1d-4/Ebf_m1_150_q_0.7_p_20_Pextra1p5/slurm-%j.out

echo "=== Job $SLURM_JOB_ID on $(hostname), start: $(date) ==="
export OMP_NUM_THREADS=16
cd "/groups/mrenzo/Neev/Projects/Pair-Instability/binary_ppisn/runs/z_1d-4/Ebf_m1_150_q_0.7_p_20_Pextra1p5" || exit 1
make clean && ./rn
echo "=== End: $(date) ==="
