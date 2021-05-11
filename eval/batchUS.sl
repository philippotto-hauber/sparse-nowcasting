#!/bin/bash
#SBATCH --job-name=evalUS
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=180gb
#SBATCH --time=02:00:00
#SBATCH --output=evalUS.out
#SBATCH --error=evalUS.err
#SBATCH --partition=cluster

module load matlab/2020b

for i in `seq 1 24`
do
	./batch_eval_hpc $i 1&
done
wait


