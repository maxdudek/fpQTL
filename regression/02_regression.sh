#!/bin/bash
#SBATCH --job-name="regression_32"
#SBATCH -c 32
#SBATCH --mem=800G
#SBATCH -t 240:00:00
#SBATCH -o job_out/slurm.%x.%j.out
#SBATCH -e job_out/slurm.%x.%j.out

module load R/4.4.0

# for FP_METHOD in PRINT_no_gaussian PRINT; do
for FP_METHOD in PRINT_no_gaussian; do
    echo "Running regression for $FP_METHOD..."
    echo ""
    Rscript regression.R $FP_METHOD 32
    echo ""
done
