#!/bin/bash
#SBATCH -N 1    # nodes
#SBATCH -n 1    # cores
#SBATCH --mem=16G    # memory per node
##SBATCH --mem-per-cpu=6G    # memory per cpu
#SBATCH -p shared    # partition
#SBATCH -t 0-01:00    # runtime d-hh:mm
#SBATCH -o slurm/snakemake.out
#SBATCH -e slurm/snakemake.err
##SBATCH --mail-user=xxxxxxxx@xxxxx.com
##SBATCH --mail-type=END
##SBATCH --dependency=after:jobid[:jobid...] # job can begin after specified jobs have started
##SBATCH --dependency=afterany:jobid[:jobid...] # job can begin after specified jobs have terminated
##SBATCH --dependency=afterok:jobid[:jobid...] # job can begin after specified jobs have completed with exit code zero
##SBATCH --dependency=afternotok:jobid[:jobid...] # job can begin after specified jobs have failed

# Use snakemake environment with snakemake 7.18.2 installed (cannot see any env modules on cluster with snakemake)

# No need to load conda, as the Miniconda in my ~ is available from within a batch job.
#   module load Anaconda3/2022.05

# Note use of source activate (conda activate became preferred from conda v4.4)
    source activate snakemake

snakemake --use-conda -c 1 demultiplex_16S

# snakemake -c 8 --use-conda some_rule --forcerun upstream_rule

# snakemake -c 1 --use-conda -rerun-incomplete --unlock some_rule
# snakemake -c 8 --use-conda -rerun-incomplete some_rule
