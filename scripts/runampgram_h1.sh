#!/bin/bash
#PBS -j oe
#PBS -N runampgram_homo1
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=8:mem=80gb
#PBS -m ae
#PBS -J 1-100:1

cd $PBS_O_WORKDIR
shopt -s expand_aliases
source /etc/profile.d/modules.sh
echo "Job identifier is $PBS_JOBID"
echo "Working directory is $PBS_O_WORKDIR"

module load anaconda3
source $CONDA_PROF/conda.sh
conda activate R-403

echo $PBS_ARRAY_INDEX
echo $PBS_O_WORKDIR

cd /home/user/ampgram_homo

Rscript runampgram_h1.R ampin${PBS_ARRAY_INDEX}.fasta