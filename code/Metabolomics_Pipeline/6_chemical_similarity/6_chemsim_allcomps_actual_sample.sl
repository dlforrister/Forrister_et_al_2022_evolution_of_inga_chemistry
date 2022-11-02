#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=notchpeak
#SBATCH --account=coley
#SBATCH --time=6:00:00

echo "Start Date:`date`"
echo "Program:$0"
echo "#Arguments:$#"
export SPECIES_ID=$1
echo " My species id is: $SPECIES_ID "

export WORK_DIR=/uufs/chpc.utah.edu/common/home/inga-group1/4_directories_chem_evolution_2019_03_26/
export FILENAME=/uufs/chpc.utah.edu/common/home/inga-group1/Figures_Evol_Inga_Chemistry/Evolving_Samples/6_chemsim_allcomps_actual_sample.R

cd $WORK_DIR
ml RStudio/1.1.463

Rscript $FILENAME $SPECIES_ID > setup.$SPECIES_ID.$SLURM_JOBID.out

echo " Return statement:$?"
echo "End Date:`date`"
