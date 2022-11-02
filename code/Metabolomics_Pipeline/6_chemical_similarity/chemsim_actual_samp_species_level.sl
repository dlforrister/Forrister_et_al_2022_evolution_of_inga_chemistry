#!/bin/bash
#SBATCH --account=coley
#SBATCH --partition=lonepeak
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=abrsoule@gmail.com
#SBATCH --job-name=chemsim_nocutoff_species


export WORK_DIR=/uufs/chpc.utah.edu/common/home/inga-group1/Figures_Evol_Inga_Chemistry/Evolving_Samples/
export FILENAME=6_chemsim_allcomps_actual_sample_species_AS.R

ml RStudio/1.1.463

cd $WORK_DIR

Rscript $FILENAME > setup.$SLURM_JOBID.out

