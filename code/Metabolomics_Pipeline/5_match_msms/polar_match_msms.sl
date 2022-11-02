#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --partition=notchpeak
#SBATCH --account=coley
#SBATCH --time=2:00:00

echo "Start Date:`date`"
echo "Program:$0"
echo "#Arguments:$#"
export SPECIES_ID=$1
echo " My species id is: $SPECIES_ID "
cd K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/
ml RStudio/1.1.463
K:/Lab_Map/git_repositories/Evolution_Of_Inga_Chemistry-master/code/5_match_compounds_to_msms_spectra_polar.R $SPECIES_ID
echo " Return statement:$?"
echo "End Date:`date`"
