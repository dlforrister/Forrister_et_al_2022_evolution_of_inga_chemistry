#!/bin/bash

for i in {1..98}
do
    echo "job $i" #Here you put sbatch mysl.sl
    sbatch ../6_chemsim_allcomps_null_model.sl $i
done
echo "DONE"
