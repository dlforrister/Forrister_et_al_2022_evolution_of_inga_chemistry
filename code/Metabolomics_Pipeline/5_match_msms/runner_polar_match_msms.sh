#!/bin/bash

for i in {1..177}
do
    echo "job $i" #Here you put sbatch mysl.sl
    sbatch polar_match_msms.sl $i
done
echo "DONE"
