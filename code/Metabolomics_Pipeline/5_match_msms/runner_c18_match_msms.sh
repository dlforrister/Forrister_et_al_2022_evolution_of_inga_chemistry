#!/bin/bash

for i in {1..177}
do
    echo "job $i" #Here you put sbatch mysl.sl
    sbatch c18_match_msms.sl $i
done
echo "DONE"
