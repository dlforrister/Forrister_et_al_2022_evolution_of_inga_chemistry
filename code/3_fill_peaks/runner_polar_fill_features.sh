#!/bin/bash

for i in {1..812}

do
    echo "job $i" #Here you put sbatch mysl.sl
    sbatch ../polar_fill_features.sl $i
done
echo "DONE"
