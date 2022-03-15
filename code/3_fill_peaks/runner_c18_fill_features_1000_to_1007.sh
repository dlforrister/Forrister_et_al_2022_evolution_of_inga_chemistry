#!/bin/bash

for i in {1000..1007}

do
    echo "job $i" #Here you put sbatch mysl.sl
    sbatch ../c18_fill_features.sl $i
done
echo "DONE"
