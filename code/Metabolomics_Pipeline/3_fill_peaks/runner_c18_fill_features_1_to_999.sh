#!/bin/bash

for i in {1..999}

do
    echo "job $i" #Here you put sbatch mysl.sl
    sbatch ../C18_fill_featuress.sl $i
done
echo "DONE"
