#!/bin/bash

folder="/mnt/d/Project/GTSR/benchmark/stability/mlp/pacman_data/"

for mof in "$folder"*.cif; do
    echo "----------------------------------------"
    echo "run: $mof"

    network -ha -volpo 0 0 5000 "$mof"
    network -ha -res "$mof"
done