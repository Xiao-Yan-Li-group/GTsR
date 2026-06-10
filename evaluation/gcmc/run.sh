#!/bin/bash

BASE_DIR="$1"
NP="${2:-1}"

echo "BASE_DIR = $BASE_DIR"
echo "NP = $NP"

running=0

for job_dir in "$BASE_DIR"/*/; do
    if [[ -f "$job_dir/run.py" ]]; then
        (
            cd "$job_dir" || exit
            echo "Running: $job_dir"
            python run.py
        ) &

        ((running++))

        if (( running >= NP )); then
            wait -n
            ((running--))
        fi
    fi
done

wait
