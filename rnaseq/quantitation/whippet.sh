#!/bin/bash

set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 <input directory>"
fi

WHIPPET_DIR="$HOME/projects/Whippet.jl/bin"

TARGET_DIRECTORY=$1

regex="./(.+_DM_.+)_L003_R1_001.fastq.gz"
num_procs=0

for R1 in $TARGET_DIRECTORY/*; do
    if [[ $R1 =~ $regex ]]; then
        name=${BASH_REMATCH[1]}
        julia "$WHIPPET_DIR/whippet-quant.jl" \
              "${name}_L003_R1_001.fastq.gz" "${name}_L003_R2_001.fastq.gz" \
              -o $name &> "$name.log" &
        pids[$((num_procs++))]=$!
    fi
done

# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done

echo "All done."
