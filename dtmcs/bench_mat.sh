#!/bin/bash

# Usage: ./bench_mat.sh <input_file>
# Example: ./bench_mat.sh brp/models

if [ $# -lt 1 ]; then
    echo "Usage: $0 <input_file>"
    echo "Example: $0 brp/models"
    exit 1
fi

input_file="$1"

if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found"
    exit 1
fi

while IFS="" read -r line || [ -n "$line" ]
do
    # Skip empty lines and comment lines (starting with #)
    if [ -n "$line" ] && [ "${line#\#}" = "$line" ]; then
        # Extract first and third words from the line
        set -- $line
        model_name="${1//.pm/_}"
        const_str="${3//=/_}"
        const_str="${const_str//,/_}"
        drn_fname="$model_name$const_str".drn
        
        echo storm --prism "$line" --buildfull --prismcompat --engine sparse --exportbuild $drn_fname
        
    fi
done < "$input_file"

