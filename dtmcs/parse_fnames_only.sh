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
        set -- $line
	dir_name=$(echo $1 | sed 's/[0-9].*//')
        drn_fname="${1//pm/drn}"
 	drdd_fname="${1//pm/drdd}"       
        storm --prism "$dir_name"/$line --buildfull --prismcompat --engine sparse --exportbuild $dir_name/$drn_fname
        storm --prism "$dir_name"/$line --buildfull --prismcompat --engine dd --exportbuild $dir_name/$drdd_fname
    fi
done < "$input_file"

