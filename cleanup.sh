#!/bin/bash

# Define list of SRR IDs
srr_ids=(
    # SRR31810741
    # SRR31810742
    # SRR31814601
    # SRR31814602
    # SRR31814603
    # SRR31814604
    # SRR31814607
    # SRR31814608
    # SRR31815266
    # SRR31815269
    # SRR31885140
    # SRR31885141
    # SRR31885142
    # SRR31885143
    # SRR31885144
)

# Base path
base_path="/storage2/egusmao/projects/Bloom/data/raw/GSE285812"

# Iterate over the SRR IDs
for SRRID in "${srr_ids[@]}"; do
    echo "Processing $SRRID"

    # Define file paths
    script_file="$base_path/prefetch_${SRRID}.sh"
    log_file="$base_path/prefetch_${SRRID}.log"
    # data_file_pattern="$base_path/${SRRID}_index.html?view=run_browser&acc=${SRRID}&display=data-access*"

    # Try to delete each file, if it exists
    for file in "$script_file" "$log_file"; do
        if [ -e "$file" ]; then
            echo "Found: $file"
            rm "$file"
            if [ ! -e "$file" ]; then
                echo "Deleted: $file"
            else
                echo "Failed to delete: $file"
            fi
        else
            echo "Not found: $file"
        fi
    done

    # Delete matching files with wildcard (like the malformed output file)
    for file in $data_file_pattern; do
        if [ -e "$file" ]; then
            echo "Found: $file"
            rm "$file"
            if [ ! -e "$file" ]; then
                echo "Deleted: $file"
            else
                echo "Failed to delete: $file"
            fi
        else
            echo "Not found: $file"
        fi
    done

    echo ""
done

