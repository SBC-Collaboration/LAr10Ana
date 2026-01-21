#!/bin/bash

# Quit if error happens, but don't kill the entire terminal in jupyter.
(
set -e
# Directory where run data are copied to for the grid job
TEMP_DIR="/pnfs/coupp/scratch/users/${USER}/temp_data"
# Directory where job output are saved to
OUT_DIR="/pnfs/coupp/scratch/users/${USER}/grid_output"
# Directory where finished job outputs are copied to
RECON_DIR="/exp/e961/data/SBC-25-recon/dev-output"
LOG_DIR="/exp/e961/data/SBC-25-recon/dev-logs"
# File containing all submitted jobs
JOBS_LIST="${HOME}/.cache/sbc_job_list.csv"

# Process each run
for folder in "$OUT_DIR"/20*_*-*_*; do
    [ -d "$folder" ] || continue
    
    folder_name=$(basename "$folder")
    # Extract run number and job ID (format: 20251120_0-38273648.0)
    run_num=$(echo "$folder_name" | cut -d'-' -f1)
    job_id=$(echo "$folder_name" | cut -d'-' -f2)
    
    log_file=$(find "$folder" -name "*.log" -type f | head -n1)
    echo "Processing output for run ${run_num}, job ${job_id}..."

    # Check exit code in log file
    if [ -n "$log_file" ] && [ -f "$log_file" ]; then
        if grep -q "EXIT CODE 0" "$log_file"; then
            # Success: check version before moving
            dest_folder="$RECON_DIR/$run_num"
            should_move=true
            
            if [ -d "$dest_folder" ]; then
                # Compare versions
                old_version=$(cat "$dest_folder/version.txt" 2>/dev/null || echo "")
                new_version=$(cat "$folder/version.txt" 2>/dev/null || echo "")
                
                if [ -z "$new_version" ]; then
                    should_move=false
                    echo "Skipped $folder_name (no version info)"
                elif [ "$new_version" = "$old_version" ]; then
                    should_move=true
                    echo "Moving $folder_name (same version)"
                else
                    # Compare versions using sort -V
                    latest=$(printf "%s\n%s" "$old_version" "$new_version" | sort -V | tail -n1)
                    if [ "$latest" = "$new_version" ]; then
                        should_move=true
                        echo "Upgrading from $old_version to $new_version"
                    else
                        should_move=false
                        echo "Skipping $folder_name (existing: $old_version, incoming: $new_version)"
                    fi
                fi
            fi
            
            if [ "$should_move" = true ]; then
                rm -rf "$dest_folder"
                mv "$folder" "$dest_folder"
                echo "Moved $folder_name to $dest_folder"
            else
                rm -rf "$folder"
            fi

            # Delete corresponding tar file
            tar_file="$TEMP_DIR/${run_num}.tar"
            if [ -f "$tar_file" ]; then
                rm "$tar_file"
                echo "Deleted $tar_file"
            fi

        else
            # Failed: delete the folder
            rm -rf "$folder"
            echo "Deleted $folder_name (non-zero exit code)"
        fi
    else
        echo "Warning: Log file not found for $folder_name"
    fi

    # Remove job from CSV list
    if [ -f "$JOBS_LIST" ]; then
        sed -i "/, ${run_num}, /d" "$JOBS_LIST"
        echo "Removed ${run_num} from job list"
    fi
done

)