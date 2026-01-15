#!/bin/bash

# Directory containing subdirectories with analysis data
data_dir="/home/pwlodzimierz/Rhynchospora/TRASH_2025"

# Initialize counters for each step
step0_count=0
step1_count=0
step2_count=0
step3_count=0
step4_count=0
step5_count=0

# Print table header
printf "%-70s %-8s %-8s %-8s %-8s %-8s %-8s\n" "Subdirectory" "TRASH" "Step 1" "Step 2" "Step 3" "Step 5" "Step 4"

# Loop through subdirectories
for subdir in "$data_dir"/*; do
    # Extract subdirectory name
    subdirname=$(basename "$subdir")

    # Check if step 0 file exists
    if [[ -f "$subdir"/"$subdirname"_repeats_with_seq.csv ]]; then
        step0=true
        ((step0_count++))
    else
        step0=false
    fi

    # Check if step 1 file exists
    if [[ -f "$subdir"/"${subdirname%.fa*}"_repeats_filtered.csv ]]; then
        step1=true
        ((step1_count++))
    else
        step1=false
    fi

    # Check if step 2 file exists
    if [[ -f "$subdir"/"$subdirname"_edta_modified.csv ]]; then
        step2=true
        ((step2_count++))
    else
        step2=false
    fi
	step3=false
    # Check if step 3 file exists
    if [[ -f "$subdir"/"${subdirname%.fasta}"_edta_filtered.csv ]]; then
        step3=true
        ((step3_count++))
    fi
	
    if [[ -f "$subdir"/"${subdirname%.fa}"_edta_filtered.csv ]]; then
        step3=true
        ((step3_count++))
    fi

	step5=false
    # Check if step 5 file exists
    if [[ -f "$subdir"/"${subdirname%.fasta}"_helixer_filtered.csv ]]; then
        step5=true
        ((step5_count++))
    fi
	
    if [[ -f "$subdir"/"${subdirname%.fa}"_helixer_filtered.csv ]]; then
        step5=true
        ((step5_count++))
    fi
    


    # Check if step 4 file exists
    if [[ -f "$subdir"/genome_summary_"$subdirname".csv ]]; then
        step4=true
        ((step4_count++))
    else
        step4=false
    fi

    # Print row for each subdirectory
    printf "%-70s %-8s %-8s %-8s %-8s %-8s %-8s\n" "$subdirname" "$step0"  "$step1" "$step2" "$step3" "$step5" "$step4"
done

# Print summary row
printf "%-70s %-8s %-8s %-8s %-8s %-8s %-8s\n" "==================" "====="  "=====" "=====" "=====" "=====" "====="
printf "%-70s %-8s %-8s %-8s %-8s %-8s %-8s\n" "" "TRASH"  "1.Reps" "2.EDTA" "3.EDTA" "5.helixer" "4.nothing"
printf "%-70s %-8s %-8s %-8s %-8s %-8s %-8s\n" "" "done"  "filtered" "modified" "filtered" "filtered" "here"

printf "%-70s %-8s %-8s %-8s %-8s %-8s %-8s\n" "Total:" "$step0_count"  "$step1_count" "$step2_count" "$step3_count" "$step5_count" "$step4_count"
