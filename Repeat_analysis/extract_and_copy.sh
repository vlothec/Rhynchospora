#!/bin/bash

# Ensure correct usage
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 SOURCE_DIRECTORY TARGET_DIRECTORY"
  exit 1
fi

# Assign arguments to variables
SOURCE_DIR=$1
TARGET_DIR=$2

# Check if source directory exists
if [ ! -d "$SOURCE_DIR" ]; then
  echo "Source directory does not exist: $SOURCE_DIR"
  exit 1
fi

# Create the target directory if it does not exist
mkdir -p "$TARGET_DIR"

# Find and copy files, preserving the directory structure
find "$SOURCE_DIR" -type f -name '*_helixer_filtered.csv' | while read -r file; do # *_repeats_filtered.csv *_edta_filtered.csv *_helixer_filtered.csv *_centromeric_arrays.csv cen_locations_plot_* (can be meta holo or other) *_genome_metadata.csv *satellite_te_transition_data_matrix.csv *_repeats_filtered_array_assigned.csv
  # Extract the relative path of the file with respect to the source directory
  rel_path=$(realpath --relative-to="$SOURCE_DIR" "$file")
  
  # Determine the directory structure to create in the target directory
  target_dir="$TARGET_DIR"
  
  # Create the necessary subdirectory in the target directory
  mkdir -p "$target_dir"
  
  # Copy the file to the target directory
  cp "$file" "$target_dir"
  
  echo "Copied: $file -> $target_dir"
done

echo "Files have been copied to $TARGET_DIR, preserving the directory structure."