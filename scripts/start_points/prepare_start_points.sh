#!/bin/bash

# A robust shell script to reformat a specific CSV file format.
# It selects specific columns, renames the header, and saves a new file.
#
# Usage:
#   ./reformat_csv.sh /path/to/your/inputfile.csv

# --- 1. Input Validation ---
# Check if a file path was provided as an argument.
if [ -z "$1" ]; then
    echo "Error: No input file provided."
    echo "Usage: $0 <path_to_input_csv>"
    exit 1
fi

INPUT_FILE="$1"

# Check if the provided file actually exists.
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File not found at '$INPUT_FILE'"
    exit 1
fi

# --- 2. Define Output File Name ---
# Create a new filename based on the input filename.
OUTPUT_FILE="${INPUT_FILE%.csv}_reformatted.csv"

# --- 3. Process the file with awk ---
# awk is a powerful tool for column-based data manipulation.
awk '
BEGIN {
    # Set the input Field Separator (FS) to a comma, surrounded by any
    # number of spaces. This correctly handles headers like ", col1 , col2".
    FS = " *, *"
    # Set the Output Field Separator (OFS) to a single comma.
    OFS = ","
}
# This block runs only for the first line (NR==1), which is the header.
NR == 1 {
    # Manually print the new, correct header.
    printf "time"
    for (i = 1; i <= 200; i++) {
        printf "%s", OFS "f" i
    }
    printf "\n"
    next # Tells awk to stop processing this line and move to the next one.
}
# This block runs for all other lines (the data rows).
{
    # Print the first field ($1), which corresponds to "Total Simulation Time".
    printf "%s", $1
    # Loop from the 5th field to the 204th, which correspond to the "freq at X" columns.
    for (j = 5; j <= 204; j++) {
        printf "%s", OFS $j
    }
    printf "\n"
}
' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Reformatting complete."
echo "Clean data saved to: $OUTPUT_FILE"


