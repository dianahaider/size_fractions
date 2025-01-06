#!/bin/bash

# Define input and output directories
INPUT_DIR="bib_refs"
OUTPUT_DIR="bib_refs"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through all .ris files in the input directory
for ris_file in "$INPUT_DIR"/*.ris; do
    # Extract filename without extension
    filename=$(basename "$ris_file" .ris)
    
    # Define output .bib file path
    bib_file="$OUTPUT_DIR/$filename.bib"
    
    # Convert .ris to .bib using pandoc
    pandoc "$ris_file" -f ris -t bibtex -o "$bib_file"
    
    # Check if the conversion was successful
    if [[ $? -eq 0 ]]; then
        # Replace "and" with ", " in the author field
        sed -i.bak -E 's/author\s*=\s*\{([^}]*)\}/author = {\1}/; s/ and /, /g' "$bib_file"
        rm "$bib_file.bak"  # Remove backup created by sed
        echo "Converted and reformatted: $ris_file -> $bib_file"
    else
        echo "Failed to convert: $ris_file"
    fi
done

echo "Conversion completed."
