#!/bin/bash

# install ncbi_datasets: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/

# Download NCBI table of genomes you want from https://www.ncbi.nlm.nih.gov/datasets/genome/
# and give it as the first argument

table=$1

if [ -z "$table" ] || [ ! -f "$table" ]; then
    echo "Usage: $0 <ncbi_genome_table.tsv>"
    exit 1
fi

# Detect column indices from header
header=$(head -n 1 "$table")

get_col_index() {
    # Returns 1-based column index matching the given string
    echo "$header" | tr '\t' '\n' | grep -in "$1" | head -1 | cut -d: -f1
}

genbank_col=$(get_col_index "Assembly Accession")
species_col=$(get_col_index "Organism Name")
level_col=$(get_col_index "Assembly Level")

# Validate that required columns were found
if [ -z "$genbank_col" ] || [ -z "$species_col" ] || [ -z "$level_col" ]; then
    echo "Error: could not find required columns in header."
    echo "  Assembly Accession column : ${genbank_col:-NOT FOUND}"
    echo "  Organism Name column      : ${species_col:-NOT FOUND}"
    echo "  Assembly Level column     : ${level_col:-NOT FOUND}"
    echo ""
    echo "Header detected:"
    echo "$header" | tr '\t' '\n' | nl
    exit 1
fi

# Output subdirectories
CHR_DIR="chr_level_genomes"
ADD_DIR="additional_genomes"
mkdir -p "$CHR_DIR" "$ADD_DIR"

WORKDIR=$(pwd)

# Process data lines (skip header)
tail -n +2 "$table"|sed -E "s/\t\t/\tNA\t/g" | while IFS=$'\t' read -ra fields; do

    genbank="${fields[$((genbank_col - 1))]}"
    species_name="${fields[$((species_col - 1))]}"
    level="${fields[$((level_col - 1))]}"

    # Skip empty lines if any
    [ -z "$genbank" ] && continue

    # Reformat genome name (replace all spaces with underscores)
    sp_GCA="${species_name// /_}_${genbank}"

    if [ "$level" = "Chromosome" ]; then
        subdir="$CHR_DIR"
    else
        subdir="$ADD_DIR"
    fi

    dest="${WORKDIR}/${subdir}/${sp_GCA}"

    echo "Downloading ${sp_GCA} in ${subdir}/"

    if [ ! -d "$dest" ]; then
        mkdir -p "$dest"
        cd "$dest"

        datasets download genome accession "$genbank" --include genome
        sleep 5
        unzip ncbi_dataset.zip
        sleep 5
        mv ncbi_dataset/*/*/*.fna .
        rm -rf ncbi* md5sum.txt README.md

        cd "$WORKDIR"
    else
        echo "${subdir}/${sp_GCA} already exists, skipping."
    fi

done

echo ""
echo "Done."
echo "  Chromosome-level genomes : ${CHR_DIR}/"
echo "  Other genomes            : ${ADD_DIR}/"
