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

# Remove duplicate assemblies sharing the same GenBank identifier
# When both GCA and GCF are present we keep the GCA record and drop the GCF
filtered_table=$(mktemp)
trap 'rm -f "$filtered_table"' EXIT

awk -v gb="$genbank_col" 'BEGIN{FS=OFS="\t"}
NR==1 { print; next }                       # keep header untouched
{
    core = $gb
    sub(/^GC[AF]_/, "", core)               # strip GCA_/GCF_ prefix
    sub(/\..*$/,    "", core)               # strip .version suffix
    lines[NR] = $0
    accs[NR]  = $gb
    cores[NR] = core
    if (substr($gb, 1, 3) == "GCA") has_gca[core] = 1
    order[++n] = NR
}
END {
    for (i = 1; i <= n; i++) {
        r = order[i]
        if (substr(accs[r], 1, 3) == "GCF" && has_gca[cores[r]]) {
            printf("WARNING: duplicate identifier %s -> dropping RefSeq record %s (keeping the GCA record)\n", \
                   cores[r], accs[r]) > "/dev/stderr"
            continue
        }
        print lines[r]
    }
}' "$table" > "$filtered_table"

# Output subdirectories
CHR_DIR="chr_level_genomes"
ADD_DIR="additional_genomes"
mkdir -p "$CHR_DIR" "$ADD_DIR"

WORKDIR=$(pwd)

# Process data lines (skip header)
awk 'BEGIN{FS=OFS="\t"} {for(i=1;i<=NF;i++) if($i=="") $i="NA"; print}' "$filtered_table" | tail -n +2 | while IFS=$'\t' read -ra fields; do

    genbank="${fields[$((genbank_col - 1))]}"
    species_name="${fields[$((species_col - 1))]}"
    level="${fields[$((level_col - 1))]}"

    # Skip empty lines if any
    [ -z "$genbank" ] && continue

    # Reformat genome name (replace all spaces with underscores)
    sp_GCA="${species_name// /_}_${genbank//./_}"

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
        unzip -q ncbi_dataset.zip
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
