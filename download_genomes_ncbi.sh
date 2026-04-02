#!/bin/bash

# install ncbi_datasets: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/

# Download NCBI table of genomes you want from https://www.ncbi.nlm.nih.gov/datasets/genome/
# and give it as the first argument

table=$1

# Read the table
cat $table| while IFS=$'\t' read -r genbank assembly species_name _ _ _ _ _ _ _ _ _; do

    if [ "$genbank" = "Assembly Accession" ]; then
        continue
    fi

    # Reformat genome name
    sp_GCA="${species_name/ /_}_$genbank"

    echo "Downloading ${sp_GCA}"
    if [ ! -d "${sp_GCA}" ];then
        mkdir -p "$sp_GCA"
        cd "$sp_GCA"

        # Download and extract .fna
        datasets download genome accession "$genbank" --include genome
        sleep 5
        unzip ncbi_dataset.zip
        sleep 5
        mv ncbi_dataset/*/*/*.fna .
        rm -rf ncbi* md5sum.txt README.md
        cd ..
    else 
        echo "${sp_GCA} already exists"
    fi
done

