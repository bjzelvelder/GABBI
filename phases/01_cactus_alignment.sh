#!/usr/bin/env bash
# Phase 1 — Whole-genome alignment with Cactus
# Inputs:  $CHR_GENOMES, $GUIDE_TREE
# Outputs: cactus_chr_level.hal

source /opt/gabbi/utils/gabbi_functions.sh
source /opt/gabbi/utils/checkpoint.sh

if checkpoint_done "step1.1_cactus_alignment"; then
    echo "[GABBI] Skipping step 1 — checkpoint found."
else
    if checkpoint_fail_exists "step1.1_cactus_alignment"; then
        rm -r 01_cactus_alignment
    fi  

    echo "[GABBI] Preparing Cactus input files..."

    mkdir -p 01_cactus_alignment
    cd 01_cactus_alignment

    # Collect absolute genome paths and species names
    find "$CHR_GENOMES" -type f \( -name "*.fasta" -o -name "*.fna" -o -name "*.fas" \) -exec realpath {} \; > cactus_genomes.paths
    
    echo > cactus_genomes.spp
    echo > cactus_genomes.paths
    for i in "$CHR_GENOMES"/*/; do
        basename "$i" >> cactus_genomes.spp
        find "$i" -type f \( -name "*.fasta" -o -name "*.fna" -o -name "*.fas" \) -exec realpath {} \; >> cactus_genomes.paths
    done

    verbose "Genome paths:"
    verbose "$(cat cactus_genomes.paths)"

    cat "$GUIDE_TREE" \
        <(echo) \
        <(paste cactus_genomes.spp cactus_genomes.paths) \
        > cactus_input.txt

    if [[ $(grep -o -f cactus_genomes.spp "$GUIDE_TREE"|wc -l) -ne $N_CHR_TAXA ]]; then
        echo "[GABBI] ERROR: Unable to find input genomes in guide tree. See ${OUT:-GABBI_out}/01_cactus_alignment/cactus_input.txt"
        checkpoint_fail "step1.1_cactus_alignment"
    fi

    echo "[GABBI] Step 1.1: Running Cactus whole-genome alignment..."

    mkdir -p cactus_tmp

    cactus ./js \
	cactus_input.txt \
	${PRE}.hal \
	--workDir cactus_tmp \
	--maxCores "$CACTUS_MAXCORES" \
	--maxDisk "$CACTUS_MAXDISK" \
	--maxMemory "$CACTUS_MAXMEM" \
	--batchLogsDir cactus_logs

#    echo "[GABBI] Cactus HAL tree:"
#    halStats --tree cactus_chr_level.hal
    if [ ! -f ${PRE}.hal ]; then
        checkpoint_fail "step1.1_cactus_alignment"
    fi
    
    verbose "Genome names = $(halStats --genomes ${PRE}.hal)"

    checkpoint_mark "step1.1_cactus_alignment"
fi

echo "[GABBI] ============================================================"
echo "[GABBI] PHASE 1: Completed at $(date '+%Y/%m/%d %H:%M:%S')"
echo "[GABBI] Cactus input:             ${OUT:-GABBI_out}/01_cactus_alignment/cactus_input.txt"
echo "[GABBI] Genome alignment:         ${OUT:-GABBI_out}/01_cactus_alignment/${PRE}.hal"
echo "[GABBI] ============================================================"

cd "$GABBI_WORKDIR"
