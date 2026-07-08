#!/usr/bin/env bash
# Phase 2 — Conserved loci extraction with Maffilter and phastCons
# Inputs:  $HAL or 01_cactus_alignment/${PRE}.hal and $CHR_GENOMES
# Outputs: 02_conserved_loci/conserved_loci/*.buff${BLOCK_LENGTH}.merge.fasta

source /opt/gabbi/utils/gabbi_functions.sh
source /opt/gabbi/utils/checkpoint.sh

mkdir -p 02_conserved_loci
cd 02_conserved_loci

if [[ -z "$HAL" ]]; then
    HAL="$GABBI_WORKDIR/01_cactus_alignment/${PRE}.hal"
fi

# List of genomes to base MAFs on (default all)
if [[ -z "$MAF_REFS" ]] ;then 
    halStats "$HAL"|awk -F"[ ,]" '$3==0{ print $1 }' > maf_references.txt
    MAF_REFS="$(realpath maf_references.txt)"
else
    cat "$MAF_REFS" > maf_references.txt
fi

# ---------------------------------------------------------------------------
# Step 2.1 — hal2maf conversion
# ---------------------------------------------------------------------------
if checkpoint_done "step2.1_hal2maf"; then
    echo "[GABBI] Skipping step 2.1 — checkpoint found."
else
    if checkpoint_fail_exists "step2.1_hal2maf"; then
        rm -r cactus_logs maf
    fi

    echo "[GABBI] Step 2.1: Converting HAL to MAF and filtering duplicates..."

    mkdir -p cactus_logs maf


    parallel -j "$PARALLEL_JOBS" \
        cactus-hal2maf \
            cactus_logs/js_hal2maf_{} \
            "$HAL" \
            maf/${PRE}.{}.maf.gz \
            --filterGapCausingDupes \
            --refGenome {} \
            --chunkSize 500000 \
            --logFile maf/${PRE}.{}.log \
            --batchParallelHal2maf "$CORES_PER_JOB" \
            --batchCores "$CORES_PER_JOB" \
            --batchLogsDir cactus_logs/js_hal2maf_{}/batch-logs-hal2maf \
        :::: "$MAF_REFS" \
        || checkpoint_fail "step2.1_hal2maf"

    # Filter paralogous loci
    for genome in $(cat "$MAF_REFS"); do
        zcat "maf/${PRE}.${genome}.maf.gz" | mafDuplicateFilter -km - > "maf/${PRE}.${genome}.single-copy.maf"
    done

    checkpoint_mark "step2.1_hal2maf"
fi

# ---------------------------------------------------------------------------
# Step 2.2 — Split MAF by chromosome
# ---------------------------------------------------------------------------
if checkpoint_done "step2.2_split_maf"; then
    echo "[GABBI] Skipping step 2.2 — checkpoint found."
else
    if checkpoint_fail_exists "step2.2_split_maf"; then
        rm -rf bed
    fi

    echo "[GABBI] Step 2.2: Splitting MAF files by chromosome..."

    mkdir -p bed

    for genome in $(cat "$MAF_REFS"); do

        halStats --bedSequences "$genome" "$HAL" > "bed/${PRE}.${genome}.bed"

        mkdir -p "split_maf/${genome}"
        python3 /opt/gabbi/scripts/extract_chromosomes_from_maf.py \
            "maf/${PRE}.${genome}.single-copy.maf" \
            "bed/${PRE}.${genome}.bed" \
            "split_maf/${genome}"

        parallel gzip {} ::: split_maf/${genome}/*.maf

    done || checkpoint_fail "step2.2_split_maf"

    checkpoint_mark "step2.2_split_maf"
fi

# ---------------------------------------------------------------------------
# Step 2.3 — Maffilter
# ---------------------------------------------------------------------------
if checkpoint_done "step2.3_maffilter"; then
    echo "[GABBI] Skipping step 2.3 — checkpoint found."
else
    if checkpoint_fail_exists "step2.3_maffilter"; then
        rm -r maffilter
    fi

    echo "[GABBI] Step 2.3: Filtering alignments with Maffilter, allowing blocks of alignments >=$BLOCK_LENGTH nt and >=$BLOCK_SIZE genomes..."

    param="/opt/gabbi/config/maffilter.optionfile"

    for genome in $(cat "$MAF_REFS"); do
        mkdir -p "maffilter/${genome}"
    done

    parallel --plus -j "$THREADS" \
        maffilter SP='$(basename {//})' DATA={/..} BLOCK_SIZE="$BLOCK_SIZE" BLOCK_LENGTH="$BLOCK_LENGTH" param="$param" \
        ::: $(find split_maf -type f -name "*.maf.gz") \
        || checkpoint_fail "step2.3_maffilter"
    
    verbose "$(ls maffilter/*/*.maf*|head)"
    # Remove empty maf files (only two header lines) to avoid phastcons errors
    for i in $(find maffilter -type f -name "*.maf.gz");do
        if [[ $(zcat $i|wc -l) == 2 ]] ;then
            verbose "Removing empty $i file"
            rm $i
        fi
    done || checkpoint_fail "step2.3_maffilter"
    
    checkpoint_mark "step2.3_maffilter"
fi

# ---------------------------------------------------------------------------
# Step 2.4 — PhastCons
# ---------------------------------------------------------------------------
if checkpoint_done "step2.4_phastcons"; then
    echo "[GABBI] Skipping step 2.4 — checkpoint found."
else
    if checkpoint_fail_exists "step2.4_phastcons"; then
        rm -r phastcons
        find maf -type f -name "*single-copy.maf.gz" 2>/dev/null | grep -q "." \
            && parallel -j "$THREADS" gunzip {} ::: $(find maf -type f -name "*single-copy.maf.gz")
    fi

    echo "[GABBI] Step 2.4: Computing conservation scores with PhastCons..."

    TREE=$(halStats --tree "$HAL")

    for genome in $(cat "$MAF_REFS"); do
        mkdir -p "phastcons/${genome}"
    done
        
    parallel -j "$PARALLEL_JOBS" phyloFit \
        --tree \""$TREE"\" \
        --subst-mod REV \
        --out-root "{}/{/}" \
        "maf/${PRE}.{/}.single-copy.maf" \
        ::: phastcons/* \
        || checkpoint_fail "step2.4_phastcons"

    parallel -j "$PARALLEL_JOBS" gzip {} \
        ::: $(find maf -type f -name "*.single-copy.maf")

    find maffilter/ -type f -name "*.maf.gz" 2>/dev/null | grep -q "." \
            && parallel -j "$THREADS" gunzip {} ::: $(find maffilter/ -type f -name "*.maf.gz" )

    # parallelize on parallel jobs instead of threads to avoid OOM kill on clusters
    parallel --plus -j "$PARALLEL_JOBS" '
        genome=$(basename {//})
	phastCons \
	    --msa-format MAF \
	    --most-conserved "phastcons/${genome}/{/.}.phastcons.bed" \
	    {} \
	    "phastcons/${genome}/${genome}.mod" > "phastcons/${genome}/{/.}.phastcons.wig" \
	' ::: $(find maffilter -type f -name "*.maf") \
	|| checkpoint_fail "step2.4_phastcons"

    parallel -j "$THREADS" gzip {} \
        ::: $(find maffilter -type f "*.maf")

    # Normalise BED output format for Phyluce
    for genome in $(cat "$MAF_REFS"); do
	cat "phastcons/${genome}/"*.phastcons.bed \
	    | sed -E -e "s/\.b[0-9]+//g" -e "s/\.l[0-9]+//g" \
	    | awk -F'\t' '{print $1"\t"$2"\t"$3}' \
	    > "phastcons/${genome}/${PRE}.${genome}.bed"
    done

    checkpoint_mark "step2.4_phastcons"
fi

# ---------------------------------------------------------------------------
# Step 2.5 — Extract conserved loci sequences
# ---------------------------------------------------------------------------
if checkpoint_done "step2.5_conserved_loci"; then
    echo "[GABBI] Skipping step 2.5 — checkpoint found."
else
    if checkpoint_fail_exists "step2.5_conserved_loci"; then
        rm -rf 2bit_genomes conserved_loci phyluce_logs
    fi
    echo "[GABBI] Step 2.5: Extracting conserved loci sequences..."

    mkdir -p 2bit_genomes conserved_loci phyluce_logs
    
    source /opt/miniconda3/etc/profile.d/conda.sh
    conda init --all
    conda activate phyluce

    for genome in $(cat "$MAF_REFS"); do
        
        if [[ -n "$CHR_GENOMES" ]]; then
            genome_fasta=$(find "$CHR_GENOMES/${genome}" \( -name "*.fasta" -o -name "*.fas" -o -name "*.fna" \) -type f)
        else
            echo "[GABBI] Getting $genome genome from $HAL..."
            hal2fasta "$HAL" "$genome" > "2bit_genomes/${genome}.fasta" \
            || checkpoint_fail "step2.5_conserved_loci"
            genome_fasta="2bit_genomes/${genome}.fasta"
        fi

        verbose "genome_fasta=$genome_fasta"

        faToTwoBit \
            "$genome_fasta" \
            "2bit_genomes/${genome}.2bit" \
            || checkpoint_fail "step2.5_conserved_loci"

        # Discard loci shorter than 20 bp
        awk '{if ($3-$2 >= 20) print $0}' \
            "phastcons/${genome}/${PRE}.${genome}.bed" \
            > "conserved_loci/${PRE}.${genome}.20.bed"

        phyluce_probe_get_genome_sequences_from_bed \
            --bed "conserved_loci/${PRE}.${genome}.20.bed" \
            --twobit "2bit_genomes/${genome}.2bit" \
            --buffer-to "$BLOCK_LENGTH" \
            --output "conserved_loci/${PRE}.${genome}.20.buff${BLOCK_LENGTH}.fasta" \
            --log-path phyluce_logs \
            || checkpoint_fail "step2.5_conserved_loci"

        verbose "Merging overlapping sequences from bed files..."
        python3 /opt/gabbi/scripts/merge_overlapping_seq_in_fasta.py \
            "conserved_loci/${PRE}.${genome}.20.buff${BLOCK_LENGTH}.fasta" \
            "conserved_loci/${PRE}.${genome}.20.buff${BLOCK_LENGTH}.merge.fasta" \
            || checkpoint_fail "step2.5_conserved_loci"
    done

    for fasta in $(find conserved_loci -type f -name "*buff${BLOCK_LENGTH}.merge.fasta"); do
        verbose "Converting $fasta to one-liner FASTA..."
        cat "$fasta" \
            | sed -e '/>/s/$/#/g' -e '/>/s/^/#/g' \
            | tr -d '\n' | tr '#' '\n' | sed '/^$/d' \
            > "${fasta/.fasta}.ok.fasta"
    done

    checkpoint_mark "step2.5_conserved_loci"

fi

echo "[GABBI] Removing temporary MAF files..."
find maf -type f -name "*.maf.gz" -print -delete 
find maffilter -type f -name "*.maf.gz" -print -delete
rm -r split_maf

echo "[GABBI] Removing temporary phastCons files..."
find phastcons -type f -name "*phastcons*" -print -delete

echo "[GABBI] ============================================================"
echo "[GABBI] PHASE 2: Completed at $(date '+%Y/%m/%d %H:%M:%S')"
echo "[GABBI] Minimum alignment block size:                   $BLOCK_SIZE out of ${N_NODES} including ancestral genomes"
echo "[GABBI] Minimum alignment block length:                 $BLOCK_LENGTH"
echo "[GABBI] Conserved loci found by phastcons:              ${OUT:-GABBI_out}/02_conserved_loci/conserved_loci/*.merge.fasta"
echo "[GABBI] Note: block size and length can be changed with --block-size and --block-length options."
echo "[GABBI] ============================================================"

cd "$GABBI_WORKDIR"
