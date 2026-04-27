#!/usr/bin/env bash
# Phase 3 — Cross-BLASTn and provisional SHR loci identification
# Inputs:  02_conserved_loci/conserved_loci/*.buff${BLOCK_LENGTH}.merge.ok.fasta
# Outputs: 03_cross_blast/cactus_chr_level.raw_shr_from_blastn.filtered.list

source /opt/gabbi/utils/gabbi_functions.sh
source /opt/gabbi/utils/checkpoint.sh

mkdir -p 03_cross_blast
cd 03_cross_blast

MAF_REFS="$(realpath "$GABBI_WORKDIR"/02_conserved_loci/maf_references.txt)"

# ---------------------------------------------------------------------------
# Step 3.1 — Build BLAST databases and run all-vs-all cross-BLASTn
# ---------------------------------------------------------------------------
if checkpoint_done "step3.1_cross_blast"; then
    echo "[GABBI] Skipping step 3.1 — checkpoint found."
else
    if checkpoint_fail_exists "step3.1_cross_blast"; then
        rm -rf blast_db cross_blast
    fi

    echo "[GABBI] Step 3.1: Running cross-BLASTn between phastcons results..."

    mkdir -p blast_db
    find "$GABBI_WORKDIR"/02_conserved_loci/conserved_loci/ -type f -name "*.merge.ok.fasta" -exec ln -s {} blast_db \;

    parallel --plus -j $PARALLEL_JOBS \
        makeblastdb -dbtype nucl -in {} \
        ::: $(find -L blast_db/ -type f -name "*.merge.ok.fasta") \
        || checkpoint_fail "step3.1_cross_blast"

    mkdir -p cross_blast

    # To catch blast error
    set +e

    for ref in $(cat "$MAF_REFS"); do
        for query in $(cat "$MAF_REFS"); do
            [[ "$query" == "$ref" ]] && continue
            debug "BLASTn: $query on $ref"
            {
            blastn \
                -db "blast_db/${PRE}.${ref}.20.buff${BLOCK_LENGTH}.merge.ok.fasta" \
                -query "blast_db/${PRE}.${query}.20.buff${BLOCK_LENGTH}.merge.ok.fasta" \
                -word_size "$CROSS_BLAST_WS" \
                -qcov_hsp_perc "$CROSS_BLAST_QC" \
                -outfmt 6 \
                > "cross_blast/${query}_on_${ref}.w${CROSS_BLAST_WS}.blastn" 
            } || checkpoint_fail "step3.1_cross_blast" 
            # Filter hits by e-value
            awk -v ev="$CROSS_BLAST_EV" \
                '{if (($11 + 0) < ev) print $0}' \
                "cross_blast/${query}_on_${ref}.w${CROSS_BLAST_WS}.blastn" \
                > "cross_blast/${query}_on_${ref}.w${CROSS_BLAST_WS}.ev${CROSS_BLAST_EV}.blastn"
            # Merge all filtered hits into a single table
            awk -v q="$query" -v r="$ref" \
                '{print $1 "\t" q "\t" $2 "\t" r}' \
                "cross_blast/${query}_on_${ref}.w${CROSS_BLAST_WS}.ev${CROSS_BLAST_EV}.blastn" \
                >> cross_blast/cross_blast.w${CROSS_BLAST_WS}.ev${CROSS_BLAST_EV}.table
        done
    done

    set -e
    
    checkpoint_mark "step3.1_cross_blast"
fi

# ---------------------------------------------------------------------------
# Step 3.2 — Cluster BLAST hits into provisional SHR loci
# ---------------------------------------------------------------------------
if checkpoint_done "step3.2_shr_clustering"; then
    echo "[GABBI] Skipping step 3.2 — checkpoint found."
else
    if checkpoint_fail_exists "step3.2_shr_clustering"; then
        rm -rf shr_clustering
    fi
    
    echo "[GABBI] Step 3.2: Clustering cross-BLASTn hits into provisional SHR loci..."

    mkdir -p shr_clustering

    python3 /opt/gabbi/scripts/regroup_matches_from_blastn.py \
        cross_blast/cross_blast.w${CROSS_BLAST_WS}.ev${CROSS_BLAST_EV}.table \
        shr_clustering/${PRE}.raw_shr_from_blastn.list \
        || checkpoint_fail "step3.2_shr_clustering"

    # Annotate loci with per-locus taxon count and duplication count
    awk '
    {
        shr = $3; species = $2
        lines[NR] = $0
        shr_species[shr][species]++
        shr_lines[NR] = shr
    }
    END {
        for (i = 1; i <= NR; i++) {
            shr = shr_lines[i]
            distinct = 0; duplicates = 0
            for (s in shr_species[shr]) {
                distinct++
                if (shr_species[shr][s] > 1) duplicates++
            }
            print lines[i] "\t" distinct "\t" duplicates
        }
    }
    ' shr_clustering/${PRE}.raw_shr_from_blastn.list \
        > shr_clustering/${PRE}.raw_shr_from_blastn.stat.list

    debug "Total raw SHR loci: $(wc -l < shr_clustering/${PRE}.raw_shr_from_blastn.stat.list)"

    # Apply taxon count and duplication filters
    # Column 4: distinct taxa count; Column 5: duplicated taxa count
    awk -v mt="$MIN_TAXA_ABS" -v md="$TEMP_ALLOW_DUPES" \
        '{if (($4 >= mt) && ($5 <= md)) print $0}' \
        shr_clustering/${PRE}.raw_shr_from_blastn.stat.list \
        > shr_clustering/${PRE}.shr_from_blastn.mintax${MIN_TAXA_ABS}.dupes${TEMP_ALLOW_DUPES}.list \
        || checkpoint_fail "step3.2_shr_clustering"

    echo "[GABBI] SHR loci after filtering: $(wc -l < shr_clustering/${PRE}.shr_from_blastn.mintax${MIN_TAXA_ABS}.dupes${TEMP_ALLOW_DUPES}.list)"

    checkpoint_mark "step3.2_shr_clustering"
fi

echo "[GABBI] ============================================================"
echo "[GABBI] PHASE 3: Completed at $(date '+%Y/%m/%d %H:%M:%S')"
echo "[GABBI] Minimum number of taxa to keep a SHR:        $MIN_TAXA_ABS (${TEMP_TAX_THRESHOLD}% of ${N_CHR_TAXA})"
echo "[GABBI] Maximum number of duplications allowed:      ${TEMP_ALLOW_DUPES}"
echo "[GABBI] Number of temporary SHR before filtering:    $(egrep -o "uce-[0-9]+" shr_clustering/${PRE}.raw_shr_from_blastn.list|sort -u|wc -l)"
echo "[GABBI] Number of temporary SHR after filtering:     $(wc -l < shr_clustering/${PRE}.shr_from_blastn.mintax${MIN_TAXA_ABS}.dupes${TEMP_ALLOW_DUPES}.list)"
echo "[GABBI] Note: these thresholds can be changed with --temp-tax-threshold and --temp-allow-dupes options."
echo "[GABBI] ============================================================"

cd "$GABBI_WORKDIR"
