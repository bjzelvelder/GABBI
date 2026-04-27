#!/usr/bin/env bash
# Phase 4 — SHR sequence extraction and ancestral sequence reconstruction
# Inputs:  03_cross_blast/shr_clustering/${PRE}.shr_from_blastn.mintax${MIN_TAXA_ABS}.dupes${TEMP_ALLOW_DUPES}.list
# Outputs: 04_shr_extraction/${PRE}.anc.loci.fasta

source /opt/gabbi/utils/gabbi_functions.sh
source /opt/gabbi/utils/checkpoint.sh

mkdir -p 04_shr_extraction
cd 04_shr_extraction

# ---------------------------------------------------------------------------
# Step 4.1 — Extract per-locus FASTA sequences
# ---------------------------------------------------------------------------
if checkpoint_done "step4.1_shr_extraction"; then
    echo "[GABBI] Skipping step 4.1 — checkpoint found."
else
    if checkpoint_fail_exists "step4.1_shr_extraction"; then
        rm -rf shr
    fi

    echo "[GABBI] Step 4.1: Extracting SHR sequences..."

    mkdir -p shr

    cut -f1-3 "$GABBI_WORKDIR"/03_cross_blast/shr_clustering/${PRE}.shr_from_blastn.mintax${MIN_TAXA_ABS}.dupes${TEMP_ALLOW_DUPES}.list \
    | while read slice sp shr; do
        cat "$GABBI_WORKDIR"/03_cross_blast/blast_db/${PRE}.${sp}.20.buff${BLOCK_LENGTH}.merge.ok.fasta \
            | grep -A1 "${slice} " \
            | sed "s/>/>${sp}|/g" \
            | sed "/>/ s/ /|/g" \
            | sed "/^--$/d" \
            >> "shr/${shr}.temp.fasta"
    done \
    || checkpoint_fail "step4.1_shr_extraction"

    debug "SHR loci extracted: $(find shr/ -name "*.fasta" -type f | wc -l)"

    # Merge all per-locus FASTA into a single file
    parallel --plus -j "$THREADS" '\
        sed -E "s/>/>{..}|/g" {} \
        ' ::: $(find shr -type f) > ${PRE}.temp.loci.fasta \
    || checkpoint_fail "step4.1_shr_extraction"

    checkpoint_mark "step4.1_shr_extraction"
fi


# ---------------------------------------------------------------------------
# Step 4.2 — Compute gene trees
# ---------------------------------------------------------------------------
if checkpoint_done "step4.2_alignments"; then
    echo "[GABBI] Skipping step 4.2 — checkpoint found."

else
    if checkpoint_fail_exists "step4.2_alignments"; then
        rm -rf alignments
    fi

    echo "[GABBI] Step 4.2: Aligning temporary loci and building gene trees..."

    # Remove variation held by only one taxon (fix the AMAS trim threshold to 2 taxa out of total taxa)
    trim=$(awk -v n="$N_CHR_TAXA" 's=200/n { print int(s) }' <(echo ))
    debug "trim=$trim"

    phylomera \
        --input shr \
        --output alignments \
        --prefix ${PRE} \
        --trim $trim \
        --drop 0 \
        --genetrees MFP \
        --nogenepart \
        --perc 0 \
        --threads "$THREADS" \
        --debug \
        --continue \
        --config /opt/gabbi/config/phylomera.conf \
    || checkpoint_fail "step4.2_alignments"
    # debug option prevents TMP folder to be removed

    checkpoint_mark "step4.2_alignments"
fi

silence=false
# ---------------------------------------------------------------------------
# Step 4.3 — Reconstruct ancestral sequences
# ---------------------------------------------------------------------------
if checkpoint_done "step4.3_ancestral_seqs"; then
    echo "[GABBI] Skipping step 4.3 — checkpoint found."

else
    # Check the number of taxa before making ancestral seqs
    if [[ "$N_CHR_TAXA" -le 3 ]]; then
	echo "[GABBI] WARNING: Unable to compute ancestral sequences on fewer than 4 taxa. Keeping clean sequences for temporary probes."
	parallel --plus -j "$THREADS" '
	    sed -E -e "/>/ s/>/>{/...}|/g" -e "/>/! s/[-N]//g" -e "s/\.temp\.mafft//g" {}
	' ::: $(find alignments/TMP.phylomera.${PRE}/FASTA/ -type f -name "*dropped0") \
	> ${PRE}.temp.anc.loci.fasta

	echo "[GABBI] Skipping steps 4.3 and 4.4"
	checkpoint_mark "step4.3_ancestral_seqs" > /dev/null
	checkpoint_mark "step4.4_temp_loci" > /dev/null
        silence=true
    else
	if checkpoint_fail_exists "step4.3_ancestral_seqs"; then
	    rm -rf ancestral_seqs
	fi
        echo "[GABBI] Step 4.3: Computing ancestral sequences of temporary loci..."
	mkdir -p ancestral_seqs

	ancestral_seq() {
	    local phylip="$1"
	    local base="${phylip%.phylip}"
	    local shr
	    shr=$(basename "$base")
	    local model=$(grep "Best" $base.nogenepart.MFP.iqtree |cut -f2 -d':'|tr -d ' ')

	    iqtree3 -s $phylip -m $model -asr \
		-te $base.nogenepart.MFP.treefile \
		-pre ancestral_seqs/$shr.anc
	}
	export -f ancestral_seq

	parallel --plus -j "$THREADS" ancestral_seq {} ::: $(find alignments/GENETREES/ -type f -name "*.phylip") \
	    || checkpoint_fail "step4.3_ancestral_seqs"

	find ancestral_seqs/ -name "*.gz" -delete

	checkpoint_mark "step4.3_ancestral_seqs"
    fi
fi

# ---------------------------------------------------------------------------
# Step 4.4 — Create temporary loci file
# ---------------------------------------------------------------------------
if checkpoint_done "step4.4_temp_loci"; then
    if ! $silence;then
        echo "[GABBI] Skipping step 4.4 — checkpoint found."
    fi
else
    if checkpoint_fail_exists "step4.4_temp_loci"; then
        rm ${PRE}.temp.anc.loci.fasta
    fi

    echo "[GABBI] Step 4.4: Merging extent and ancestral sequences of temporary loci..."

    parallel -j "$THREADS" convert_state {} ::: $(find ancestral_seqs/ -type f -name "*.state") \
        || checkpoint_fail "step4.4_temp_loci"

    # Merge extent and ancestral sequences
    parallel --plus -j "$THREADS" '
        cat <(sed -E -e "/>/ s/>/>{/...}|/g" -e "/>/! s/[-N]//g" alignments/TMP.phylomera.*/FASTA/{/...}.*dropped0 ) \
            {} \
            <(echo)
    ' ::: $(find ancestral_seqs/ -type f -name "*.ok.fasta") \
    > ${PRE}.temp.anc.loci.fasta

    echo "[GABBI] Number of temporary probes: $(grep -c ">" ${PRE}.temp.anc.loci.fasta)"

    checkpoint_mark "step4.4_temp_loci"
fi

echo "[GABBI] ============================================================"
echo "[GABBI] PHASE 4: Completed at $(date '+%Y/%m/%d %H:%M:%S')"
echo "[GABBI] Temporary probe set (with ancestral sequences):    ${OUT:-GABBI_out}/04_shr_extraction/${PRE}.anc.loci.fasta"
echo "[GABBI] Number of temporary targeted loci:                 $(ls alignments/TMP.phylomera.${PRE}/FASTA/*dropped0|wc -l)"
echo "[GABBI] Number of temporary probes:                        $(grep -c ">" ${PRE}.temp.anc.loci.fasta)"
echo "[GABBI] ============================================================"


cd "$GABBI_WORKDIR"
