#!/usr/bin/env bash
# Phase 6 — Get targeted loci and ancestral sequences from final phyluce probe set
# Inputs:  05_adding_genomes/mapping/${PRE}-genome-fasta/
#          05_adding_genomes/conserved_loci/${PRE}.${SHR_THRESHOLD}.phyluce.loci.cons.fasta
# Outputs: 06_final_targeted_loci/${PRE}.final.anc.loci.fasta

source /opt/gabbi/utils/gabbi_functions.sh
source /opt/gabbi/utils/checkpoint.sh

mkdir -p 06_final_targeted_loci
cd 06_final_targeted_loci

# ---------------------------------------------------------------------------
# Step 6.1 — Extract targeted loci obtained with phyluce
# ---------------------------------------------------------------------------
if checkpoint_done "step6.1_extract_targeted_loci"; then
    echo "[GABBI] Skipping step 6.1 — checkpoint found."
else
    if checkpoint_fail_exists "step6.1_extract_targeted_loci"; then
        rm -rf split_loci targeted_loci
    fi  

    echo "[GABBI] Step 6.1: Extracting targeted loci from lastz mapping..."

    mkdir -p split_loci targeted_loci

    genome_loci="$GABBI_WORKDIR/05_adding_genomes/mapping/${PRE}-genome-fasta/"

    # Regroup per-taxon sequences into per-locus FASTA files
    for genome in "$genome_loci"/*.fasta; do
        sp=${genome##*/}
        sp=${sp/.fasta/}
        seqkit split -i $genome --out-dir split_loci --by-id-prefix "${sp}." \
            || checkpoint_fail "step6.1_extract_targeted_loci"
        find split_loci/ -type f -name "${sp}.*" | xargs sed -E -i "s/>.*/>${sp}/g"
    done

    parallel -j "$THREADS" " \
        cat split_loci/*__{}__* > targeted_loci/{}.fasta
    " ::: $(tail -n +4 "$GABBI_WORKDIR"/05_adding_genomes/final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.conf)
    
	for i in $(find targeted_loci -type f -name "*.fasta"); do
        mv $i ${i/uce/shr}
	done
	
    checkpoint_mark "step6.1_extract_targeted_loci"
fi

# ---------------------------------------------------------------------------
# Step 6.2 — Align and clean targeted loci to compute genetrees
# ---------------------------------------------------------------------------
if checkpoint_done "step6.2_align_targeted_loci"; then
    echo "[GABBI] Skipping step 6.2 — checkpoint found."
else
    if checkpoint_fail_exists "step6.2_align_targeted_loci"; then
        rm -rf final_alignments
    fi  

    echo "[GABBI] Step 6.2: Aligning targeted loci and building gene trees..."

    # Remove variation held by only one taxon (fix the AMAS trim threshold to 2 taxa out of total taxa)
    N_TAXA=$(ls "$GABBI_WORKDIR"/05_adding_genomes/mapping/${PRE}-genome-fasta/|wc -l)
    trim=$(awk -v n="$N_TAXA" 's=200/n { print int(s) }' <(echo ))

    phylomera \
        --input targeted_loci \
        --output final_alignments \
        --prefix ${PRE} \
        --trim $trim \
        --drop 0 \
        --genetrees MFP \
        --nogenepart \
        --perc 0 \
        --threads "$THREADS" \
        --debug  \
        --continue \
        --config /opt/gabbi/config/phylomera.conf \
        || checkpoint_fail "step6.2_align_targeted_loci"

    checkpoint_mark "step6.2_align_targeted_loci"
fi

# ---------------------------------------------------------------------------
# Step 6.3 — Reconstruct ancestral sequences
# ---------------------------------------------------------------------------
if checkpoint_done "step6.3_final_ancestral_seqs"; then
    echo "[GABBI] Skipping step 6.3 — checkpoint found."
else
    if checkpoint_fail_exists "step6.3_final_ancestral_seqs"; then
        rm -rf ancestral_seqs
    fi
    
    echo "[GABBI] Step 6.3: Reconstructing ancestral sequences with IQ-TREE..."

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

    parallel --plus -j "$THREADS" ancestral_seq {} ::: $(find final_alignments/GENETREES/ -type f -name "*.phylip") \
	|| checkpoint_fail "step6.3_final_ancestral_seqs"

    rm -rf ancestral_seqs/*.gz

    checkpoint_mark "step6.3_final_ancestral_seqs"
fi

# ---------------------------------------------------------------------------
# Step 6.4 — Create final targeted loci file
# ---------------------------------------------------------------------------
if checkpoint_done "step6.4_final_targeted_loci"; then
    echo "[GABBI] Skipping step 6.4 — checkpoint found."
else
    if checkpoint_fail_exists "step6.4_final_targeted_loci"; then
        rm -rf final_consensus_loci ${PRE}.final.anc.loci.fasta
    fi

    echo "[GABBI] Step 6.4: Merging ancestral sequences with final targeted loci..."

    parallel -j "$THREADS" convert_state {} ::: $(find ancestral_seqs/ -type f -name "*.state") \
        || checkpoint_fail "step6.4_final_targeted_loci"

    parallel --plus -j "$THREADS" '
        cat <(sed -E -e "/>/ s/>/>{/...}|/g" -e "/>/! s/[-N]//g" final_alignments/TMP.phylomera.*/FASTA/{/...}.*dropped0 ) \
            {}
    ' ::: $(find ancestral_seqs/ -type f -name "*.ok.fasta") \
    > ${PRE}.final.anc.loci.fasta
	sed -i 's/\([^>]\)>/\1\n>/g' ${PRE}.final.anc.loci.fasta

    # Make a consensus for phylomera references
    mkdir -p final_consensus_loci
    for i in $(find final_alignments/TMP.phylomera.${PRE}/FASTA/ -type f -name "*dropped0"); do
        base=${i##*/}
        base=${base%%.*}
        cp $i final_consensus_loci/$base.mafft.fasta
    done || checkpoint_fail "step6.4_final_targeted_loci"

    cd final_consensus_loci/
        Rscript /opt/gabbi/scripts/make_consensus_from_mafft_v3.R all iupac 1 \
            || checkpoint_fail "step6.4_final_targeted_loci"
    cd ..
    cat $(find final_consensus_loci/ -type f -name "*.cons") > ${PRE}.final.anc.loci.cons.fasta
    sed -i "s/_mafft//g" ${PRE}.final.anc.loci.cons.fasta
	
    echo "[GABBI] Final number of targeted loci: $(grep -c ">" ${PRE}.final.anc.loci.cons.fasta)"

    checkpoint_mark "step6.4_final_targeted_loci"

fi

echo "[GABBI] ============================================================"
echo "[GABBI] PHASE 6: Completed at $(date '+%Y/%m/%d %H:%M:%S')"
echo "[GABBI] Final targeted SHR and ancestral sequences:       ${OUT:-GABBI_out}/06_final_targeted_loci/${PRE}.final.anc.loci.fasta"
echo "[GABBI] Consensus sequences of targeted SHR:              ${OUT:-GABBI_out}/06_final_targeted_loci/${PRE}.final.anc.loci.cons.fasta"
echo "[GABBI] Number of targeted loci:                          $(grep -c ">" ${PRE}.final.anc.loci.cons.fasta)"
echo "[GABBI] Number of targeted sequences:                     $(grep -c ">" ${PRE}.final.anc.loci.fasta)"
echo "[GABBI] ============================================================"

cd "$GABBI_WORKDIR"
