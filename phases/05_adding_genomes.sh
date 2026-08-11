#!/usr/bin/env bash
# Phase 5 — Temporary probe design and mapping to additional genomes
# Inputs:  04_shr_extraction/${PRE}.temp.anc.loci.fasta, $ADD_GENOMES
# Outputs: 05_adding_genomes/mapping/${PRE}-genome-fasta/
#          05_adding_genomes/consensus_loci/${PRE}.${SHR_THRESHOLD}.phyluce.loci.cons.fasta

source /opt/gabbi/utils/gabbi_functions.sh
source /opt/gabbi/utils/checkpoint.sh

mkdir -p 05_adding_genomes
cd 05_adding_genomes

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate /opt/miniconda3/envs/phyluce

# ---------------------------------------------------------------------------
# Step 5.1 — Generate tiled temporary probes
# ---------------------------------------------------------------------------
if checkpoint_done "step5.1_temp_probes"; then
    echo "[GABBI] Skipping step 5.1 — checkpoint found."
else
    if checkpoint_fail_exists "step5.1_temp_probes"; then
        rm -rf temp_probes
    fi

    echo "[GABBI] Step 5.1: Generating tiled temporary probes..."

    mkdir -p temp_probes

    python3 /opt/gabbi/scripts/get_tiled_probes.py \
        "$GABBI_WORKDIR"/04_shr_extraction/${PRE}.temp.anc.loci.fasta \
        temp_probes/${PRE}.temp.anc.probes.fasta

    # Replace gap characters with ambiguous bases (required by Phyluce, in case any remains)
    sed -i "/>/! s/-/N/g" temp_probes/${PRE}.temp.anc.probes.fasta

    echo "[GABBI] Temporary probes generated: $(grep -c ">" temp_probes/${PRE}.temp.anc.probes.fasta)"
    checkpoint_mark "step5.1_temp_probes"
fi

# ---------------------------------------------------------------------------
# Step 5.2 — Prepare additional genome assemblies
# ---------------------------------------------------------------------------
if checkpoint_done "step5.2_add_genomes_prep"; then
    echo "[GABBI] Skipping step 5.2 — checkpoint found."
else
    if checkpoint_fail_exists "step5.2_add_genomes_prep"; then
        rm -rf add_genomes_2bit
    fi  

    echo "[GABBI] Step 5.2: Preparing additional genomes..."

    mkdir -p add_genomes_2bit

    # Get 2bit files for each additional genome
    for genome in "$ADD_GENOMES"/*;do
        base=$(basename "$genome")
        genome_fasta=$(find "$genome" \( -name "*.fasta" -o -name "*.fas" -o -name "*.fna" \) -type f)
        verbose "Running faToTwoBit on genome_fasta=$genome_fasta"
        mkdir -p add_genomes_2bit/${base}
        faToTwoBit "$genome_fasta" "add_genomes_2bit/${base}/${base}.2bit" \
            || checkpoint_fail "step5.2_add_genomes_prep"
    done

    # Get 2bit files for each chromosomal genome
    echo "[GABBI] Preparing chromosomal genomes..."
    if [[ -n "$CHR_GENOMES" ]]; then
        for chr in "$CHR_GENOMES"/*; do
            genome=$(basename "$chr")
            mkdir -p add_genomes_2bit/${genome}
            genome_fasta=$(find "$CHR_GENOMES/$genome" \( -name "*.fasta" -o -name "*.fas" -o -name "*.fna" \) -type f)
            verbose "Running faToTwoBit on genome_fasta=$genome_fasta"
            faToTwoBit "$genome_fasta" "add_genomes_2bit/${genome}/${genome}.2bit" \
                || checkpoint_fail "step5.2_add_genomes_prep"
        done
    else
        for genome in $(halStats "$HAL"|awk -F"[ ,]" '$3==0{ print $1 }'); do
            mkdir -p add_genomes_2bit/${genome}
            echo "[GABBI] Getting $genome genome from $HAL..."
            hal2fasta "$HAL" "$genome" > "add_genomes_2bit/${genome}/${genome}.fasta" \
                || checkpoint_fail "step5.2_add_genomes_prep"
            verbose "Running faToTwoBit on genome_fasta=$genome_fasta"
            faToTwoBit "add_genomes_2bit/${genome}/${genome}.fasta" "add_genomes_2bit/${genome}/${genome}.2bit" \
                || checkpoint_fail "step5.2_add_genomes_prep"
            rm "add_genomes_2bit/${genome}/${genome}.fasta"
        done
    fi

    checkpoint_mark "step5.2_add_genomes_prep"
fi

# ---------------------------------------------------------------------------
# Step 5.3 — Map probes to genomes with LASTZ via Phyluce
# ---------------------------------------------------------------------------
if checkpoint_done "step5.3_lastz"; then
    echo "[GABBI] Skipping step 5.3 — checkpoint found."
else
    if checkpoint_fail_exists "step5.3_lastz"; then
        rm -rf mapping
    fi  
    echo "[GABBI] Step 5.3: Mapping temporary probes to genomes with LASTZ..."

    mkdir -p mapping/phyluce_logs

    echo "[scaffolds]" > mapping/assembled_genomes.conf
    for i in $(realpath add_genomes_2bit/*/*2bit); do
        base="${i##*/}"
        echo "${base%.*}:${i}" >> mapping/assembled_genomes.conf
    done

    verbose "Genome config:"
    verbose "$(cat mapping/assembled_genomes.conf)"

    ulimit -n 8192

    phyluce_probe_run_multiple_lastzs_sqlite \
        --probefile temp_probes/${PRE}.temp.anc.probes.fasta \
        --scaffoldlist $(cut -f1 -d':' mapping/assembled_genomes.conf | tail -n +2) \
        --genome-base-path add_genomes_2bit \
        --identity 50 \
        --cores "$THREADS" \
        --db mapping/${PRE}.temp.anc.probes.sqlite \
        --output mapping/${PRE}-genome-lastz \
         --log-path mapping/phyluce_logs \
        || checkpoint_fail "step5.3_lastz"

    phyluce_probe_slice_sequence_from_genomes \
        --contig_orient \
        --conf mapping/assembled_genomes.conf \
        --lastz mapping/${PRE}-genome-lastz/ \
        --probes 180 \
        --name-pattern "${PRE}.temp.anc.probes.fasta_v_{}.lastz.clean" \
        --output mapping/${PRE}-genome-fasta \
        --log-path mapping/phyluce_logs \
        || checkpoint_fail "step5.3_lastz"

    checkpoint_mark "step5.3_lastz"
fi

# ---------------------------------------------------------------------------
# Step 5.4 — Build multi-FASTA table and identify shared SHR loci
# ---------------------------------------------------------------------------
if checkpoint_done "step5.4_multifasta_table"; then
    echo "[GABBI] Skipping step 5.4 — checkpoint found."
else
    if checkpoint_fail_exists "step5.4_multifasta_table"; then
        rm -rf multifasta_table
    fi  

    echo "[GABBI] Step 5.4: Building multi-FASTA table..."

    mkdir -p multifasta_table/phyluce_logs

    # Find the genome capturing the most SHR loci
    egrep -o "INFO - .* written" mapping/phyluce_logs/phyluce_probe_slice_sequence_from_genomes.log \
        | awk -F'[: ]' '{print $3 "\t" $17}' \
        | sort -n -k2 \
        | tail -n 1 \
        | cut -f1 \
        > multifasta_table/base_taxon.txt

    verbose "Base taxon = $(cat multifasta_table/base_taxon.txt)"

    # Use this genome to generate the sqlite database
    phyluce_probe_get_multi_fasta_table \
        --fastas mapping/${PRE}-genome-fasta \
        --output multifasta_table/${PRE}.multifastas.sqlite \
        --base-taxon "$(cat multifasta_table/base_taxon.txt)" \
        --log-path multifasta_table/phyluce_logs \
        || checkpoint_fail "step5.4_multifasta_table"

    phyluce_probe_query_multi_fasta_table \
        --db multifasta_table/${PRE}.multifastas.sqlite \
        --base-taxon "$(cat multifasta_table/base_taxon.txt)" \
        --output multifasta_table/${PRE}.conf \
        --log-path multifasta_table/phyluce_logs \
        || checkpoint_fail "step5.4_multifasta_table"

    egrep -o "Loci shared.*" multifasta_table/phyluce_logs/phyluce_probe_query_multi_fasta_table.log \
        | sed -e "s/,//g" -e "s/\.0//g" \
        > multifasta_table/${PRE}.table

    SHR_THRESHOLD_TAXA=$(tail -n +2 multifasta_table/${PRE}.table|awk -v thr="$SHR_THRESHOLD" 'END { print NR * ( thr / 100 ) }')
    SHR_THRESHOLD_NUMBER=$(awk -v thr_taxa="$SHR_THRESHOLD_TAXA" '$4 >= thr_taxa { print $6 }' multifasta_table/${PRE}.table | head -n 1)

    echo "[GABBI] SHR loci shared by >= ${SHR_THRESHOLD}% of taxa: $SHR_THRESHOLD_NUMBER"
    echo "[GABBI] Note: You can change this threshold with --shr-threshold based on 05_adding_genomes/multifasta_table/${PRE}.table"
    echo "[GABBI] and restart GABBI with --restart 5.5"
    sleep 10

    checkpoint_mark "step5.4_multifasta_table"
fi

if checkpoint_done "step5.5_final_phyluce_probes"; then
    echo "[GABBI] Skipping step 5.5 — checkpoint found."
else
    if checkpoint_fail_exists "step5.5_final_phyluce_probes"; then
        rm -rf final_phyluce_probes
    fi  

    mkdir -p final_phyluce_probes/phyluce_logs
    SHR_THRESHOLD_TAXA=$(tail -n +2 multifasta_table/${PRE}.table|awk -v thr="$SHR_THRESHOLD" 'END { print int( NR * ( thr / 100 ) ) }')
    SHR_THRESHOLD_NUMBER=$(awk -v thr_taxa="$SHR_THRESHOLD_TAXA" '$4 >= thr_taxa { print $6 }' multifasta_table/${PRE}.table | head -n 1)

    echo "[GABBI] Step 5.5: Getting final phyluce probes targeting $SHR_THRESHOLD_NUMBER loci shared by $SHR_THRESHOLD_TAXA taxa..."

    phyluce_probe_query_multi_fasta_table \
        --db multifasta_table/${PRE}.multifastas.sqlite \
        --base-taxon $(cat multifasta_table/base_taxon.txt) \
        --output final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.conf \
        --specific-counts $SHR_THRESHOLD_TAXA \
        --log-path final_phyluce_probes/phyluce_logs \
        || checkpoint_fail "step5.5_final_phyluce_probes"

    phyluce_probe_get_tiled_probe_from_multiple_inputs \
        --fastas mapping/${PRE}-genome-fasta \
        --multi-fasta-output final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.conf \
        --probe-prefix "shr-" --designer GABBI --design ${PRE} \
        --tiling-density "$FP_TILING_DENSITY" \
        --overlap middle \
        --masking "$FP_MASKING" \
        --remove-gc --two-probes \
        --output final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.probe_list.fasta \
        --log-path final_phyluce_probes/phyluce_logs \
        || checkpoint_fail "step5.5_final_phyluce_probes"

    # Clean probe set
    phyluce_probe_easy_lastz \
        --target final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.probe_list.fasta \
        --query final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.probe_list.fasta \
        --identity 50 --coverage 50 \
        --output final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.probe_list.to_self.lastz \
        --log-path final_phyluce_probes/phyluce_logs \
        || checkpoint_fail "step5.5_final_phyluce_probes"

    # probe prefix is left as "uce-" because previous attempts to modify probe prefix have no effect
    phyluce_probe_remove_duplicate_hits_from_probes_using_lastz \
        --fasta final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.probe_list.fasta \
        --lastz final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.probe_list.to_self.lastz \
        --probe-prefix "uce-" \
        --log-path final_phyluce_probes/phyluce_logs \
        || checkpoint_fail "step5.5_final_phyluce_probes"

    final_probes=$(egrep -o "Kept [0-9]+" final_phyluce_probes/phyluce_logs/phyluce_probe_remove_duplicate_hits_from_probes_using_lastz.log|awk '{print $2}')
    echo "[GABBI] Final number of probes using phyluce: $final_probes"

    checkpoint_mark "step5.5_final_phyluce_probes"

fi

if checkpoint_done "step5.6_consensus_loci"; then
    echo "[GABBI] Skipping step 5.6 — checkpoint found."

else
    if checkpoint_fail_exists "step5.6_consensus_loci"; then
        rm -rf consensus_loci
    fi  

    echo "[GABBI] Reconstructing consensus SHR from probes using phyluce..."

    mkdir -p consensus_loci/phyluce_logs

    phyluce_probe_reconstruct_uce_from_probe \
        --input final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.probe_list-DUPE-SCREENED.fasta \
        --output consensus_loci/${PRE}.${SHR_THRESHOLD}.phyluce.loci.cons.fasta \
        --log-path consensus_loci/phyluce_logs \
        || checkpoint_fail "step5.6_consensus_loci"

    # Renaming UCE into SHR didn't work so far, so apply these changes on this final file
    sed -i "s/uce-/shr-/g" consensus_loci/${PRE}.${SHR_THRESHOLD}.phyluce.loci.cons.fasta

    echo "[GABBI] Final number of loci: $(grep -c ">" consensus_loci/${PRE}.${SHR_THRESHOLD}.phyluce.loci.cons.fasta)"
    checkpoint_mark "step5.6_consensus_loci"
fi

SHR_THRESHOLD_TAXA=$(tail -n +2 multifasta_table/${PRE}.table|awk -v thr="$SHR_THRESHOLD" 'END { print NR * ( thr / 100 ) }')

echo "[GABBI] ============================================================"
echo "[GABBI] PHASE 5: Completed at $(date '+%Y/%m/%d %H:%M:%S')"
echo "[GABBI] Minimum number of taxa to keep a locus:            $SHR_THRESHOLD_TAXA (${SHR_THRESHOLD}%)"
echo "[GABBI] Final probe set (without ancestral sequences):     ${OUT:-GABBI_out}/05_adding_genomes/final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.probe_list-DUPE-SCREENED.fasta"
echo "[GABBI] Consensus sequences of targeted loci:              ${OUT:-GABBI_out}/05_adding_genomes/consensus_loci/${PRE}.${SHR_THRESHOLD}.phyluce.loci.cons.fasta"
echo "[GABBI] Number of targeted loci:                           $(grep -c ">" consensus_loci/${PRE}.${SHR_THRESHOLD}.phyluce.loci.cons.fasta)"
echo "[GABBI] Number of probes:                                  $(grep -c ">" final_phyluce_probes/${PRE}.${SHR_THRESHOLD}.probe_list-DUPE-SCREENED.fasta)"
echo "[GABBI] Note: You can change this threshold with --shr-threshold based on 05_adding_genomes/multifasta_table/${PRE}.table"
echo "[GABBI]       and restart GABBI with --restart 5.5"
echo "[GABBI] ============================================================"

cd "$GABBI_WORKDIR"
