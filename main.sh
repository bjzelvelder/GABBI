#!/usr/bin/env bash
# called by %runscript with exec /opt/gabbi/main.sh $@

set -euo pipefail

source /opt/gabbi/utils/gabbi_functions.sh

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in

	# -- Required arguments ---------------------------------------------
	--chr-genomes)          export CHR_GENOMES="$2";           shift 2 ;;
	--guide-tree)           export GUIDE_TREE="$2";            shift 2 ;;
        --hal)                  export HAL="$2";                   shift 2 ;;
	--add-genomes)          export ADD_GENOMES="$2";           shift 2 ;;

	# -- Pipeline control -----------------------------------------------
        --prefix)               export PRE="$2";                   shift 2 ;;
        -o|--out-dir)           export OUT="$2";                   shift 2 ;;
	-t|--threads)           export THREADS="$2";               shift 2 ;;
        --restart)              export RESTART="$2";               shift 2 ;;
	--stop-before)          export STOP_BEFORE="$2";           shift 2 ;;
        --verbose)              export GABBI_VERBOSE=1;              shift ;;

	# -- Genome alignment options ---------------------------------------
        --cactus-maxDisk)       export CACTUS_MAXDISK="$2";        shift 2 ;;
        --cactus-maxCores)      export CACTUS_MAXCORES="$2";       shift 2 ;;
	--cactus-maxMemory)     export CACTUS_MAXMEM="$2";         shift 2 ;;
        --maf-references)       export MAF_REFS="$2";              shift 2 ;;
	--block-size)           export BLOCK_SIZE="$2";            shift 2 ;;
	--block-length)         export BLOCK_LENGTH="$2";          shift 2 ;;

	# -- Temporary probe options ----------------------------------------
	--cross-blast-ev)       export CROSS_BLAST_EV="$2";        shift 2 ;;
	--cross-blast-ws)       export CROSS_BLAST_WS="$2";        shift 2 ;;
	--cross-blast-qc)       export CROSS_BLAST_QC="$2";        shift 2 ;;
	--temp-tax-threshold)   export TEMP_TAX_THRESHOLD="$2";    shift 2 ;;
	--temp-allow-dupes)     export TEMP_ALLOW_DUPES="$2";      shift 2 ;;

	# -- Final probe design options -------------------------------------
	--shr-threshold)        export SHR_THRESHOLD="$2";         shift 2 ;;
        --final-probes-tiling)  export FP_TILING_DENSITY="$2";     shift 2 ;;
        --final-probes-masking) export FP_MASKING="$2";            shift 2 ;;

	*) echo "[GABBI] ERROR: Unknown argument: $1" >&2; exit 1 ;;

    esac
done

# ---------------------------------------------------------------------------
# Validate required arguments
# ---------------------------------------------------------------------------

[[ -n "${OUT:-}" && ! -w "$(dirname "${OUT:-}")" ]] && echo "[GABBI] ERROR: --out-dir '${OUT}' parent directory is not writable." >&2 && exit 1
[[ -d "${OUT:-GABBI_out}" ]] && echo "[GABBI] WARNING: output directory '${OUT:-GABBI_out}' already exists. Checking existing checkpoints." >&2
GABBI_WORKDIR="${OUT:-GABBI_out}"

[[ -z "${HAL:-}" && -z "${CHR_GENOMES:-}" ]] && echo "[GABBI] ERROR: --chr-genomes OR --hal is required. Please provide a path to the directory containing chromosome-level genomes with --chr-genomes or a HAL cactus alignment file with --hal." && exit 1

[[ -n "${CHR_GENOMES:-}" && ! -d "$CHR_GENOMES" ]] && echo "[GABBI] ERROR: --chr-genomes '${CHR_GENOMES}' does not exist or is not a directory." >&2 && exit 1
[[ -n "${CHR_GENOMES:-}" && -n "$(find "${CHR_GENOMES}" -maxdepth 1 -mindepth 1 -type d -name '*.*')" ]] && echo "[GABBI] ERROR: Chromosome-level genome folders should not contain any \".\"." >&2 && exit 1

[[ -n "${HAL:-}" && ! -f "${HAL:-}" ]] && echo "[GABBI] ERROR: --hal '${HAL}' does not exist or is not a file." >&2 && exit 1

if [[ -z "${HAL:-}" ]] ;then 
    export N_CHR_TAXA=$(ls -d "$CHR_GENOMES"/*/ 2>/dev/null |wc -l)
    if [[ -n "${CHR_GENOMES:-}" && "$N_CHR_TAXA" -eq 0 ]];then
        echo "[GABBI] ERROR: --chr-genomes '${CHR_GENOMES}' contains no subdirectory." >&2
        exit 1
    fi
else
    export N_CHR_TAXA=$(halStats "$HAL"|awk -F"[ ,]" '$3==0{ print $1 }'|wc -l)
fi

verbose "N_CHR_TAXA = $N_CHR_TAXA"

[[ -z "${GUIDE_TREE:-}" && -z "${HAL:-}" ]] && echo "[GABBI] ERROR: --guide-tree is required. Please provide a path to a Newick-format species tree of chromosome-level genomes to provide a guide tree for Cactus." >&2 && exit 1
[[ -n "${GUIDE_TREE:-}" && ! -f "${GUIDE_TREE:-}" ]] && echo "[GABBI] ERROR: --guide-tree '${GUIDE_TREE}' does not exist or is not a file." >&2 && exit 1

[[ -z "${ADD_GENOMES:-}" ]] && echo "[GABBI] WARNING: No additional genomes provided. GABBI will stop before phase '05_adding_genomes'. To finish running GABBI with additional genomes, rerun the same command with --add-genomes option." >&2 && STOP_BEFORE="05"
[[ -n "${ADD_GENOMES:-}" && ! -d "$ADD_GENOMES" ]] && echo "[GABBI] ERROR: --add-genomes '${ADD_GENOMES}' does not exist or is not a directory." >&2 && exit 1
export N_ADD_TAXA=$(ls -d "$ADD_GENOMES"/*/ 2>/dev/null |wc -l)

if [[ -n "${ADD_GENOMES:-}" && "$N_ADD_TAXA" -eq 0 ]];then
    echo "[GABBI] ERROR: --add-genomes '${ADD_GENOMES}' contains no subdirectory." >&2
    exit 1
fi
[[ -n "${ADD_GENOMES:-}" && -n "$(find "${ADD_GENOMES}" -maxdepth 1 -mindepth 1 -type d -name '*.*')" ]] && echo "[GABBI] ERROR: Additional genome folders should not contain any \".\"." >&2 && exit 1

verbose "N_ADD_TAXA = $N_ADD_TAXA"

# -- Genome alignment options -----------------------------------------------
[[ -n "${THREADS:-}" && ! "$THREADS" =~ ^[1-9][0-9]*$ ]] && echo "[GABBI] ERROR: --threads '${THREADS}' must be a positive integer." >&2 && exit 1

[[ -n "${CACTUS_MAXDISK:-}" && ! "$CACTUS_MAXDISK" =~ ^[1-9][0-9]*[MG]$ ]] && echo "[GABBI] ERROR: --cactus-maxDisk '${CACTUS_MAXDISK}' must be a positive integer followed by 'M or G' (Mega or Gigabytes)." >&2 && exit 1
[[ -n "${CACTUS_MAXCORES:-}" && ! "$CACTUS_MAXCORES" =~ ^[1-9][0-9]*$ ]] && echo "[GABBI] ERROR: --cactus-maxCores '${CACTUS_MAXCORES}' must be a positive integer." >&2 && exit 1
[[ -n "${CACTUS_MAXMEM:-}" && ! "$CACTUS_MAXMEM" =~ ^[1-9][0-9]*[MG]$ ]] && echo "[GABBI] ERROR: --cactus-maxMemory '${CACTUS_MAXMEM}' must be a positive integer followed by 'M or G' (Mega or Gigabytes)." >&2 && exit 1

[[ -n "${MAF_REFS:-}" && ! -f "${MAF_REFS:-}" ]] && echo "[GABBI] ERROR: --maf-references '${GUIDE_TREE}' does not exist or is not a file." >&2 && exit 1
if [[ -n "${MAF_REFS:-}" ]]; then
    if [[ -z "${HAL:-}" ]] ;then 
	if [[ $(grep -o -f "$MAF_REFS" <(ls -d "$CHR_GENOMES"/*/ 2>/dev/null) | wc -l ) -ne $(cat "$MAF_REFS"|wc -l) ]]; then
	    echo "[GABBI] ERROR: Some genomes listed in '${MAF_REFS:-}' cannot be found in chromosome-level genomes directory." >&2
	    verbose "Matching genomes: $(grep -o -f "$MAF_REFS" <(ls -d "$CHR_GENOMES"/*/ 2>/dev/null))"
	    exit 1
	fi
    else
	if [[ $(grep -o -f "$MAF_REFS" <(halStats --genomes "$HAL") | wc -l ) -ne $(cat "$MAF_REFS"|wc -l) ]]; then
	    echo "[GABBI] ERROR: Some genomes listed in '${MAF_REFS:-}' cannot be found in '${HAL}' genome alignment file." >&2
	    verbose "Matching genomes: $(grep -o -f "$MAF_REFS" <(halStats --genomes "$HAL"))"
	    exit 1
	fi
    fi
fi

[[ -n "${BLOCK_SIZE:-}" && ! "$BLOCK_SIZE" =~ ^[1-9][0-9]*$ ]] && echo "[GABBI] ERROR: --block-size '${BLOCK_SIZE}' must be a positive integer (minimum number of taxa in a block of alignment)." >&2 && exit 1
[[ -n "${BLOCK_LENGTH:-}" && ! "$BLOCK_LENGTH" =~ ^[1-9][0-9]*$ ]] && echo "[GABBI] ERROR: --block-length '${BLOCK_LENGTH}' must be a positive integer (minimum alignment block length in nucleotides)." >&2 && exit 1

# -- Temporary probe options -------------------------------------------------
[[ -n "${CROSS_BLAST_EV:-}" && ! "$CROSS_BLAST_EV" =~ ^[0-9]+([.][0-9]+)?(E|e)-?[0-9]+$ ]] && echo "[GABBI] ERROR: --cross-blast-ev '${CROSS_BLAST_EV}' must be a valid scientific notation e-value (e.g. 1E-6)." >&2 && exit 1
[[ -n "${CROSS_BLAST_WS:-}" && ! "$CROSS_BLAST_WS" =~ ^[1-9][0-9]*$ ]] && echo "[GABBI] ERROR: --cross-blast-ws '${CROSS_BLAST_WS}' must be a positive integer." >&2 && exit 1
[[ -n "${CROSS_BLAST_QC:-}" && ! "$CROSS_BLAST_QC" =~ ^([1-9][0-9]?|100)$ ]] && echo "[GABBI] ERROR: --cross-blast-qc '${CROSS_BLAST_QC}' must be an integer between 1 and 100 (percentage)." >&2 && exit 1

[[ -n "${TEMP_TAX_THRESHOLD:-}" && ! "$TEMP_TAX_THRESHOLD" =~ ^([1-9][0-9]?|100)$ ]] && echo "[GABBI] ERROR: --temp-tax-threshold '${TEMP_TAX_THRESHOLD}' must be an integer between 1 and 100 (percentage)." >&2 && exit 1
[[ -n "${TEMP_ALLOW_DUPES:-}" && "$TEMP_ALLOW_DUPES" > $N_CHR_TAXA && ! "$TEMP_ALLOW_DUPES" =~ ^[0-9]+$ ]] && echo "[GABBI] ERROR: --temp-allow-dupes '${TEMP_ALLOW_DUPES}' must be a non-negative integer." >&2 && exit 1

# -- Final probe design options ---------------------------------------------
[[ -n "${SHR_THRESHOLD:-}" && ! "$SHR_THRESHOLD" =~ ^([1-9][0-9]?|100)$ ]] && echo "[GABBI] ERROR: --shr-threshold '${SHR_THRESHOLD}' must be an integer between 1 and 100 (percentage)." >&2 && exit 1
[[ -n "${FP_TILING_DENSITY:-}" && ! "$TEMP_ALLOW_DUPES" =~ ^[1-9][0-9]*$ ]] && echo "[GABBI] ERROR: --final-probes-tiling '${FP_TILING_DENSITY}' must be a positive integer." >&2 && exit 1
[[ -n "${FP_MASKING:-}" && ! "$FP_MASKING" =~ ^0(\.[0-9]+)?$|^1(\.0+)?$ ]] && echo "[GABBI] ERROR: --final-probes-masking '${FP_MASKING}' must be a number between 0 and 1." >&2 && exit 1

# -- General options --------------------------------------------------------
[[ -n "${PRE:-}" && "${PRE:-}" =~ / ]] && echo "[GABBI] ERROR: --prefix should not contain any \"/\"." >&2 && exit 1

# Validate --stop-before against the list of known phase names
VALID_PHASES=(
    "1" "01" "01_cactus_alignment" "cactus_alignment" "cactus"
    "2" "02" "02_conserved_loci" "conserved_loci"
    "3" "03" "03_cross_blast" "cross_blast" "blast"
    "4" "04" "04_shr_extraction" "shr_extraction" "extraction"
    "5" "05" "05_adding_genomes" "adding_genomes" "add_genomes"
    "6" "06" "06_final_targeted_loci" "final_targeted_loci"
)

if [[ -n "${STOP_BEFORE:-}" ]]; then
    valid=0
    for phase in "${VALID_PHASES[@]}"; do
        [[ "$STOP_BEFORE" == "$phase" ]] && valid=1 && break
    done
    if [[ $valid -eq 0 ]]; then
        echo "[GABBI] ERROR: --stop-before '${STOP_BEFORE}' is not a recognised phase name." >&2
        echo "[GABBI] Valid values are: $(echo "${VALID_PHASES[@]}"|sed "s/ /, /g")" >&2
        exit 1
    fi
fi

export STOP_BEFORE="${STOP_BEFORE:-}"   # empty by default = run all phases

# Validate RESTART from a list of valid steps
VALID_STEPS=(
    "1.1" "1.1_cactus_alignment"
    "2.1" "2.1_hal2maf" "2.2" "2.2_split_maf" "2.3" "2.3_maffilter" "2.4" "2.4_phastcons" "2.5" "2.5_conserved_loci"
    "3.1" "3.1_cross_blast" "3.2" "3.2_shr_clustering"
    "4.1" "4.1_shr_extraction" "4.2" "4.2_alignments" "4.3" "4.3_ancestral_seqs" "4.4" "4.4_temp_loci"
    "5.1" "5.1_temp_probes" "5.2" "5.2_add_genomes_prep" "5.3" "5.3_lastz" "5.4" "5.4_multifasta_table"
    "5.5" "5.5_final_phyluce_probes" "5.6" "5.6_consensus_loci"
    "6.1" "6.1_extract_targeted_loci" "6.2" "6.2_align_targeted_loci" "6.3" "6.3_final_ancestral_seqs" "6.4" "6.4_final_targeted_loci"
)
if [[ -n "${RESTART:-}" ]]; then
    valid=0
    for step in "${VALID_STEPS[@]}"; do
        [[ "$RESTART" == "$step" ]] && valid=1 && break
    done
    if [[ $valid -eq 0 ]]; then
        echo "[GABBI] ERROR: --restart '${RESTART}' is not a recognised step name." >&2
        echo "[GABBI] Valid values are: $(echo "${VALID_STEPS[@]}"|sed "s/ /, /g")" >&2
        exit 1
    fi
fi

# Resolve all input paths to absolute before any cd can invalidate them
[[ -n "${CHR_GENOMES:-}" ]] && export CHR_GENOMES=$(realpath "$CHR_GENOMES")
[[ -n "${GUIDE_TREE:-}" ]] && export GUIDE_TREE=$(realpath "$GUIDE_TREE")
[[ -n "${HAL:-}" ]] && export HAL=$(realpath "$HAL")
[[ -n "${MAF_REFS:-}" ]] && export MAF_REFS=$(realpath "$MAF_REFS")
[[ -n "${ADD_GENOMES:-}" ]] && export ADD_GENOMES=$(realpath "$ADD_GENOMES")

verbose "CHR_GENOMES=${CHR_GENOMES:-}"
verbose "GUIDE_TREE=${GUIDE_TREE:-}"
verbose "HAL=${HAL:-}"
verbose "MAF_REFS=${MAF_REFS:-}"
verbose "ADD_GENOMES=${ADD_GENOMES:-}"

export GABBI_WORKDIR=$(realpath "$GABBI_WORKDIR")

# Get core number to run multiple parallel jobs based on the number of taxa
CORES_PER_JOB=$(( THREADS / N_CHR_TAXA ))
export CORES_PER_JOB=$(( CORES_PER_JOB < 1 ? 1 : CORES_PER_JOB ))
export PARALLEL_JOBS=$(( THREADS / CORES_PER_JOB ))

verbose "CORES_PER_JOB=$CORES_PER_JOB"
verbose "PARALLEL_JOBS=$PARALLEL_JOBS"

# ---------------------------------------------------------------------------
# Main GABBI script that sources utilities, sets arguments and default parameters and calls phases in order
# ---------------------------------------------------------------------------

GABBI_ROOT="/opt/gabbi"

# Load utilities and defaults
source "${GABBI_ROOT}/utils/checkpoint.sh"
source "${GABBI_ROOT}/utils/defaults.sh"


####################################################################


echo "[GABBI] Pipeline starting at $(date '+%Y/%m/%d %H:%M:%S')"

mkdir -p "$GABBI_WORKDIR"
cd "$GABBI_WORKDIR"
echo "[GABBI] Working directory: $GABBI_WORKDIR"

# Remove checkpoints after specified restart step
verbose "RESTART=$RESTART"
checkpoint_restart

# Run phases sequentially
for phase in "${GABBI_ROOT}/phases/"*.sh; do
    phase_name=$(basename "$phase" .sh)
    verbose "phase=$phase_name"
    if [[ -n "${STOP_BEFORE:-}" && "$phase_name" =~ "${STOP_BEFORE:-}" ]]; then
        echo "[GABBI] Stopping before phase '${STOP_BEFORE}' as requested."
        break
    fi
    if [[ -n "${HAL:-}" && "$phase_name" =~ "01_cactus_alignment" ]]; then
        echo "[GABBI] Skipping Cactus whole-genome alignment and using '${HAL}' as input for phase 2."
        continue
    fi
    
    bash "$phase"
done

echo "[GABBI]"
echo "[GABBI] ============================================================"
echo "[GABBI] GABBI Pipeline: Completed at $(date '+%Y/%m/%d %H:%M:%S')"
echo "[GABBI] GABBI output: $GABBI_WORKDIR"
echo "[GABBI] ============================================================"
