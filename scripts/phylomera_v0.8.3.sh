#!/bin/bash

VERSION="0.8.3"
# Fixed a bug that did not process markers when the reference was reversed in the alignment
# Added --continue/-y option to say "yes" to all interactive user questions (leaving -p 0 by default)
# Added --restart option to specify which step to restart without interactive prompt
# Added some statistics to the STDOUT

# Help
print_help() {
    echo "Phylomera v$VERSION"
    echo
    echo "Usage: $0 -i <fasta_folder> -o <output_folder> -pre <output_prefix> [options]"
    echo
    echo "DESCRIPTION:"
    echo "  This pipeline was developped to run phylogenomic analyses from a set of separate genomic markers."
    echo "  Important steps can be skipped or restarted interactively when running the same command again."
    echo "  ..."
    echo
    echo "REQUIRED ARGUMENTS:"
    echo "  --input, -i <fasta_folder>       Directory containing FASTA files to process (marker names must be the prefix of each file; e.g. uce_100.fasta)"
    echo "  --output, -o <output_folder>     Path to output directory"
    echo "  --prefix, -pre <output_prefix>   Prefix of output files (not a path)"
    echo
    echo "OPTIONS:"
    echo "  --refs, -r <path>                Multifasta file containing one reference sequence per marker (e.g. consensus sequence from tiled probes for UCE data)"
    echo "                                      References are used to split alignments into core and flanking regions that will be processed and partitionned separately"
    echo "                                      Marker names must match exactly with headers id (delimited by regex: \"^(\\S+)\\s?\")"
    echo "  --annotation, -a <auto|path>     Process coding markers more accurately with OMM_MACSE and make codon partitions"
    echo "                                      - auto: Define coding regions automatically based on OMM_MACSE output (! could lead to data loss)"
    echo "                                      - annotation file: only process core regions of coding markers (see example ${0%/*}/annotation_file.txt)"
    echo "  --ntax, -n <number>              Maximum number of taxa per file (default: $MAXTAX)"
    echo "  --trim <%>                       Trim alignment columns with less than this % characters (default: $GAP)"
    echo "  --drop <%>                       Drop sequences with less than this % characters (default: $DROP)"
    echo "  --perc, -p <%>                   Minimum percentage of taxa required in each alignment to build a supermatrix (e.g. 70, \"80 90\"; default: interactive unless --continue or --force option is specified)"
    echo "  --sptree, -s <model>             Ask Phylomera to build species trees with IQ-TREE 3 using the generated partition file (prefix.perc.nex) and the specified model (default: MFP+MERGE)"
    echo "  --genetrees, -g <model>          Ask Phylomera to build gene trees with IQ-TREE 3 using the specified model (default: GTR+I+F+G)"
    echo "  --nogenepart                     Ignore partitions when building gene trees (use this option with -g GTR+G8 for an accurate alpha parameter estimation)"
    echo "  --config, -c <conf>              Config file with path to external scripts (default: PATH) and internal scripts (default: ${0%/*}/)"
    echo "  --threads, -t <threads>          Number of threads to use (default: $THREADS)"
    echo "  --continue, -y                   Continue Phylomera without asking user (similar to answering yes to all interactive steps)"
    echo "  --restart <step>                 Restart Phylomera at the specified step (all, alignments, cleaning, stats, perc_spp, sptree, genetrees)"
    echo "  --verbose                        Print all subcommands to STDOUT (default: only pipeline steps)"
    echo "  --debug                          Do not remove temporary directory TMP.phylomera.prefix/"
    echo "  --help, -h                       Show this help message and exit"
    echo
    echo
    echo "Developed by Benjamin Zelvelder"
    echo "Last update: 16/03/2026"
    echo
    echo "If you use Phylomera, please cite:"
    echo "  -  The first article in which it is described:  "
    echo "       Zelvelder et al. 202X A new method based on genome alignments provides a highly resolutive target enrichment set for weevils (Coleoptera, Curculionoidea)"
    echo "  -  MAFFT: "
    echo "       Katoh K, Standley DM. 2013 MAFFT Multiple Sequence Alignment Software Version 7: Improvements in performance and usability. Molecular Biology and Evolution 30, 772–780. (doi:10.1093/molbev/mst010)"
    echo "  -  The OMM_MACSE pipeline: "
    echo "       Scornavacca C, Belkhir K, Lopez J, Dernat R, Delsuc F, Douzery EJP, Ranwez V. 2019 OrthoMaM v10: Scaling-Up Orthologous Coding Sequence and Exon Alignments with More than One Hundred Mammalian Genomes. Molecular Biology and Evolution 36, 861–862. (doi:10.1093/molbev/msz015)"
    echo "  -  The MACSE v2 program: "
    echo "       Ranwez V, Douzery EJP, Cambon C, Chantret N, Delsuc F. 2018 MACSE v2: Toolkit for the Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons. Molecular Biology and Evolution 35, 2582–2584. (doi:10.1093/molbev/msy159)"
    echo "  -  GNU Parallel: "
    echo "       Tange O. 2011 GNU Parallel - The Command-Line Power Tool. login: The USENIX Magazine, 42-47."
    echo "  -  HmmCleaner: "
    echo "       Di Franco A, Poujol R, Baurain D, Philippe H. 2019 Evaluating the usefulness of alignment filtering methods to reduce the impact of errors on evolutionary inferences. BMC Evol Biol 19, 21. (doi:10.1186/s12862-019-1350-2)"
    echo "  -  AMAS:"
    echo "       Borowiec ML. 2016 AMAS: a fast tool for alignment manipulation and computing of summary statistics. PeerJ 4, e1660. (doi:10.7717/peerj.1660)"
    echo "  -  IQTREE3"
    echo "       Wong T et al. 2025 IQ-TREE 3: Phylogenomic Inference Software using Complex Evolutionary Models. (doi:10.32942/X2P62N)"
    echo "  -  SeqKit:"
    echo "       Shen W, Le S, Li Y, Hu F. 2016 SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS ONE 11, e0163962. (doi:10.1371/journal.pone.0163962)"
    exit $1
}

log() {
    echo "$(date '+%Y/%m/%d %H:%M:%S') - $*"  | tee -a "$LOGFILE"
}

logquiet() {
    echo "$(date '+%Y/%m/%d %H:%M:%S') - $*"  >> "$LOGFILE"
}

logcat() {
    awk -v type="$1" -v date="$(date '+%Y/%m/%d %H:%M:%S')" '{print date " - " type " - " $0 }' $2 | tee -a "$LOGFILE"
}

# Initialize variables
CONF=false
DIR=false
OUT=false
PRE=false
CONTINUE=false
RESTART="NA"
MAXTAX=500
GAP=50
DROP=50
REFS=false
ANN=false
ann_format=false
IQTREE_THRESHOLDS=(NA)
SPTREE_MODEL=false
GENETREES_MODEL=false
NOGENEPART=false
THREADS=$(($(nproc) / 2))
DEBUG=false
CMDLINE=$(echo "$0 $*")
skip_tree=false

# Get options arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input)
	    shift
	    [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after -i|--input"; exit 1; }
	    DIR=$(realpath "$1"); shift ;;
        -o|--output)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after -o|--output"; exit 1; }
            OUT=$(realpath "$1"); shift ;;
        -pre|--prefix)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after -pre|--prefix"; exit 1; }
            PRE=$1; shift ;;
        -n|--ntax)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after -n|--ntax"; exit 1; }
            MAXTAX=$1; shift ;;
        --trim)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after --trim"; exit 1; }
            GAP=$1; shift ;;
        --drop)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after --drop"; exit 1; }
            DROP=$1; shift ;;
        -r|--refs)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after -r|--refs"; exit 1; }
            REFS=$(realpath "$1"); shift ;;
        -a|--annotation)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after -a|--annotation"; exit 1; }
            if [ "$1" = "auto" ];then ann_format="auto";else ANN=$(realpath "$1");fi
            shift ;;
        -p|--perc)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after -p|--perc"; exit 1; }
            IQTREE_THRESHOLDS="$1"; shift ;;
        -s|--sptree)
	    shift
	    if [[ -n "$1" && "$1" != -* ]]; then
		SPTREE_MODEL="$1"
		shift
	    else
                SPTREE_MODEL="MFP+MERGE"
            fi
	    ;;
        -g|--genetrees)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after -g|--genetrees"; exit 1; }
            GENETREES_MODEL="$1"; shift ;;
        --nogenepart)
            NOGENEPART=true
            shift ;;
        --restart)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after --restart"; exit 1; }
            RESTART="$1"; shift ;;
        -y|--continue)
            CONTINUE=true
            shift ;;
        -t|--threads)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after -t|--threads"; exit 1; }
            THREADS=$1; shift ;;
        -c|--config)
            shift
            [[ -z "$1" || "$1" == -* ]] && { echo "[ERROR] - Missing argument after -c|--config"; exit 1; }
            CONF=$(realpath "$1"); shift ;;
        --verbose)
            LOGFILE="/dev/stdout"
            shift ;;
        --debug)
            DEBUG=true
            shift ;;
        -h|--help)
            print_help 0;;
        *)
            echo "[ERROR] - Unknown option: $1"
            print_help 1
            shift ;;
    esac
done

# Check arguments validity
if [[ $DIR == false || $OUT == false || $PRE == false ]]; then
    echo "[ERROR] - Arguments --input, --output and --prefix are required"
    exit 1
fi
if [[ $ANN != false && $REFS == false ]];then
    echo "[ERROR] - Providing an annotation file without a reference file is not currently supported"
    exit 1
fi
if [[ "$PRE" =~ / ]]; then
    echo "[ERROR] - Prefix: \"$PRE\" should not contain any \"/\""
    exit 1
fi
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || (( THREADS < 1 || THREADS > $(nproc) )); then
    echo "[ERROR] - Argument -t|--threads must be an integer between 1 and $nproc"
    exit 1
fi
if ! [[ "$GAP" =~ ^[0-9]+$ ]] || (( GAP < 0 || GAP > 100 )); then
    echo "[ERROR] - Argument --gap must be an integer between 0 and 100"
    exit 1
else
    gap_threshold=$(awk "BEGIN {printf \"%.2f\n\", $GAP/100}")
fi
if ! [[ "$DROP" =~ ^[0-9]+$ ]] || (( DROP < 0 || DROP > 100 )); then
    echo "[ERROR] - Argument --drop must be an integer between 0 and 100"
    exit 1
fi
if ! [[ "$MAXTAX" =~ ^[0-9]+$ ]] || (( MAXTAX < 1 )); then
    echo "[ERROR] - Argument --ntax must be an integer greater or equal to 1"
    exit 1
fi
RESTART_values=("NA" "all" "alignments" "cleaning" "stats" "perc_spp" "sptree" "genetrees")
if [[ ! " ${RESTART_values[*]} " =~ "$RESTART" ]]; then
    echo "[ERROR] - Argument --restart must be one of: ${RESTART_values[*]}"
    exit 1
fi

if [[ "$RESTART" = "all" ]]; then
    echo "$(date '+%Y/%m/%d %H:%M:%S') - [WARNING] - Restarting from scratch, removing output folder"
    sleep 5
    rm -rf "$OUT"
fi

# Initialize output and tmp folder
mkdir -p "$OUT" ; cd "$OUT"

TMP="TMP.phylomera.$PRE"
mkdir -p $TMP

# Create log file
LOGFILE="$PRE.phylomera.log"
if ! compgen -G "$LOGFILE" > /dev/null ; then
    echo "$CMDLINE" > "$LOGFILE"
elif [ "$RESTART" != "NA" ]; then
    log "[INFO] - Appending output to $LOGFILE"
elif $CONTINUE; then
    log "[INFO] - Appending output to $LOGFILE"
else
    while true; do
        read -p "$(echo "$(date '+%Y/%m/%d %H:%M:%S') - [WARNING] - Log file $LOGFILE already exists. Continue writing to the existing log file? (y/n) ")" ans
        case $ans in
            [Yy] )
                log "[INFO] - Appending output to $LOGFILE"
                break
                ;;
            [Nn] )
                echo "$(date '+%Y/%m/%d %H:%M:%S') - [WARNING] - Overwriting $LOGFILE"
                echo "$CMDLINE" > "$LOGFILE"
                break
                ;;
            * )
                echo "$(date '+%Y/%m/%d %H:%M:%S') - [ERROR] - Valid answers: y (yes) or n (no)."
                ;;
        esac
    done
fi
log "[INFO] - Running Phylomera v$VERSION in $OUT with prefix $PRE"

# Get script paths
if [ ! "$CONF" = false ];then 
    if ! compgen -G "$CONF" > /dev/null; then
        log "[ERROR] - $CONF doesn't exist"
        exit 1
    fi
    if ! source $CONF > /dev/null;then
        log "[ERROR] - Error reading config file $CONF"
        exit 1
    fi
    # Replace ~ by $HOME in conf file
    sed -i 's/~/$HOME/g' "$CONF"
    # Check paths
    command -v $MAFFT > /dev/null 2>&1 || { log "[ERROR] - Could not execute mafft from \"$MAFFT\"";exit 1; }
    command -v $HMMCLEANER > /dev/null 2>&1 || { log "[ERROR] - Could not execute HmmCleaner from \"$HMMCLEANER\"";exit 1; }
    command -v $OMM_MACSE > /dev/null 2>&1 || { log "[ERROR] - Could not execute MACSE from \"$OMM_MACSE\"";exit 1; }
    command -v $IQTREE3 > /dev/null 2>&1 || { log "[ERROR] - Could not execute iqtree3 from \"$IQTREE3\"";exit 1; }
    command -v $PARALLEL > /dev/null 2>&1 || { log "[ERROR] - Could not execute gnu-parallel from \"$PARALLEL\"";exit 1; }
    command -v $AMAS > /dev/null 2>&1 || { log "[ERROR] - Could not execute AMAS from \"$AMAS\"";exit 1; }
    command -v $SEQKIT > /dev/null 2>&1 || { log "[ERROR] - Could not execute SeqKit from \"$SEQKIT\"";exit 1; }
    for int_scripts in make_partition_from_probe2.py drop_short_seq.py clean_amas_partition.py;do
        command -v $SCRIPTS/$int_scripts > /dev/null 2>&1 || { log "[ERROR] - Could not execute $SCRIPTS/$int_scripts";exit 1; }
    done
else
    # External scripts must be accessible from PATH and phylomera scripts in the same folder as phylomera
    for ext_scripts in mafft HmmCleaner.pl omm_macse_v11.05b.sif iqtree3 parallel AMAS.py seqkit;do
        if ! which $ext_scripts > /dev/null;then
            log "[ERROR] - Command \"$ext_scripts\" not found in PATH. Consider adding --config <config file> to provide custom paths to external commands"
            exit 1
        else
            VAR=$(echo "${ext_scripts%%.*}"|tr a-z A-Z|sed -E "s/_V[0-9]+//")
            declare "$VAR=$ext_scripts"
        fi
    done
    SCRIPTS=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
    for int_scripts in make_partition_from_probe2.py drop_short_seq.py clean_amas_partition.py;do
        if ! which $SCRIPTS/$int_scripts > /dev/null;then
            log "[ERROR] - Script \"$int_scripts\" not found in $SCRIPTS. Consider adding --config <config file> to provide custom paths to scripts"
            exit 1
        fi
    done
fi
        
# Check input folder
if ! compgen -G "$DIR" > /dev/null; then
    log "[ERROR] - $DIR doesn't exist"
    exit 1
else
    if ! compgen -G "$DIR"/*.fas > /dev/null && ! compgen -G "$DIR"/*.fasta > /dev/null; then
	log "[ERROR] - No FASTA file found in: $DIR (.fas or .fasta expected)"
	exit 1
    else
        # If an output alredy exists, check if inital files are the same (same number of total headers)
        if compgen -G "initial_files.txt" > /dev/null ;then
            out_ntax=$(grep ">" $(cat initial_files.txt)|cut -f2 -d':'|sort -u |wc -l)
            # add -L to follow symlink destination
            cur_ntax=$(find -L "$DIR/" -type f \( -name "*.fas" -o -name "*.fasta" \) -exec grep ">" {} \; |sort -u|wc -l)
            if [ $out_ntax -ne $cur_ntax ];then
                echo "$(date '+%Y/%m/%d %H:%M:%S') - [ERROR] - Initial files that were used by phylomera to produce existing outputs (initial_files.txt) are not the same as input files ($DIR). If the taxon sampling changed, please use a different output name to restart alignments."
                exit 1
            fi
        fi
        
	# Get initial fasta paths
	find -L "$DIR" \( -name "*.fas" -o -name "*.fasta" \) > initial_files.txt
        mkdir -p FASTA
        # check marker names
        sed -E "s:.*/([^\.]+)\..*:\1:g" initial_files.txt > marker_names.txt
        if [ "$(sort -u marker_names.txt|wc -l)" -ne "$(cat marker_names.txt|wc -l)" ];then
            log "[ERROR] - Initial files prefix are not unique. Make sure each FASTA file corresponds to one unique marker ID before first \".\""
            exit 1
        fi
    fi
fi

# Chek ref file
if [ ! "$REFS" = false ]; then
    if ! compgen -G "$REFS" > /dev/null; then
	log "[ERROR] - $REFS doesn't exist"
	exit 1
    fi
    # Check format of ref file
    if ! grep -qE '^>' "$REFS" ; then 
	log "[ERROR] - $REFS is not a valid FASTA"
	exit 1
    fi
    log "[INFO] - Preparing references: $REFS"
    # Check if all markers in input sequences are unique in ref file
    grep ">" "$REFS"|egrep -o "^(\\S+)\\s?"|sed 's/>//g'|sort|uniq -c > $TMP/markers_on_ref.tmp
    grep -Fix -f <(sort marker_names.txt) <(sort $TMP/markers_on_ref.tmp|awk '{ print $2 }') > $TMP/matching_markers_ref.tmp
    if [ "$(cat $TMP/markers_on_ref.tmp|wc -l)" -lt "$(cat marker_names.txt|wc -l)" ] || \
       [ "$(cat $TMP/matching_markers_ref.tmp|wc -l)" -ne  "$(cat marker_names.txt|wc -l)" ];then
        diff -y <(sort marker_names.txt) <(sort $TMP/markers_on_ref.tmp|awk '{ print $2 }') > mismatch_ref.txt
        log "[ERROR] - Some sequences in FASTA files are missing from the reference file. Make sure marker names match exactly between input files and reference headers in $OUT/mismatch_ref.txt"
        exit 1
    fi
    if [ "$(awk '{print $1}' $TMP/markers_on_ref.tmp |sort -u|wc -l)" -eq 1 ];then
        if [ "$(awk '{print $1}' $TMP/markers_on_ref.tmp |sort -u)" -eq 1 ];then
            $SEQKIT split -i "$REFS" --out-dir "$TMP/" --by-id-prefix "ref_" >> "$LOGFILE" 2>&1
            find $TMP -name "ref_*" > $TMP/reference_list.tmp
            xargs sed -E -i "s/>.*/>REF/g" < $TMP/reference_list.tmp
        else
            log "[ERROR] - Sequences from reference file match with multiple marker names. Check $TMP/markers_on_ref.txt to see what went wrong"
            exit 1
        fi
    else
        log "[ERROR] - Some sequences in reference file are not unique to marker names. Check $TMP/markers_on_ref.txt to see what went wrong"
        exit 1
    fi
fi

# Check annotation file
if [ ! "$ANN" = false ]; then
    if ! compgen -G "$ANN" > /dev/null; then
	log "[ERROR] - $ANN doesn't exist"
	exit 1
    fi
    log "[INFO] - Reading annotation file: $ANN"
    # Check if all markers in input sequences are unique in ann file
    grep -wio -f marker_names.txt "$ANN"|sort|uniq -c > $TMP/markers_on_ann.tmp
    if [ "$(cat $TMP/markers_on_ann.tmp|wc -l)" -eq 0 ];then
        log "[ERROR] - No markers were found in annotation file. Make sure $OUT/marker_names.txt match exactly with the annotation file"
        exit 1
    fi
    grep -Fix -f <(sort marker_names.txt) <(sort $TMP/markers_on_ann.tmp|awk '{ print $2 }') > $TMP/matching_markers_ann.tmp
    if [ "$(cat $TMP/markers_on_ann.tmp|wc -l)" -lt "$(cat marker_names.txt|wc -l)" ] || \
       [ "$(cat $TMP/matching_markers_ann.tmp|wc -l)" -ne  "$(cat marker_names.txt|wc -l)" ];then
        diff -y <(sort marker_names.txt) <(sort $TMP/markers_on_ann.tmp|awk '{ print $2 }') > mismatch_annotation.txt
        log "[ERROR] - Some sequences in FASTA files are missing from the annotation file. Make sure marker names match exactly between input files and annotation file in $OUT/mismatch_annotation.txt"
        exit 1
    fi
    if [ "$(awk '{print $1}' $TMP/markers_on_ann.tmp |sort -u|wc -l)" -eq 1 ];then
        if [ "$(awk '{print $1}' $TMP/markers_on_ann.tmp |sort -u)" -eq 1 ];then
	    ann_format=$(awk '
		/^[[:space:]]*$/ || /^[[:space:]]*#/ { next }
		NF < 2 { print "unknown"; exit }
		NF == 2 || $3 == "" || $3 ~ /^[[:space:]]*$/ { print "probe"; exit }
		$2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && NF >= 4 && $4 != "" { print "bed"; exit }
		{ print "unknown"; exit }
		' "$ANN")
	    if [ "$ann_format" = "unknown" ];then
		log "[ERROR] - Unrecognized annotation file format"
		exit 1
	    elif [ "$ann_format" = "probe" ];then
		# If it's a fasta file, it can be an error or annotations are in headers. trim sequences from file and check if there is any mention of cds in headers
		if [ "$(grep -c ">" "$ANN")" -gt 0 ];then
		    log "[WARNING] - Annotation file is a FASTA file"
		    grep ">" "$ANN" > $TMP/ann.tmp
		    ANN=$TMP/ann.tmp
		fi
		if [ "$(egrep -ic "(cds)" "$ANN")" -eq 0 ];then
		    log "[ERROR] - No CDS detected in annotation file"
		    exit 1
		fi
	    fi
            logquiet "[LOG] - Type of annotation file: $ann_format "
        else
            log "[ERROR] - Markers in annotation file match with multiple marker names. Check $TMP/markers_on_ann.tmp to see what went wrong"
            exit 1
        fi
    else
        log "[ERROR] - Some markers in annotation file are not unique to marker names. Check $TMP/markers_on_ann.tmp to see what went wrong"
        exit 1
    fi
fi

#### Check IQTREE models from a list of valid models?

############ Begin script

# current wd: $OUT

log "[INFO] - Number of markers: $(cat marker_names.txt|wc -l)"

# If aligned files exist, ask user.
if ! compgen -G "FASTA/*.mafft.fasta" > /dev/null ; then # 
    # Check if input fasta is already aligned
    head -n 20 initial_files.txt |while read file;do
	cat "$file"|sed -e '/>/s/$/#/g' -e '/>/s/^/#/g'|tr -d '\n'|tr '#' '\n'|sed '/^$/d'|grep -v '^>'|awk '{print length}'|sort -u|wc -l
    done > $TMP/check_ali.tmp
    if [[ "$(sort -u $TMP/check_ali.tmp|wc -l)" != 1 ]] ;then
	# Files are not aligned
        skip_align=false
    elif [[ "$RESTART" = "alignments" ]]; then
        log "[INFO] - Restarting alignments for already aligned files"
        skip_align=false
    elif $CONTINUE; then
       log "[INFO] - Continuing with currently aligned files"
       cp initial_files.txt aligned_files.txt
       skip_align=true
    else
	while true; do
	    read -p "$(log "[WARNING] - Input files already aligned. Skip MAFFT alignments? (y/n) ")" ans
	    case $ans in
		[Yy] )
		    log "[INFO] - Continuing with currently aligned files"
	            cp initial_files.txt aligned_files.txt
		    skip_align=true
		    break
		    ;;
		[Nn] )
		    skip_align=false
		    break
		    ;;
		* )
		    log "[ERROR] - Valid answers: y (yes) or n (no)."
		    ;;
	    esac
	done
    fi
elif [[ "$RESTART" = "alignments" ]]; then
    log "[INFO] - Restarting alignments"
    rm -r FASTA; mkdir -p FASTA
    skip_align=false
elif $CONTINUE; then
    log "[INFO] - Continuing with existing alignment files"
    find FASTA -name "*.mafft.fasta" > aligned_files.txt
    skip_align=true
else
    while true; do
        read -p "$(log "[WARNING] - Alignments found in $OUT/FASTA/. Continue with existing files? (y/n) ")" ans
        case $ans in
            [Yy] )
                log "[INFO] - Continuing with existing alignment files"
                find FASTA -name "*.mafft.fasta" > aligned_files.txt
                skip_align=true
                break
                ;;
            [Nn] )
                log "[INFO] - Overwriting alignment files"
                rm -r FASTA; mkdir -p FASTA
                skip_align=false
                break
                ;;
            * )
                log "[ERROR] - Valid answers: y (yes) or n (no)."
                ;;
        esac
    done
fi

# Align input files
if [ "$skip_align" = false ];then 
    if [ ! "$REFS" = false ]; then
	log "[INFO] - Aligning input FASTAs with MAFFT linsi model"
	$PARALLEL --jobs "$THREADS" "
	    $MAFFT --adjustdirectionaccurately --localpair --maxiterate 1000 --out $TMP/{/.}.mafft.tmp.fasta {}
	" :::: initial_files.txt >> $LOGFILE 2>&1
	# Get marker number, corresponding reference sequence and align reference sequence on it
        log "[INFO] - Aligning reference sequences to MAFFT alignments"
	$PARALLEL --jobs "$THREADS" "
	    $MAFFT --adjustdirectionaccurately --localpair --maxiterate 1000 --addfragments $TMP/ref_{2}.fasta --out FASTA/{1/.}.mafft.fasta $TMP/{1/.}.mafft.tmp.fasta
	" :::: initial_files.txt ::::+ marker_names.txt >> $LOGFILE 2>&1
    else
	log "[INFO] - Aligning input FASTAs with MAFFT linsi model"
	$PARALLEL --jobs "$THREADS" "
	    $MAFFT --adjustdirectionaccurately --localpair --maxiterate 1000 --out FASTA/{/.}.mafft.fasta {}
	" :::: initial_files.txt >> $LOGFILE 2>&1
    fi
    find FASTA -type f -name "*.mafft.fasta" > aligned_files.txt
    # remove _R_ prefix to reverse sequences 
    $PARALLEL --jobs "$THREADS" sed -i "s/_R_//g" {} :::: aligned_files.txt
fi

log "[INFO] - Files aligned: $(cat aligned_files.txt|wc -l)"

# If phylomera files exist, ask user
if ! compgen -G "raw_phylomera_files.txt" > /dev/null || [ "$skip_align" = false ]; then
    skip_phylomera=false
elif [[ "$RESTART" = "cleaning" ]]; then
    log "[INFO] - Restarting alignment cleaning, removing phylomera files"
    find FASTA -type f ! -name "*.mafft.fasta" -delete
    skip_phylomera=false
elif $CONTINUE; then
    log "[INFO] - Continuing with existing phylomera files"
    skip_phylomera=true
else
    while true; do
        read -p "$(log "[WARNING] - Phylomera files already exist. Continue with existing files? (y/n) ")" ans
        case $ans in
            [Yy] )
                log "[INFO] - Continuing with existing phylomera files"
                skip_phylomera=true
                break
                ;;
            [Nn] )
                log "[INFO] - Removing phylomera files"
                find FASTA -type f ! -name "*.mafft.fasta" -delete
                skip_phylomera=false
                break
                ;;
            * )
                log "[ERROR] - Valid answers: y (yes) or n (no)"
                ;;
        esac
    done
fi

## Cleaning alignments
if [ "$skip_phylomera" = false ];then
    # If available, split alignment based on ref
    if [ ! "$REFS" = false ]; then
	log "[INFO] - Spliting alignments into core and flanking regions"
	$PARALLEL -j "$THREADS" "
	    python $SCRIPTS/make_partition_from_probe2.py {} REF 20 $TMP/{/.}.part
	    $AMAS split -i {} --split-by $TMP/{/.}.part -f fasta -d dna --remove-empty
	    mv {= s/\..*// =}_part1-out.fas FASTA/{/.}.f1.fasta 
            mv {= s/\..*// =}_part2-out.fas FASTA/{/.}.core.fasta
            mv {= s/\..*// =}_part3-out.fas FASTA/{/.}.f2.fasta 
	" :::: aligned_files.txt >> $LOGFILE 2>&1
###  || { echo [WARNING] - No left flanking data was found for {/.} ; warn=true ; }

        # Remove empty files
        find FASTA -size 0 -print -delete > /dev/null
        # save a tmp list of files for future removal
	find FASTA -type f \( -name "*.f1.fasta" -o -name "*.core.fasta" -o -name "*.f2.fasta" \) > $TMP/split_files.txt
	# Replace aligned files with core and flanking sequences
        cp $TMP/split_files.txt aligned_files.txt
    fi
    # keeping resulting alignment from previous linsi with ref and not realigning all with einsi seems more conservative and better after HMMCleaner
   
 
    # Make a one liner fasta and remove _R_ flag from reverse sequences and REF sequence
    log "[INFO] - Reformatting alignments"
    $PARALLEL -j "$THREADS" "cat {}  | sed -e '/>/s/$/#/g' -e '/>/s/^/#/g'| tr -d '\n'| tr '#' '\n'| sed '/^$/d' | sed -E '/>/s/>[ _]R[ _]/>/g'|sed '/>REF/{N;d;}' > FASTA/{/.}.ok.fasta" :::: aligned_files.txt
    find FASTA -type f -name "*.ok.fasta" > aligned_files.txt

    # Annotate coding and non coding regions
    if [ "$ann_format" = "probe" ];then
        # Probe annotation only
	log "[INFO] - Processing coding markers with OMM_MACSE based on annotations file"
	mkdir -p $TMP/MACSE
	grep -i "CDS" "$ANN" |grep -wio -f marker_names.txt > $TMP/coding_markers.tmp
        log "[INFO] - Coding markers: $(cat $TMP/coding_markers.tmp|wc -l)"
        grep -vwi -f $TMP/coding_markers.tmp marker_names.txt > $TMP/non_coding_markers.tmp
	grep -wi -f $TMP/coding_markers.tmp aligned_files.txt | grep "core" > $TMP/macse_files.txt
	# If no output is given by macse, flag the sequence so that it is processed as non-coding
	$PARALLEL --plus -j "$THREADS" "
	    $OMM_MACSE --out_dir $TMP/MACSE/{/..} --out_file_prefix {/..} --in_seq_file {} --replace_FS_by_gaps --min_percent_NT_at_ends $gap_threshold
	    if [ -s $TMP/MACSE/{/..}/{/..}_final_align_NT.aln ];then
		mv $TMP/MACSE/{/..}/{/..}_final_align_NT.aln {..}-cds.ok.macse.fasta
	    else
		mv {} {..}-nc.ok.fasta
                logquiet "[WARNING] - MACSE retrieved no results for {} and will be treated as non-coding"
	    fi
	" :::: $TMP/macse_files.txt >> "$LOGFILE" 2>&1
	# Update aligned_files to exclude cds that were already cleaned by macse
        grep -wi -f $TMP/non_coding_markers.tmp aligned_files.txt | grep "core" > $TMP/core_nc_files.txt
        # add core files that retrieved no results after MACSE
        find FASTA -type f -name "*core-nc.ok.fasta" >> $TMP/core_nc_files.txt
        log "[INFO] - Coding markers with no MACSE output: $(cat $TMP/core_nc_files.txt|wc -l)"
        # get flanking files
	find FASTA -type f \( -name "*.f1.ok.fasta" -o -name "*.f2.ok.fasta" \) > $TMP/flanking_files.txt
        cat $TMP/core_nc_files.txt $TMP/flanking_files.txt > aligned_files.txt
        
    elif [ "$ann_format" = "auto" ];then
	# Detect coding alignments with OMM_MACSE: non coding regions shouldn't output any final alignment
	log "[INFO] - Inferring coding and non-coding regions de novo with OMM_MACSE"
	mkdir -p $TMP/MACSE
        ### Beter with splited files, maybe implement a way to split the alignments to lose as few nt as possible
	$PARALLEL --plus -j "$THREADS" "
	    $OMM_MACSE --out_dir $TMP/MACSE/{/..} --out_file_prefix {/..} --in_seq_file {} --replace_FS_by_gaps --min_percent_NT_at_ends $gap_threshold
	    if [ -s $TMP/MACSE/{/..}/{/..}_final_align_NT.aln ];then
                mv $TMP/MACSE/{/..}/{/..}_final_align_NT.aln {..}-cds.ok.macse.fasta
            else
                mv {} {..}-nc.ok.fasta
                logquiet "[WARNING] - MACSE retrieved no results for {} and will be treated as non-coding"
            fi
	" :::: aligned_files.txt >> $LOGFILE 2>&1
        # Update aligned_files to exclude cds files
        find FASTA -type f -name "*-nc.ok.fasta" > aligned_files.txt

    elif [ "$ann_format" = "bed" ];then 
	# bed : trim and elongate, then split according to CDS boundaries
	log "[INFO] - Getting coding regions from bed annotation file (not implemented yet)"
        exit 1
    fi
	
    
    # Process non-coding / unannotated alignments with HMMCleaner
    
    ####  Consider making multiple rounds?
    log "[INFO] - Cleaning alignments with HMMcleaner"
    $PARALLEL --plus -j "$THREADS" "
	echo Processing {.}_hmm.fasta
	sed -i '/>/! s/-/N/g' {}
	$HMMCLEANER -profile=leave-one-out --specificity {}
	sed -i '/>/s/ /_/g' {.}_hmm.fasta
    " :::: aligned_files.txt >> $LOGFILE 2>&1

    log "[INFO] - Removing alignment columns with less than $GAP% characters"
    # Remove empty files that cause unnecesssary errors
    find FASTA -size 0 -print -delete  > /dev/null
    $PARALLEL --plus -j "$THREADS" "
	$AMAS trim -i {.}_hmm.fasta -t $gap_threshold -f fasta -d dna
	mv trimmed_{/.}_hmm.fasta-out.fas {.}_hmm.t$GAP.fasta
    " :::: aligned_files.txt >> $LOGFILE 2>&1

    # Apply following command to all alignments, coding or non-coding
    find FASTA -type f \( -name "*.ok_hmm.t$GAP.fasta" -o -name "*.ok.macse.fasta" \) > cleaned_files.txt

    log "[INFO] - Removing sequences with less than $DROP% characters in alignments"
    $PARALLEL --plus -j "$THREADS" "python $SCRIPTS/drop_short_seq.py {} {}_dropped$DROP $DROP" :::: cleaned_files.txt >> $LOGFILE 2>&1
    
    # Only keep specimen names in headers
    log "[INFO] - Cleaning alignments and temporary files"
    mkdir -p $TMP/FASTA
    clean_ali() {
	i=$1
	init=$2
	if [ ! -f "$i" ];then echo "$i does not exist";return 0 ;fi
	sed -E "/>/s/__.*$//g" "$i" > "${i}.cleannames.fasta" # output from aTRAM pipeline regroup_uce.sh
	sed -E -i "s/[| ]+.*$//g" "${i}.cleannames.fasta" # other output, e.g. from IBA
	sed -i "/>/! s/.*/\U&/g" "${i}.cleannames.fasta" # make all dna uppercase
	mv "${i}.cleannames.fasta" ${init/.mafft/.phylomera} # rename the final alignment files
    }
    export -f clean_ali
    $PARALLEL --plus -j "$THREADS" "clean_ali {}_dropped$DROP {...}.fasta" :::: cleaned_files.txt >> $LOGFILE 2>&1

    # Clean FASTA dir
    find FASTA -name '*.ok*' -print0 | xargs -0 mv -t $TMP/FASTA
    # Remove empty files that make following scripts crash
    find FASTA -size 0 -print -delete  > /dev/null
    # Get list of phylomera files
    find FASTA -name "*.phylomera.*fasta" > raw_phylomera_files.txt
    if [ ! -s raw_phylomera_files.txt ];then
	log "[ERROR] - Something went wrong during alignment cleaning, no phylomera files found"
	exit 1
    fi
    if [[ "$(egrep -o -f marker_names.txt raw_phylomera_files.txt |sort -u|wc -l)" -ne "$(cat marker_names.txt|wc -l)" ]]; then
        log "[WARNING] - Some markers are missing from phylomera files: $(egrep -o -f marker_names.txt raw_phylomera_files.txt |sort -u|wc -l) remaining"
    fi
fi    

# Check stats
if ! compgen -G "taxa.txt" > /dev/null || [ "$skip_phylomera" = false ]; then
    # Stats don't exist yet
    skip_stats=false
elif [[ "$RESTART" = "stats" ]]; then
    log "[INFO] - Restarting alignments statistics"
    skip_stats=false
elif $CONTINUE; then
    log "[INFO] - Skipping alignments statistics"
    skip_stats=true
else
    while true; do
	read -p "$(log "[WARNING] - Alignments statistics already exist. Continue with existing files? (y/n) ")" ans
	case $ans in
	    [Yy] )
                log "[INFO] - Skipping alignments statistics"
                skip_stats=true
		break
		;;
	    [Nn] )
		skip_stats=false
                break
		;;
	    * )
		log "[ERROR] - Valid answers: y (yes) or n (no)"
		;;
	esac
    done
fi

if [ "$skip_stats" = false ]; then
    log "[INFO] - Writing raw alignments summary into $PRE.raw_summary.txt using AMAS"
    # Overall summary with AMAS
    $AMAS summary -d dna -f fasta -c "$THREADS" -i $(cat raw_phylomera_files.txt) -o $PRE.raw_summary.txt >> "$LOGFILE" 2>&1

    # Remove sequences < 20 bp and with less than 4 sequences
    tail -n +2 $PRE.raw_summary.txt| awk '$3>=20 && $2>=4 {print "FASTA/" $1}' > phylomera_files.txt

    if [[ ! -s phylomera_files.txt ]] ;then
        log "[WARNING] - There are less than 4 taxa in this dataset. Statistics will be computed on all alignments."
        tail -n +2 $PRE.raw_summary.txt| awk '{print "FASTA/" $1}' > phylomera_files.txt
        skip_tree=true
    else
        log "[INFO] - Removed alignments shorter than 20 bp and with less than 4 taxa"
    fi

    # Total number of taxa based on phylomera files headers
    if [[ -s phylomera_files.txt ]] ;then
        cat $(cat phylomera_files.txt)|grep ">"|sort -u|sed "s/>//g" > taxa.txt
    else
        log "[ERROR] - No alignments found."
        exit 1
    fi
    ntax=$(cat taxa.txt|wc -l)

    # Check the number of taxa found
    if [[ "$ntax" -gt "$MAXTAX" ]]; then
	log "[ERROR] - More than $MAXTAX individuals detected when combining markers"
	log "[ERROR] - If that's not supposed to be the case, phylomera files are probably wrong. Check taxa.txt and phylomera files headers to see what went wrong"
	log "[ERROR] - If the dataset has more than $MAXTAX individuals, you can specify the exact number of individuals with --ntax <number>"
	exit 1
    fi
    log "[INFO] - Number of taxa remaining after alignment cleaning: $ntax"

    $PARALLEL -j "$THREADS" 'echo -e "{}\t$(cat $(egrep -wi "{}" phylomera_files.txt)|grep ">"|sort -u|wc -l)"' :::: marker_names.txt > taxa_per_marker.txt
    cut -f1-2 $PRE.raw_summary.txt |tail -n +2 > ntaxa_per_file.txt
    grep -f taxa.txt $(cat phylomera_files.txt) > $TMP/marker_per_taxa.tmp
    paste <(grep -wio -f marker_names.txt $TMP/marker_per_taxa.tmp) <(cut -f2 -d'>' $TMP/marker_per_taxa.tmp) |sort -u|cut -f2|sort|uniq -c|sort -k1,1 -nr > marker_per_taxa.txt
    log "[INFO] - Wrote taxa_per_marker.txt and marker_per_taxa.txt"
else
    ntax=$(cat taxa.txt|wc -l)
fi

# Run phylogenetic analyses
if ! compgen -G "perc_spp.txt" > /dev/null || [ "$skip_stats" = false ]; then
    skip_perc=false
elif [[ "$RESTART" = "perc_spp" ]]; then
    log "[INFO] - Restarting supermatrix completeness computation"
    skip_perc=false
elif $CONTINUE; then
    skip_perc=true
else
    while true; do
	read -p "$(log "[WARNING] - perc_spp.txt already exists. Continue with existing file? (y/n) ")" ans
	case $ans in
	    [Yy] )
                skip_perc=true
		break
		;;
	    [Nn] )
		skip_perc=false
                break
		;;
	    * )
		log "[ERROR] - Valid answers: y (yes) or n (no)."
		;;
	esac
    done
fi

if [ "$skip_perc" = "false" ]; then
    log "[INFO] - Writing marker/taxa matrix in perc_spp.txt"
    VALID_THRESHOLDS=$(echo {0..100..5})

    echo -e "Perc\tTaxa\tRemaining_taxa\tMarkers\tDiff_to_next" > perc_spp.txt
    for perc in $VALID_THRESHOLDS;do
	ntaxa_perc=$(awk -v ntax="$ntax" -v perc="$perc" '{ print int(perc/100*ntax) }' <(echo ""))
        # Get all files for that percentage
        awk -v ntaxa_perc="$ntaxa_perc" '$2 >= ntaxa_perc { print "FASTA/" $1 }' ntaxa_per_file.txt > $TMP/list_${perc}spp_files.tmp
        if [ -s $TMP/list_${perc}spp_files.tmp ];then
            # list of file is not empty
	    cat $(cat $TMP/list_${perc}spp_files.tmp)|grep ">"|sort -u|sed "s/>//g" > $TMP/taxa_${perc}spp.tmp
	    ntaxa_remain=$(cat $TMP/taxa_${perc}spp.tmp|wc -l)
	    n_marker=$(grep -wio -f marker_names.txt $TMP/list_${perc}spp_files.tmp |sort -u| wc -l)
	    echo -e "${ntaxa_remain}\t${perc}\t${n_marker}\t${ntaxa_perc}"
        else
            echo -e "0\t${perc}\t0\t${ntaxa_perc}"
        fi
    done > $TMP/perc_spp.tmp
    awk '
    {
	ind[NR] = $1
	perc[NR] = $2
	val[NR] = $3
        tax[NR] = $4
    }
    END {
	for (i = 1; i <= NR; i++) {
	    if (i < NR) {
		diff = (val[i] - val[i+1]) * -1
	    } else {
		diff = "NA"
	    }
	    print perc[i] "\t" tax[i] "\t" ind[i] "\t" val[i] "\t" diff
	}
    }' $TMP/perc_spp.tmp >> perc_spp.txt

    # If no --perc given, let the user choose according to the shape of the curve linking number of markers and number of taxa.
    if $CONTINUE; then
        if [[ "${IQTREE_THRESHOLDS[@]}" == "NA" ]]; then
            IQTREE_THRESHOLDS=0
        fi
    else
	if [[ "${IQTREE_THRESHOLDS[@]}" == "NA" ]]; then
	    # interactive
	    logcat "[INFO]" perc_spp.txt
	    log "[INFO] - Higher percentages minimize biases in phylogenetic reconstruction"
	    log "[INFO] - The highest \"Diff_to_next\" value is usually a good tradeoff between taxon sampling and the number of markers"
	    log "[INFO] - You can also take into account the number of taxa remaining after filtration. Taxa lost during analyses will be written to lost_taxa.txt"
	    while true;do
		read -p "$(log "[INFO] - Please select a percentage of taxa to build the supermatrix on. Multiple values will generate multiple supermatrices: ")" IQTREE_THRESHOLDS
		IQTREE_THRESHOLDS=$(echo "$IQTREE_THRESHOLDS" | tr ',' ' ')
		ALL_VALID=true
		for val in $IQTREE_THRESHOLDS; do
		    if ! [[ " $VALID_THRESHOLDS " =~ (^|[[:space:]])"$val"($|[[:space:]]) ]]; then
			ALL_VALID=false
			break
		    fi
		done

		if $ALL_VALID; then
		    break
		else
		    log "[ERROR] - Unknown value. Allowed values are: $VALID_THRESHOLDS"
		fi
	    done
	fi
    fi
fi

if [[ "$skip_tree" = false ]];then
    # Loop over percentages to compute species trees
    for pct in $IQTREE_THRESHOLDS; do

	log "[INFO] - Creating nexus and phylip files for $pct% species per marker"
	export pct ntax
	$PARALLEL -j "$THREADS" '
	    if [[ $(grep -c ">" {}) -ge $(awk -v ntax="$ntax" -v perc="$pct" "BEGIN{print int(perc/100*ntax)}") ]]; then
		echo {} 
	    fi
	' :::: phylomera_files.txt > $TMP/list_${pct}perc_spp.tmp
	sort $TMP/list_${pct}perc_spp.tmp > list_${pct}perc_spp.txt

	# Merge files from that list into one nexus supermatrix
	$AMAS concat -i $(cat list_${pct}perc_spp.txt) -t $PRE.$pct.phylip -u phylip -p $TMP/$PRE.$pct.nex.tmp -y nexus -d dna -f fasta >> "$LOGFILE" 2>&1
	paste -d'_' <(grep -wio -f marker_names.txt list_${pct}perc_spp.txt) <(sed -E 's:.*phylomera\.(.*)\.fasta:\1:g' list_${pct}perc_spp.txt|tr "-" "_") > $TMP/list_${pct}perc_spp.parts
        log "[INFO] - Number of taxa in $pct% supermatrix: $(head -n 1 $PRE.$pct.phylip|cut -f1 -d' ')"
        log "[INFO] - Size of $pct% supermatrix: $(head -n 1 $PRE.$pct.phylip|cut -f2 -d' ') nt"
        
	# Modify AMAS nexus file
	python $SCRIPTS/clean_amas_partition.py -i $TMP/$PRE.$pct.nex.tmp -o $PRE.$pct.nex -c cds -l $TMP/list_${pct}perc_spp.parts >> "$LOGFILE" 2>&1
        log "[INFO] - Number of initial partitions: $(cat $TMP/list_${pct}perc_spp.parts|wc -l)"

        if [[ "$RESTART" = "sptree" ]]; then
            redo="--redo"
        else
            redo=""
        fi
	# Run IQTREE if specified
	if [ ! "$SPTREE_MODEL" = false ]; then
	    mkdir -p TREES
	    # option -rclusterf 10 should not cause any problem if no MERGE is in species tree
	    log "[INFO] - Running species tree for $pct% spp in TREES/ with command: iqtree3 -s $PRE.$pct.phylip -p $PRE.$pct.nex -m $SPTREE_MODEL -rclusterf 10 -pre TREES/$PRE.$pct.${SPTREE_MODEL/+/} -bb 1000 -alrt 1000 -bnni -nt AUTO -ntmax $THREADS -pers 0.2 -safe $redo"
	    $IQTREE3 -s "$PRE.$pct.phylip" -p "$PRE.$pct.nex" -m $SPTREE_MODEL -pre TREES/$PRE.$pct.${SPTREE_MODEL/+/} -bb 1000 -alrt 1000 -bnni -nt AUTO -ntmax $THREADS -pers 0.2 -safe $redo > /dev/null 2>&1
	else
	    log "[INFO] - Example command to build a species tree: iqtree3 -s $PRE.$pct.phylip -p $PRE.$pct.nex -m MFP+MERGE -rclusterf 10 -pre $PRE.$pct.MFPMERGE -bb 1000 -alrt 1000 -bnni -nt AUTO -ntmax $THREADS -pers 0.2 -safe $redo"
	fi
    done

    ########## Compute genetrees oriented with sptree with Generax?
    # Compute genetrees 
    if [ ! "$GENETREES_MODEL" = false ]; then
	# Compute genetrees on lowest percentage if there are multiple
	pct=$(echo "$IQTREE_THRESHOLDS" | tr ' ' '\n' | sort -n | head -n 1)
	mkdir -p GENETREES
	# Get nexus and phylip files of each marker
	grep -wio -f marker_names.txt list_${pct}perc_spp.txt|sort -u > $TMP/marker_list_$pct.tmp
        
	export AMAS SCRIPTS LOGFILE TMP pct # export variables to let parallel get access to them
	$PARALLEL -j "$THREADS" '
	    grep -wi "{}" list_${pct}perc_spp.txt > $TMP/{}.genetree_files.tmp
	    sed -E "s:FASTA/([^\./]+).*phylomera\.(.*)\.fasta:\1_\2:g" $TMP/{}.genetree_files.tmp|tr "-" "_" > $TMP/{}.genetree_files.parts
	    $AMAS concat -i $(cat $TMP/{}.genetree_files.tmp) -t GENETREES/{}.phylip -u phylip -p $TMP/{}.nex.tmp -y nexus -d dna -f fasta >> "$LOGFILE" 2>&1
	    python $SCRIPTS/clean_amas_partition.py -i $TMP/{}.nex.tmp -o GENETREES/{}.nex -c cds -l $TMP/{}.genetree_files.parts >> "$LOGFILE" 2>&1
	' :::: $TMP/marker_list_$pct.tmp
	find GENETREES -type f -name "*.phylip" > genetrees_files.txt
        log "[INFO] - Number of genetrees to compute: $(cat genetrees_files.txt|wc -l)"

        if [[ "$RESTART" = "genetrees" ]]; then
            redo="--redo"
        else
            redo=""
        fi
	
	if [ "$NOGENEPART" = false ]; then
	    log "[INFO] - Computing genetrees in GENETREES/ with partitions"
	    # make genetrees with partitions (e.g. for ASTRAL)
	    $PARALLEL -j $((THREADS / 2)) $IQTREE3 -s {} -p {.}.nex -m $GENETREES_MODEL -pre GENETREES/{/.}.${GENETREES_MODEL//+/} -bb 1000 -alrt 1000 -bnni -nt 2 -pers 0.2 -safe $redo :::: genetrees_files.txt > /dev/null 2>&1
	else
	    log "[INFO] - Computing genetrees in GENETREES/ without partitions"
	    # no gene partitions
	    $PARALLEL -j $((THREADS / 2)) $IQTREE3 -s {} -m $GENETREES_MODEL -pre GENETREES/{/.}.nogenepart.${GENETREES_MODEL//+/} -bb 1000 -alrt 1000 -bnni -nt 2 -pers 0.2 -safe $redo :::: genetrees_files.txt > /dev/null 2>&1
	fi
    fi
fi

if [ "$DEBUG" = false ];then
    rm -rf $TMP
fi

log "[INFO] - Phylomera finished running."
