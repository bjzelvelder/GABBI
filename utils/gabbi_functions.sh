#!/usr/bin/env bash

# Various functions used inside the gabbi pipeline

debug() {
    if [[ "${GABBI_DEBUG:-0}" -eq 1 ]]; then
        echo "[GABBI][DEBUG] $*" >&2
    fi
	return 0
}
export -f debug

phylomera() { /opt/gabbi/scripts/phylomera*.sh "$@"; }
export -f phylomera

convert_state() {
    local file="$1"
    local anc="${file%.state}"
    local shr
    shr=$(basename "${anc%.anc}")

    # Convert IQ-TREE state files to FASTA
    grep -v '#' "$file" \
	| awk 'NR > 1 { seq[$1] = seq[$1] $3 }
	       END { for (n in seq) print ">" n "\n" seq[n] }' \
	> "${anc}.fasta"
    sed -E "s/>/>${shr}|/g" "${anc}.fasta" \
	| sed -e '/>/s/$/#/g' -e '/>/s/^/#/g' \
	| tr -d '\n' | tr '#' '\n' | sed '/^$/d' \
	> "${anc}.ok.fasta"
}
export -f convert_state
