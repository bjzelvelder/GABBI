#!/usr/bin/env bash

# =========================================================================
# Checkpoint utility
# =========================================================================
# A lightweight checkpoint system. Each completed step writes a sentinel
# file to GABBI_WORKDIR/.gabbi_checkpoints/. On re-run, completed steps are 
# skipped and failed steps are re-run. If a RESTART value if given, steps that
# follow the restart step are re-run

CHECKPOINT_DIR="${GABBI_WORKDIR}/.gabbi_checkpoints"
mkdir -p "${CHECKPOINT_DIR}"

checkpoint_done() {
    local step="$1"
    [[ -f "${CHECKPOINT_DIR}/${step}.done" && ! -f "${CHECKPOINT_DIR}/${step}.fail" ]]
}

checkpoint_mark() {
    local step="$1"
    [[ -f "${CHECKPOINT_DIR}/${step}.fail" ]] && rm "${CHECKPOINT_DIR}/${step}.fail"
    [[ -f "${CHECKPOINT_DIR}/${step}.restart" ]] && rm "${CHECKPOINT_DIR}/${step}.restart"
    touch "${CHECKPOINT_DIR}/${step}.done"
    echo "[GABBI] Step '${step}' completed and checkpointed at $(date -Iseconds)."
}

checkpoint_fail() {
    local step="$1"
    touch "${CHECKPOINT_DIR}/${step}.fail"
    echo "[GABBI] ERROR: Step '${step}' failed (exit code $?). Exiting." >&2
    exit 1
}

checkpoint_restart() {
    [[ -z "$RESTART" ]] && return
    local triggered=0

    # Check that the target step exists in either .done or .fail
    if [[ -z "$(ls "${CHECKPOINT_DIR}"/*${RESTART}*.done 2>/dev/null)" && \
          -z "$(ls "${CHECKPOINT_DIR}"/*${RESTART}*.fail 2>/dev/null)" ]]; then
        echo "[GABBI] WARNING: --restart step ${RESTART} not found in existing checkpoints!" >&2
        return
    fi

    # Reset all .done and .fail from the target step onwards
    for f in $(ls "${CHECKPOINT_DIR}/"*.{done,fail} 2>/dev/null | sort); do
        local step=$(basename "$f")
        step="${step%.done}"; step="${step%.fail}"
        [[ "$step" =~ "$RESTART" ]] && triggered=1
        if [[ $triggered -eq 1 ]]; then
            rm -f "$f"
            touch "${CHECKPOINT_DIR}/${step}.restart"
            verbose "Checkpoint reset: $step"
        fi
    done
}

checkpoint_fail_exists() {
    local step="$1"
    if [[ -f "${CHECKPOINT_DIR}/${step}.fail" ]]; then
        echo "[GABBI] WARNING: Step ${step/step/} previously failed. Removing erroneous files."
        return 0
    fi
    if [[ -f "${CHECKPOINT_DIR}/${step}.restart" ]]; then
        echo "[GABBI] Restarting step ${step/step/}."
        return 0
    fi
    return 1
}
