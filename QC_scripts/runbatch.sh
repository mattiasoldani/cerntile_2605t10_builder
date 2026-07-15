#!/bin/bash

# wrapper to execute QC analysis on a list of FERS runs
# by M. Soldani, 2026 - developed with OpenAI Codex (GPT-5.5)

# usage: ./runbatch.sh <fers_run_list_file> [overwrite=0]
# --> examples: ./runbatch.sh ../data/input/FERS_runList_total.txt 0

set -Eeu
trap 'rc=$?; echo "ERROR: line $LINENO failed with exit code $rc: $BASH_COMMAND" >&2' ERR

usage() {
    printf "Usage: bash %s <FERS_RUN_LIST_FILE> [OVERWRITE]\n" "$0"
    printf "\n"
    printf "  FERS_RUN_LIST_FILE  Path to a text file containing FERS run numbers.\n"
    printf "                      Blank lines and lines starting with # are ignored.\n"
    printf "  OVERWRITE           Optional: 0 or 1. Default: 0.\n"
    printf "                      When 1, runs QC.sh as: OVERWRITE=1 bash QC.sh <RUN>\n"
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
    usage >&2
    exit 1
fi

RUN_LIST=$1
OVERWRITE_FLAG=${2:-0}

if [[ ! -f "${RUN_LIST}" ]]; then
    printf "Run-list file not found: %s\n" "${RUN_LIST}" >&2
    usage >&2
    exit 1
fi

if ! [[ "${OVERWRITE_FLAG}" =~ ^[01]$ ]]; then
    printf "Invalid OVERWRITE value: %s\n" "${OVERWRITE_FLAG}" >&2
    usage >&2
    exit 1
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
QC_SCRIPT="${SCRIPT_DIR}/QC.sh"

if [[ ! -f "${QC_SCRIPT}" ]]; then
    printf "QC.sh not found next to this script: %s\n" "${QC_SCRIPT}" >&2
    exit 1
fi

N_RUNS=0
N_FAILED=0

while IFS= read -r LINE || [[ -n "${LINE}" ]]; do
    LINE=${LINE%%#*}
    read -r RUN EXTRA <<< "${LINE}"

    if [[ -z "${RUN}" ]]; then
        continue
    fi

    if [[ -n "${EXTRA:-}" ]] || ! [[ "${RUN}" =~ ^[0-9]+$ ]]; then
        printf "Skipping invalid run entry: %s\n" "${LINE}" >&2
        continue
    fi

    N_RUNS=$((N_RUNS + 1))
    printf "\n=== Running QC for FERS run %s (OVERWRITE=%s) ===\n" "${RUN}" "${OVERWRITE_FLAG}"

    if ! OVERWRITE="${OVERWRITE_FLAG}" bash "${QC_SCRIPT}" "${RUN}"; then
        printf "QC failed for FERS run %s\n" "${RUN}" >&2
        N_FAILED=$((N_FAILED + 1))
    fi
done < "${RUN_LIST}"

printf "\nProcessed %d FERS run(s). Failed: %d.\n" "${N_RUNS}" "${N_FAILED}"

if [[ "${N_FAILED}" -ne 0 ]]; then
    exit 1
fi
