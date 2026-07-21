#!/usr/bin/env bash

# wrapper for batch execution
# by M. Soldani, 2026 - developed with OpenAI Codex (GPT-5.5)

# usage:
# ./buildbatch.sh <run_list_file> [--fers-dir DIR] [--digi-path PATH] [--sync-dir DIR] [--output-dir DIR] [--recreate 0|1] [--dry-run]
#
# example:
# ./buildbatch.sh ../data/input/GLOB_runList_total.txt --digi-path /eos/experiment/newtile/beamtests/26_05_t10/digi_root/splitted --sync-dir ../data/output/sync_lists --output-dir /eos/experiment/newtile/beamtests/26_05_t10/global_root/splitted --recreate 1

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"

run_list=""
fers_dir="/eos/experiment/newtile/beamtests/26_05_t10/fers_root/merged/test_align_TStamp_0"
digi_path="/eos/experiment/newtile/beamtests/26_05_t10/digi_root/splitted"
sync_dir="${repo_dir}/data/output/sync_lists"
output_dir="/eos/experiment/newtile/beamtests/26_05_t10/global_root/splitted"
recreate_outputs="1"
dry_run=0
run_count=0
fail_count=0
skip_count=0
file_skip_count=0

usage()
{
    cat <<EOF
Usage:
  $(basename "$0") RUN_LIST [options]

RUN_LIST must have a basename starting with DIGI, FERS, or GLOB.

Options:
  --fers-dir DIR       Directory containing Run<FERS_ID>.dat.root files
                       [${fers_dir}]
  --digi-path PATH     Digitiser ROOT file or directory passed to builder/builder_digi
                       [${digi_path}]
  --sync-dir DIR       Directory containing <FERS_ID>.txt sync lists
                       [${sync_dir}]
  --output-dir DIR     Output directory passed to the builders
                       [${output_dir}]
  --recreate 0|1       Forwarded recreate_outputs flag [${recreate_outputs}]
  --dry-run            Print commands without running them
  -h, --help           Show this help

Dispatch rules:
  DIGI* lists: one digitiser run ID per row -> builder_digi
  FERS* lists: one FERS run ID per row      -> builder_fers
  GLOB* lists: digitiser,FERS per row:
               FERS=-1      -> builder_digi
               digitiser=-1 -> builder_fers
               both present -> builder
EOF
}

die()
{
    echo "Error: $*" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --fers-dir)
            [[ $# -ge 2 ]] || die "--fers-dir requires a value"
            fers_dir="$2"
            shift 2
            ;;
        --digi-path)
            [[ $# -ge 2 ]] || die "--digi-path requires a value"
            digi_path="$2"
            shift 2
            ;;
        --sync-dir)
            [[ $# -ge 2 ]] || die "--sync-dir requires a value"
            sync_dir="$2"
            shift 2
            ;;
        --output-dir)
            [[ $# -ge 2 ]] || die "--output-dir requires a value"
            output_dir="$2"
            shift 2
            ;;
        --recreate)
            [[ $# -ge 2 ]] || die "--recreate requires 0 or 1"
            recreate_outputs="$2"
            shift 2
            ;;
        --dry-run)
            dry_run=1
            shift
            ;;
        -*)
            die "unknown option: $1"
            ;;
        *)
            [[ -z "${run_list}" ]] || die "multiple run lists provided"
            run_list="$1"
            shift
            ;;
    esac
done

[[ -n "${run_list}" ]] || { usage; exit 1; }
[[ -f "${run_list}" ]] || die "run list not found: ${run_list}"
[[ "${recreate_outputs}" == "0" || "${recreate_outputs}" == "1" ]] || die "--recreate must be 0 or 1"

list_name="$(basename "${run_list}")"
case "${list_name}" in
    DIGI*) list_kind="DIGI" ;;
    FERS*) list_kind="FERS" ;;
    GLOB*) list_kind="GLOB" ;;
    *) die "run-list basename must start with DIGI, FERS, or GLOB: ${list_name}" ;;
esac

require_executable()
{
    local exe="$1"
    [[ -x "${exe}" ]] || die "builder executable not found or not executable: ${exe}"
}

require_binary_marker()
{
    local exe="$1"
    local marker="$2"
    require_executable "${exe}"
    grep -aFq "${marker}" "${exe}" || die "$(basename "${exe}") does not look like it was built from ${marker}; recompile it."
}

run_cmd()
{
    printf '+'
    printf ' %q' "$@"
    printf '\n'

    ((run_count += 1))
    if [[ "${dry_run}" -eq 0 ]]; then
        if "$@"; then
            return 0
        fi

        local status=$?
        ((fail_count += 1))
        echo "Warning: command failed with exit status ${status}; continuing with next run." >&2
        return 0
    fi
}

run_cmd_for_file()
{
    local cleanup_path="$1"
    shift

    if [[ "${recreate_outputs}" == "0" && -e "${cleanup_path}" ]]; then
        ((file_skip_count += 1))
        echo "Keeping existing ${cleanup_path}; skipping before opening input ROOT file."
        return 0
    fi

    printf '+'
    printf ' %q' "$@"
    printf '\n'

    ((run_count += 1))
    if [[ "${dry_run}" -eq 0 ]]; then
        if "$@"; then
            return 0
        fi

        local status=$?
        ((fail_count += 1))
        if [[ -e "${cleanup_path}" ]]; then
            rm -f -- "${cleanup_path}" || \
                echo "Warning: failed to remove incomplete output ROOT file ${cleanup_path}." >&2
        fi
        echo "Warning: command failed with exit status ${status}; discarded ${cleanup_path} and continuing with next file." >&2
        return 0
    fi
}

strip_inline_comment()
{
    local line="$1"
    line="${line%%#*}"
    line="${line//,/ }"
    printf '%s\n' "${line}"
}

digi_files=()
expand_digi_inputs()
{
    local digi_id="$1"
    digi_files=()

    if [[ -d "${digi_path}" ]]; then
        local restore_nullglob=0
        shopt -q nullglob || restore_nullglob=1
        shopt -s nullglob
        digi_files=( "${digi_path}/${digi_id}_"*.root )
        if [[ "${restore_nullglob}" -eq 1 ]]; then
            shopt -u nullglob
        fi
        local regular_files=()
        local digi_file
        for digi_file in "${digi_files[@]}"; do
            [[ -f "${digi_file}" ]] && regular_files+=( "${digi_file}" )
        done
        digi_files=( "${regular_files[@]}" )
    else
        digi_files=( "${digi_path}" )
    fi
}

digi_split_id()
{
    local digi_file="$1"
    local digi_id="$2"
    local stem
    stem="$(basename "${digi_file}")"
    stem="${stem%.root}"
    stem="${stem#"${digi_id}_"}"
    printf '%s\n' "${stem}"
}

run_digi_file()
{
    local digi_id="$1"
    local digi_file="$2"
    local split_id
    split_id="$(digi_split_id "${digi_file}" "${digi_id}")"
    run_cmd_for_file "${output_dir}/XXXX_${digi_id}_${split_id}.root" \
        "${script_dir}/builder_digi" "${digi_file}" "${digi_id}" "${output_dir}" "${recreate_outputs}"
}

run_digi()
{
    local digi_id="$1"
    expand_digi_inputs "${digi_id}"

    if [[ "${#digi_files[@]}" -eq 0 ]]; then
        ((run_count += 1))
        ((fail_count += 1))
        echo "Warning: no digi ROOT files found for run ${digi_id} in ${digi_path}; continuing with next run." >&2
        return 0
    fi

    local digi_file
    for digi_file in "${digi_files[@]}"; do
        run_digi_file "${digi_id}" "${digi_file}"
    done
}

run_fers()
{
    local fers_id="$1"
    run_cmd "${script_dir}/builder_fers" "${fers_dir}/Run${fers_id}.dat.root" "${output_dir}" "${recreate_outputs}"
}

run_global()
{
    local digi_id="$1"
    local fers_id="$2"
    local sync_list="${sync_dir}/${fers_id}.txt"

    if [[ ! -f "${sync_list}" ]]; then
        ((skip_count += 1))
        echo "Warning: missing sync list ${sync_list}; skipping paired run DIGI=${digi_id}, FERS=${fers_id}." >&2
        return 0
    fi

    expand_digi_inputs "${digi_id}"

    if [[ "${#digi_files[@]}" -eq 0 ]]; then
        ((run_count += 1))
        ((fail_count += 1))
        echo "Warning: no digi ROOT files found for paired run DIGI=${digi_id}, FERS=${fers_id} in ${digi_path}; continuing with next run." >&2
        return 0
    fi

    local digi_file
    local split_id
    for digi_file in "${digi_files[@]}"; do
        split_id="$(digi_split_id "${digi_file}" "${digi_id}")"
        run_cmd_for_file "${output_dir}/${fers_id}_${digi_id}_${split_id}.root" \
            "${script_dir}/builder" "${sync_list}" "${fers_dir}/Run${fers_id}.dat.root" \
            "${digi_file}" "${digi_id}" "${output_dir}" "${recreate_outputs}"
    done
}

case "${list_kind}" in
    DIGI)
        require_binary_marker "${script_dir}/builder_digi" "builder_digi.cpp"
        ;;
    FERS)
        require_binary_marker "${script_dir}/builder_fers" "builder_fers.cpp"
        ;;
    GLOB)
        require_binary_marker "${script_dir}/builder" "builder.cpp"
        require_binary_marker "${script_dir}/builder_digi" "builder_digi.cpp"
        require_binary_marker "${script_dir}/builder_fers" "builder_fers.cpp"
        ;;
esac

line_no=0
while IFS= read -r raw_line || [[ -n "${raw_line}" ]]; do
    ((line_no += 1))
    line="$(strip_inline_comment "${raw_line}")"
    read -r -a fields <<< "${line}"
    [[ "${#fields[@]}" -eq 0 ]] && continue

    case "${list_kind}" in
        DIGI)
            [[ "${#fields[@]}" -ge 1 ]] || die "line ${line_no}: expected digitiser run ID"
            [[ "${fields[0]}" != "-1" ]] || die "line ${line_no}: DIGI run ID is -1"
            run_digi "${fields[0]}"
            ;;
        FERS)
            [[ "${#fields[@]}" -ge 1 ]] || die "line ${line_no}: expected FERS run ID"
            [[ "${fields[0]}" != "-1" ]] || die "line ${line_no}: FERS run ID is -1"
            run_fers "${fields[0]}"
            ;;
        GLOB)
            [[ "${#fields[@]}" -ge 2 ]] || die "line ${line_no}: expected digitiser,FERS run IDs"
            digi_id="${fields[0]}"
            fers_id="${fields[1]}"

            if [[ "${digi_id}" == "-1" && "${fers_id}" == "-1" ]]; then
                die "line ${line_no}: both digitiser and FERS IDs are -1"
            elif [[ "${fers_id}" == "-1" ]]; then
                run_digi "${digi_id}"
            elif [[ "${digi_id}" == "-1" ]]; then
                run_fers "${fers_id}"
            else
                run_global "${digi_id}" "${fers_id}"
            fi
            ;;
    esac
done < "${run_list}"

if [[ "${fail_count}" -gt 0 ]]; then
    echo "Completed ${run_count} command(s) with ${fail_count} failure(s), ${skip_count} skipped paired run(s), and ${file_skip_count} skipped existing output file(s)." >&2
    exit 1
fi

echo "Completed ${run_count} command(s) with ${skip_count} skipped paired run(s) and ${file_skip_count} skipped existing output file(s)."
