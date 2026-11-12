#!/usr/bin/env bash
# Runs one or more instances of the obnoxious facility location problem
# and saves everything that gets printed in a text file in the results folder.
#
# Usage: ./run.sh [--mode 0|1] [--plot] [--time-limit SECONDS] instance [instance ...]
#
# Examples:
#   ./run.sh ofl_2_100i
#   ./run.sh --mode 1 --plot --time-limit 600 ofl_2_100i ofl_3_100i
#   ./run.sh initialized_instances/ofl_2_100*.gms

# how to run? 'bash run.sh ofl_2_100i'


# The folder in the repository where the text files (and plots) get saved
OUT_DIR="results"

usage() {
    sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'
}

MODE=1
# 0 is LP, 1 is GS
PLOT=0
# 0 is off, 1 is on
TIME_LIMIT=""
# maximum amount of time the program can run for
INSTANCES=()


while [ $# -gt 0 ]; do
    case "$1" in
        --mode)       MODE="$2"; shift 2 ;;
        --plot)       PLOT=1; shift ;;
        --time-limit) TIME_LIMIT="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        -*)           echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
        *)            INSTANCES+=("$1"); shift ;;
    esac
done

if [ ${#INSTANCES[@]} -eq 0 ]; then
    echo "Give at least one instance." >&2
    usage >&2
    exit 1
fi

# We always work from the folder of this script, so it can be called from anywhere
cd "$(dirname "$0")" || exit 1
mkdir -p "$OUT_DIR"

PYTHON="${PYTHON:-$(command -v python3 || command -v python)}"
if [ -z "$PYTHON" ]; then
    echo "Could not find python, set it with PYTHON=/path/to/python ./run.sh ..." >&2
    exit 1
fi

FAILED=0
for INSTANCE in "${INSTANCES[@]}"; do
    NAME="$(basename "$INSTANCE" .gms)"
    BASE="${OUT_DIR}/${NAME}_mode${MODE}_$(date +%Y%m%d-%H%M%S)"

    ARGS=(--mode "$MODE")
    if [ "$PLOT" -eq 1 ]; then
        ARGS+=(--plot-file "${BASE}.png")
    fi
    if [ -n "$TIME_LIMIT" ]; then
        ARGS+=(--time-limit "$TIME_LIMIT")
    fi

    echo "Running ${NAME} (mode ${MODE}), the output is saved in ${BASE}.txt"
    "$PYTHON" -u run_OFLP.py "$NAME" "${ARGS[@]}" 2>&1 | tee "${BASE}.txt"
    if [ "${PIPESTATUS[0]}" -ne 0 ]; then
        echo "Running ${NAME} failed, see ${BASE}.txt" >&2
        FAILED=1
    fi
done

exit "$FAILED"
