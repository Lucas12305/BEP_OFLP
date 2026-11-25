#!/usr/bin/env bash
# Runs instances in every combination of mode, precutting, seperator and update, and puts the results in one csv file.
# The script can be stopped at any moment and started again with the same command, it then continues where it stopped,
# because a run only gets into the csv file when it is finished. A run that was busy when you stopped starts again.
#
# Usage: ./run_all.sh [-r REPEATS] [-j JOBS] [-t TIME_LIMIT] [-l] [-o OUT_DIR] [instance ...]
#   -r  how many times every combination is run (default 1)
#   -j  how many runs at the same time (default 1), keep this lower than the number of cores
#   -t  time limit of SCIP in seconds per run (default 3600)
#   -l  also run the LP method (default only GS, LP is slow)
#   -o  the folder for the results (default results/all), use the same folder to continue
#   without instances it runs all of the instances in initialized_instances, the smallest first
#
# To pause without losing the run that is busy: touch the file PAUSE in the OUT_DIR, the runs that are busy
# finish and no new ones start. Remove the file and start the script again to continue.
#
# Example: ./run_all.sh -r 3 -j 2 -t 600 ofl_3_100i ofl_4_100i

# bash run_all.sh -j 2 -t 600

REPEATS=1; JOBS=1; TIME_LIMIT=3600; LP=0; OUT_DIR="results/all"
while getopts "r:j:t:lo:" opt; do
    case "$opt" in
        r) REPEATS="$OPTARG" ;;
        j) JOBS="$OPTARG" ;;
        t) TIME_LIMIT="$OPTARG" ;;
        l) LP=1 ;;
        o) OUT_DIR="$OPTARG" ;;
        *) echo "Unknown option" >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

cd "$(dirname "$0")" || exit 1
PYTHON="${PYTHON:-$(command -v python3 || command -v python)}"
mkdir -p "$OUT_DIR"
CSV="${OUT_DIR}/summary.csv"
[ -f "$CSV" ] || echo "instance,mode,precutting,seperator,update,repeat,time,L" > "$CSV"

if [ $# -eq 0 ]; then
    INSTANCES=($(ls initialized_instances/*.gms | xargs -n 1 basename | sed 's/\.gms$//' | sort -V))
else
    INSTANCES=("$@")
fi

# One run, the line for the csv is only added when the run is finished
one_run() {
    INST=$1; MODE=$2; PRE=$3; SEP=$4; UPD=$5; REP=$6
    [ -f "${OUT_DIR}/PAUSE" ] && return 0
    KEY="${INST},${MODE},${PRE},${SEP},${UPD},${REP},"
    grep -q -F "$KEY" "$CSV" && return 0
    FILE="${OUT_DIR}/${INST}_mode${MODE}_pre${PRE}_sep${SEP}_upd${UPD}_rep${REP}.txt"
    echo "Running ${KEY%,}"
    "$PYTHON" -u run_OFLP.py "$INST" --mode "$MODE" --precutting "$PRE" --seperator "$SEP" --update "$UPD" --time-limit "$TIME_LIMIT" > "$FILE" 2>&1
    TIME=$(grep "total time" "$FILE" | tail -1 | sed 's/.*total time: \([0-9.]*\) s.*/\1/')
    [ -z "$TIME" ] && TIME=$(grep "Total time" "$FILE" | tail -1 | sed 's/.*Total time: \([0-9.]*\) s.*/\1/')
    L=$(grep "Objective L" "$FILE" | tail -1 | sed 's/.*: //')
    if [ -n "$TIME" ]; then
        echo "${KEY}${TIME},${L}" >> "$CSV"
    else
        echo "Run ${KEY%,} did not finish, see ${FILE}" >&2
    fi
}
export -f one_run
export OUT_DIR CSV PYTHON TIME_LIMIT

# All combinations, update 2 (from the last hull) only exists for GS
COMBOS=()
for INST in "${INSTANCES[@]}"; do
    for REP in $(seq 1 "$REPEATS"); do
        for PRE in 0 1; do
            for SEP in 0 1; do
                for UPD in 0 1 2; do
                    COMBOS+=("$INST 1 $PRE $SEP $UPD $REP")
                done
                if [ "$LP" -eq 1 ]; then
                    for UPD in 0 1; do
                        COMBOS+=("$INST 0 $PRE $SEP $UPD $REP")
                    done
                fi
            done
        done
    done
done

printf '%s\n' "${COMBOS[@]}" | xargs -P "$JOBS" -L 1 bash -c 'one_run $0 $1 $2 $3 $4 $5'
if [ -f "${OUT_DIR}/PAUSE" ]; then
    echo "Paused, remove ${OUT_DIR}/PAUSE and start the script again to continue."
else
    echo "Done, the results are in ${CSV}"
fi
