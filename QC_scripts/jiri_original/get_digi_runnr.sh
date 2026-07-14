#!/bin/bash
if ! [[ "$1" =~ ^[0-9]+$ ]]; then
  printf "Not a valid parameter.\nUsage $0 <RUN_NUMBER>\n"
  exit 1
fi

RUN=$1
DIFF_THRESHOLD=60
DIGIPATH=/eos/experiment/newtile/beamtests/26_05_t10/digi_raw
REFSPATH=/eos/experiment/newtile/beamtests/26_05_t10/fers_daq/bin/DataFiles
DT5202CONVERTER=/eos/home-j/jkvas/DRD6/newtile/tools/dt5202_raw_utils/raw2txt/dt5202txt

#DT5202CONVERTER=/home/kvas/git/dt5202_raw_utils/raw2txt/dt5202txt

#test if the FERS files exist
FERSFILE=""
if [ -f /eos/experiment/newtile/beamtests/26_05_t10/fers_daq/bin/DataFiles/Run${RUN}.0_list.dat ] ; then
    FERSFILE=/eos/experiment/newtile/beamtests/26_05_t10/fers_daq/bin/DataFiles/Run${RUN}.0_list.dat
elif [ -f /eos/experiment/newtile/beamtests/26_05_t10/fers_daq/bin/DataFiles/Run${RUN}_list.dat ] ; then
    FERSFILE=/eos/experiment/newtile/beamtests/26_05_t10/fers_daq/bin/DataFiles/Run${RUN}_list.dat
else
    echo "FERS DATA file for run ${RUN} could not be found in ${REFSPATH} folder"
    exit 1
fi

FERSTS=$(${DT5202CONVERTER} -w ${FERSFILE} | sed -n 's/^#start TS in s = //p')
echo "#FERS start TS = ${FERSTS}" >&2

mapfile -t DIGIDIRS < <(
    find "${DIGIPATH}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' \
    | grep -E '^[0-9]+$' \
    | sort -n
)

# echo "Number of digi directories: ${#DIGIDIRS[@]}"
# echo "First directory timestamp: ${DIGIDIRS[0]}"
# echo "Last directory timestamp: ${DIGIDIRS[-1]}"

NEAREST=""
NEAREST_DIFF=""
for d in "${DIGIDIRS[@]}"; do
    diff=$(( d > FERSTS ? d - FERSTS : FERSTS - d ))
    if [[ -z "$NEAREST_DIFF" || "$diff" -lt "$NEAREST_DIFF" ]]; then
        NEAREST_DIFF="$diff"
        NEAREST="$d"
    fi
done
if [[ -z "$NEAREST" ]]; then
    echo "ERROR: no digi timestamp directories found in ${DIGIPATH}" >&2
    exit 2
fi
if (( NEAREST_DIFF > DIFF_THRESHOLD )); then
    echo "ERROR: nearest digi dir differs by ${NEAREST_DIFF} s, threshold is ${DIFF_THRESHOLD} s" >&2
    echo "ERROR: FERS TS=${FERSTS}, nearest digi TS=${NEAREST}" >&2
    exit 3
fi
echo "Nearest digitizer timestamp: $NEAREST" >&2
echo "Difference: $NEAREST_DIFF seconds">&2
echo "$NEAREST"
exit 0
