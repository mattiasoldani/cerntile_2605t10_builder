#!/bin/bash
if ! [[ "$1" =~ ^[0-9]+$ ]]; then
    printf "Not a valid parameter.\nUsage bash $0 <RUN_NUMBER> (optional <DIGI_TS>)\n"
    exit 1
fi

RUN=$1
UNIXTIME_MAXDIFF=100

#file paths without trailing '/'
DIGIPATH=/eos/experiment/newtile/beamtests/26_05_t10/digi_raw
FERSPATH=/eos/experiment/newtile/beamtests/26_05_t10/fers_daq/bin/DataFiles
DT5202CONVERTER=/eos/experiment/newtile/beamtests/26_05_t10/QC/scripts/dt5202_raw_utils/raw2txt/dt5202txt #/eos/home-j/jkvas/DRD6/newtile/tools/dt5202_raw_utils/raw2txt/dt5202txt
QCPATH=/eos/experiment/newtile/beamtests/26_05_t10/QC/run${RUN}
TOOLPATH=/eos/experiment/newtile/beamtests/26_05_t10/QC/scripts #/eos/home-j/jkvas/DRD6/newtile/tools

#these parameters can be overriden by a file analysis_overrides.sh in the QC run subdirectory
DIGI_BASE_FROM=1
DIGI_BASE_TO=0 #when DIGI_BASE_TO<DIGI_VASE_FROM, the baseline will not be calculated at all
DIGI_PEAK_FROM=650
DIGI_PEAK_TO=800
DIGI_CHANNELS_TO_ANALYZE=(6 7)

FERS_BOARDS=(0 1) #for histogram
FERS_EXPRESSIONS_TO_CONVERT=( # for event conversion "<board_id>:<expression>
    "0:hg0"
    "1:hg0"
)
FERS_TOA_CHANNELS=( #BOARD:channel
    "0:0"
    "1:0"
)
CORRELATIONS=( # digi_channel:FERS_board:FERS_expression (from converted expression)
    "6:0:hg0"
    "7:1:hg0"
)
CORRELATION_ACCEPT_FACTOR=2 #fraction of events threshold for accepting correlated run. 10 = at least 1/10 events have to be correlated. Beam usage: should be 2
CORRELATION_INITIAL_OFFSET=0 #will be overwritten by automatic guess unless specified in the analysis_overrides.sh
CORRELATION_INITIAL_RANGE=500 #will be overwritten by automatic guess unless specified in the analysis_overrides.sh

set -Eeuo pipefail

echo "Creating ${QCPATH}"
mkdir -p "${QCPATH}"

mapfile -t files < <(printf '%s\n' "${FERSPATH}"/Run"${RUN}".*_list.dat | sort -V)
echo "Run ${RUN} has data in these files: ${files[@]}" >&2

###################################################################################
# check if FERS raw file exist
###################################################################################
FERSFILE=""
if [[ -f "${FERSPATH}/Run${RUN}.0_list.dat" ]] ; then 
    FERSFILE=${FERSPATH}/Run${RUN}.0_list.dat
elif [[ -f "${FERSPATH}/Run${RUN}_list.dat" ]] ; then
    FERSFILE=${FERSPATH}/Run${RUN}_list.dat
else
    echo "FERS DATA file for run ${RUN} could not be found in ${FERSPATH} folder"
    exit 1
fi

###################################################################################
# process FERS histograms
###################################################################################
for BOARD in "${FERS_BOARDS[@]}" ; do
    HISTFILE="${QCPATH}/histogram_run${RUN}_board${BOARD}.txt"
    if [[ ! -f "${HISTFILE}" ]] ; then
	echo "Creating histogram for run ${RUN}, board ${BOARD}" >&2
	${DT5202CONVERTER} -w <(cat  "${files[@]}") \
			   --print_hg_histogram \
			   --print_lg_histogram \
			   --max_bin=8191 \
			   --board_id="${BOARD}" \
			   >"${HISTFILE}"
    else
	echo "histogram for run ${RUN}, board ${BOARD} already exists. skipping" >&2
    fi
    #High gain histogram
    HISTPLOT="${QCPATH}/histogram_run${RUN}_hg_board${BOARD}.png"
    if [[ ! -f "${HISTPLOT}" ]] ; then
	echo "Creating HG histogram plot for run ${RUN}, board ${BOARD}" >&2	
	gnuplot <<EOF
set terminal pngcairo enhanced size 1920,1080 font ",6"
set output "${HISTPLOT}"
set multiplot layout 8,8 title "run ${RUN}, HG, board ${BOARD}"
set xrange [0:8200]; set yrange [0.8:*]
set logscale y; set grid
do for [i=0:63] {
   plot "${HISTFILE}" index i with histeps title "Ch".i
}
unset multiplot
EOF
    else
	echo "histogram plot for run ${RUN}, board ${BOARD} already exists. skipping" >&2
    fi

    #Low gain histogram
    HISTPLOT="${QCPATH}/histogram_run${RUN}_lg_board${BOARD}.png"
    if [[ ! -f "${HISTPLOT}" ]] ; then
	echo "Creating LG histogram plot for run ${RUN}, board ${BOARD}" >&2	
	gnuplot <<EOF
set terminal pngcairo enhanced size 1920,1080 font ",6"
set output "${HISTPLOT}"
set multiplot layout 8,8 title "run ${RUN}, LG,board ${BOARD}"
set xrange [0:1000]; set yrange [0.8:*]
set logscale y; set grid
do for [i=0:63] {
   plot "${HISTFILE}" index int(i+64) with histeps title "Ch".i
}
unset multiplot
EOF
    else
	echo "histogram plot for run ${RUN}, board ${BOARD} already exists. skipping" >&2
    fi
done

###################################################################################
# process FERS events
###################################################################################
for PAIR in "${FERS_EXPRESSIONS_TO_CONVERT[@]}"; do
    IFS=: read -r BOARD EXPR <<< "$PAIR"
    EVENTFILE="${QCPATH}/fers_run${RUN}_events_board${BOARD}_${EXPR}.txt"
    if [[ ! -f "${EVENTFILE}" ]]; then
        echo "generating event list run ${RUN}, board ${BOARD}, expression ${EXPR}" >&2
        "$DT5202CONVERTER" -w <(cat "${files[@]}") \
			   --channel=0 \
			   --board_id="$BOARD" \
			   --custom_expr="$EXPR" \
			   --print_txt_events \
            | awk '/#/' >"${EVENTFILE}"
    else
        echo "FERS events for run ${RUN}, board ${BOARD}, expression ${EXPR} already exist. Skipping" >&2
    fi
done

###################################################################################
# process FERS TOA+TOT in time
###################################################################################
for PAIR in "${FERS_TOA_CHANNELS[@]}"; do
    IFS=: read -r BOARD CHANNEL <<< "$PAIR"
    TOAPLOTFILE="${QCPATH}/TOAEvents_run${RUN}_board${BOARD}_CH${CHANNEL}.png"
    if [[ ! -f "${TOAPLOTFILE}" ]] ; then
	echo "generating TOA overview file: ${TOAPLOTFILE}"
	printf -v GNUPLOT_CMD 'cat %s | %q -w /dev/stdin --board_id %q --channel 0 --print_txt_events --custom_expr %q | grep %q' "${files[*]}" "$DT5202CONVERTER" "$BOARD" "toa$CHANNEL" "#event"
	gnuplot <<EOF
set terminal pngcairo enhanced size 1600,480
set output "${TOAPLOTFILE}"
set yrange [-100:*]
set grid
set title "Run ${RUN}, TOA in time for board ${BOARD}, ch ${CHANNEL}"
set xlabel "event number"
set ylabel "TOA values"
plot '< ${GNUPLOT_CMD}' u 5:7 w p notitle lw 0.3
EOF
    else
	echo "file ${TOAPLOTFILE} exists. skipping"
    fi    
done

###################################################################################
# FERS boards correlation
###################################################################################
# TODO

###################################################################################
# Get the starting timestamp for the FERS and find the nearest digitizer directory
###################################################################################
FERSTS=$(${DT5202CONVERTER} -w ${FERSFILE} --to_time 0.0 | sed -n 's/^#start TS in s = //p')
if ! [[ "$FERSTS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: invalid FERS timestamp: '$FERSTS'" >&2
    # exit 2
else
    echo "#FERS. start TS = ${FERSTS}" >&2
    mapfile -t DIGIDIRS < <(find "${DIGIPATH}" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | grep -E '^[0-9]+$' | sort -n)
    DIGI_TS=""
    DIGI_TS_DIFF=""
    for d in "${DIGIDIRS[@]}"; do
	diff=$(( d > FERSTS ? d - FERSTS : FERSTS - d ))
	if [[ -z "$DIGI_TS_DIFF" || "$diff" -lt "$DIGI_TS_DIFF" ]]; then
            DIGI_TS_DIFF="$diff"
            DIGI_TS="$d"
	fi
    done
    if [[ -z "$DIGI_TS" ]]; then
	echo "ERROR: no digi timestamp directories found in ${DIGIPATH}" >&2
	exit 2
    fi
    if (( DIGI_TS_DIFF > UNIXTIME_MAXDIFF )); then
	echo "ERROR: nearest digi dir differs by ${DIGI_TS_DIFF} s, threshold is ${UNIXTIME_MAXDIFF} s" >&2
	echo "ERROR: FERS TS=${FERSTS}, nearest digi TS=${DIGI_TS}" >&2
	exit 3
    fi
    CORRELATION_INITIAL_OFFSET=$(( DIGI_TS - FERSTS ))
    CORRELATION_INITIAL_RANGE=100
    echo "TRACE: Nearest digitizer run timestamp: $DIGI_TS. Difference: $DIGI_TS_DIFF s. New correlation offset guess ${CORRELATION_INITIAL_OFFSET}, range ${CORRELATION_INITIAL_RANGE}. It will be overwritten - just for info" >&2
fi
#overwrite the digitizer with a manual number if given as a parameter
if [[ $# -ge 2 ]] ; then
    if [[ "$2" =~ ^[0-9]+$ ]]; then
	echo "Digitizer run number forced as a parameter. Using $2 instead of ${DIGI_TS}"
	DIGI_TS=$2
	DIGI_TS_DIFF=$(( DIGI_TS - FERSTS ))
	echo "New time difference is ${DIGI_TS_DIFF}"
    else
	echo "second parameter $2 could not be parsed"
	exit 4
    fi
fi

###################################################################################
# waveform samples for the first N (200? events)
###################################################################################
if [[ ! -f "${QCPATH}/waveform_sample_run${RUN}.png" ]] ; then
    echo "Waveforms for first ~500 waveforms. Run ${RUN}, DIGI ${DIGI_TS}"
    gnuplot <<EOF
set terminal pngcairo size 1920,1800;
set output "${QCPATH}/waveform_sample_run${RUN}.png";
#set xlabel "digitizer time [bins]";
#set ylabel "raw digi ADC [bin]";
set multiplot layout 4,2 title "run ${RUN}"
set grid;
set xtics 100;
set xrange [0:1030];
do for [i=0:7] {
   plot '<grep -h -E "^[0-9]" ${DIGIPATH}/${DIGI_TS}/wave0000000001_0_'.i.'.txt' every ::::504000 u (int(\$0)%1030==1029?NaN:int(\$0)%1030):1 w l lc "blue" lw 0.2 title "channel ".i ;
}
unset multiplot; unset output;
EOF
else
    echo "waveforms already created before. Run ${RUN}, DIGI ${DIGI_TS}. Skipping"
fi

###################################################################################
# digitizer waveform extraction - update search parameters from analysis_overrides.sh 
###################################################################################
load_baseline() { #overrides values from a file
    local file="$1"
    local line key value
    while IFS= read -r line || [[ -n "$line" ]]; do
        # Remove possible Windows CR
        line=${line%$'\r'}
        # Skip empty lines and comments
        [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
        # Only allow NAME=integer, with no command execution
        if [[ "$line" =~ ^(DIGI_BASE_FROM|DIGI_BASE_TO|DIGI_PEAK_FROM|DIGI_PEAK_TO|CORRELATION_INITIAL_OFFSET|CORRELATION_INITIAL_RANGE)=([\-0-9\.]+)$ ]]; then
            key="${BASH_REMATCH[1]}"
            value="${BASH_REMATCH[2]}"
            case "$key" in
                DIGI_BASE_FROM) DIGI_BASE_FROM="$value" ;;
                DIGI_BASE_TO)   DIGI_BASE_TO="$value" ;;
                DIGI_PEAK_FROM) DIGI_PEAK_FROM="$value" ;;
                DIGI_PEAK_TO)   DIGI_PEAK_TO="$value" ;;
		CORRELATION_INITIAL_OFFSET) CORRELATION_INITIAL_OFFSET="$value" ;;
		CORRELATION_INITIAL_RANGE) CORRELATION_INITIAL_RANGE="$value" ;;
            esac
	elif [[ "$line" =~ ^DIGI_CHANNELS_TO_ANALYZE=\(([0-9]+([[:space:]]+[0-9]+)*)\)$ ]]; then
            read -r -a DIGI_CHANNELS_TO_ANALYZE <<< "${BASH_REMATCH[1]}"
	else
            echo "Error: invalid line in $file:" >&2
            echo "  $line" >&2
            return 1
        fi
    done < "$file"
}

if [[ -f  "${QCPATH}/analysis_overrides.sh" ]] ; then
    load_baseline "${QCPATH}/analysis_overrides.sh" || exit 1
fi

###################################################################################
# process digitizer events per channel - one line per event
###################################################################################
for CH in "${DIGI_CHANNELS_TO_ANALYZE[@]}"; do
    if [[ ! -f "${QCPATH}/digi_run${RUN}_events_ch${CH}.txt" ]] ; then
	echo "processing digitizer events for channel ${CH}, run ${RUN}"
	files=( "${DIGIPATH}/${DIGI_TS}"/wave*"${CH}".txt )
	if [[ ! -e "${files[0]}" ]]; then
	    echo "ERROR: no waveform files found for channel ${CH}" >&2
	    exit 4
	fi
	cat "${files[@]}" | python3 "${TOOLPATH}/digi-event-analysis.py" \
				    -m filtfilt \
				    --peak-from ${DIGI_PEAK_FROM} \
				    --peak-to ${DIGI_PEAK_TO} \
				    --baseline-from ${DIGI_BASE_FROM} \
				    --baseline-to ${DIGI_BASE_TO} \
				    >"${QCPATH}/digi_run${RUN}_events_ch${CH}.txt"
    else
	echo "digitizer events for channel ${CH}, run ${RUN} already created. Skipping"
    fi   
done

###################################################################################
# refinement of the offset guess
###################################################################################
# another method for a offset guess - taking the first event. This works only if the first event is the same, or within the range
FERSTSLINE=$(${DT5202CONVERTER} -w ${FERSFILE} --to_trigger 10 --board_id "${FERS_BOARDS[0]}" --print_txt_event| grep -m1 -E '^[0-9]+[[:space:]].*#event$')
read -r _ FERSFIRSTEVENTTS _ <<<${FERSTSLINE}
#echo "DEBUG: FERSFIRSTEVENTTS=${FERSFIRSTEVENTTS}"
DIGITSLINE=$(grep -m1 -E '^[0-9]+[[:space:]]' "${QCPATH}/digi_run${RUN}_events_ch${DIGI_CHANNELS_TO_ANALYZE[0]}.txt")
read -r _ _ _ _ _ DIGIFIRSTTS _ <<<${DIGITSLINE}
#echo "DEBUG DIGIFIRSTTS=${DIGIFIRSTTS}"
GUESS_OFFSET_FIRSTEVENT=$(echo "$FERSFIRSTEVENTTS - 0.000000008 * $DIGIFIRSTTS" | bc -l)
echo "first TS of FERS=${FERSFIRSTEVENTTS} s. First digi TS=${DIGIFIRSTTS} bin. offset guess=${GUESS_OFFSET_FIRSTEVENT} s"
CORRELATION_INITIAL_OFFSET=${GUESS_OFFSET_FIRSTEVENT}
CORRELATION_INITIAL_RANGE=1
# we have to rewrite the offset if overrides are provided
if [[ -f  "${QCPATH}/analysis_overrides.sh" ]] ; then
    load_baseline "${QCPATH}/analysis_overrides.sh" || exit 1
fi

###################################################################################
# correlation
###################################################################################
OFFSET=0
for CORR in "${CORRELATIONS[@]}"; do
    IFS=: read -r CH BOARD EXPR <<< "$CORR"
    MERGEFILE="${QCPATH}/merged_run${RUN}_digi${CH}_fers${BOARD}_${EXPR}_events.txt"
    if [[ ! -f "${MERGEFILE}" ]] ; then
	echo "merging of events by HW clock. Run ${RUN} using digitizer ch ${CH}, fers board ${BOARD} ${EXPR}" >&2
	python3 ${TOOLPATH}/event_synctime_merge.py \
		--file1 "${QCPATH}/digi_run${RUN}_events_ch${CH}.txt" \
		--col1 6 \
		--file2 "${QCPATH}/fers_run${RUN}_events_board${BOARD}_${EXPR}.txt" \
		--col2 2 \
		--scale1 8 \
		--scale2 1000000000 \
		--iteration-overlap 0.51 \
		--initial-range ${CORRELATION_INITIAL_RANGE} \
		--initial-offset ${CORRELATION_INITIAL_OFFSET} \
		>"${MERGEFILE}"
    else
	echo "merged file merged_run${RUN}_digi${CH}_fers${BOARD}_${EXPR}_events.txt exists. skipping" >&2
    fi

    #extract the offset and update the correlation initial range
    while IFS= read -r line; do
	if [[ "$line" =~ ^#Summary[[:space:]]+offset=([^[:space:]]+)[[:space:]]+drop1=([0-9]+)[[:space:]]+drop2=([0-9]+)[[:space:]]+synchronized=([0-9]+) ]]; then
            summary_line="$line"
            OFFSET="${BASH_REMATCH[1]}"
            DROP1="${BASH_REMATCH[2]}"
            DROP2="${BASH_REMATCH[3]}"
            SYNCHRONIZED="${BASH_REMATCH[4]}"
	    echo "#Summary offset=${OFFSET} drop1=${DROP1} drop2=${DROP2} synchronized=${SYNCHRONIZED}" >&2
	    if (( SYNCHRONIZED * CORRELATION_ACCEPT_FACTOR > ( DROP1 > DROP2 ? DROP1 : DROP2 ) )); then
		echo "correlation accepted. fixing the offset parameter"  >&2
		CORRELATION_INITIAL_OFFSET=${OFFSET}
		CORRELATION_INITIAL_RANGE=0.000000001 
	    else
		echo "correlation rejected - only ${SYNCHRONIZED} events out of $(( DROP1 > DROP2 ? DROP1 : DROP2 )) " >&2
		echo "This can be caused by providing a wrong guess for the initial offset. If this happens too oftern, consider changing the CORRELATION_INITIAL_RANGE to a higher value. For single run overrides, create a file analysis_overrides.sh un the runXXXX directory and put there two lines"
		echo "CORRELATION_INITIAL_OFFSET=<GUESS_OFFSET>"
		echo "CORRELATION_INITIAL_RANGE=100"
	    fi
	fi
    done < <(tail -n 5 "${MERGEFILE}")

    ####################################################
    # ADC correlation plot
    ####################################################
    if [[ ! -f "${QCPATH}/merged_run${RUN}_digi${CH}_fers${BOARD}_${EXPR}_ADC_correlation.png" ]] ; then
	echo "Plotting ${QCPATH}/merged_run${RUN}_digi${CH}_fers${BOARD}_${EXPR}_ADC_correlation.png" >&2
	gnuplot <<EOF
set terminal pngcairo enhanced size 800,600
set output "${QCPATH}/merged_run${RUN}_digi${CH}_fers${BOARD}_${EXPR}_ADC_correlation.png"
set xlabel "Digi ADC [bin]";
set ylabel "FERS HG ADC [bin]"
set title "run ${RUN}, ADC correlation FERS${BOARD} ${EXPR}";
set cblabel "run time [s]"
set grid
plot "${MERGEFILE}" u 12:21:16 palette ti sprintf("offset %.9f s",${OFFSET})
EOF
    else
        echo "Exists: ${QCPATH}/merged_run${RUN}_digi${CH}_fers${BOARD}_${EXPR}_ADC_correlation.png , skipping." >&2
    fi

    ####################################################
    # Time difference plot
    ####################################################
    if [[ ! -f "${QCPATH}/merged_run${RUN}_digi${CH}_fers${BOARD}_${EXPR}_events_time_difference.png" ]] ; then
	echo "Plotting ${QCPATH}/merged_run${RUN}_digi${CH}_fers${BOARD}_${EXPR}_events_time_difference.png" >&2
	gnuplot <<EOF
set terminal pngcairo enhanced size 800,600
set output "${QCPATH}/merged_run${RUN}_digi${CH}_fers${BOARD}_${EXPR}_events_time_difference.png"
set ylabel "time difference (digi - FERS) [ns]"
set title "run ${RUN}, time difference between digitizer and FERS board${BOARD}" ;
set xrange [0:*]
a=0.000000008;
b=${OFFSET};
set grid
f(x)=a*x+b;
set xlabel "FERS board ${BOARD} event number"
plot "${MERGEFILE}" u 19:((f(\$6)-\$16)*1000000000) title sprintf("offset %.9f s",${OFFSET}) w p lw 0.3
set xlabel "run time [s]";
EOF
    else
        echo "Exists: ${QCPATH}/merged_run${RUN}_digi${CH}_fers${BOARD}_${EXPR}_events_time_difference.png , skipping." >&2
    fi
done 

###################################################################################
# Time evolution plot
###################################################################################
if [[ ! -f "${QCPATH}/time_series_run${RUN}.png" ]] ; then
    echo "plotting ${QCPATH}/time_series_run${RUN}.png"
    gnuplot <<EOF
set terminal pngcairo enhanced size 1600,480
set output "${QCPATH}/time_series_run${RUN}.png"
set xlabel "time (offset applied) [s]"
set title sprintf("run ${RUN}, data time structure. offset %.9f s",${OFFSET});
set ylabel "A.U. - sub-second remainer (shifted)"
set yrange [0:3]
set grid
# Horizontal separators, full plot width
set arrow 1 from graph 0, first 1 to graph 1, first 1 nohead lc rgb "black" lw 1 back
set arrow 2 from graph 0, first 2 to graph 1, first 2 nohead lc rgb "black" lw 1 back
ptsize=0.3
# for smaller point size, use "w d", for bigger points use "w p pt 1 ps ptsize
# plot "${QCPATH}/digi_run${RUN}_events_ch${DIGI_CHANNELS_TO_ANALYZE[0]}.txt"  u (\$6*0.000000008+${OFFSET}):(0.000000008*\$6+${OFFSET}-floor(0.000000008*\$6+${OFFSET})) w p pt 1 ps ptsize title "digitizer", \
#  "${QCPATH}/fers_run${RUN}_events_board0_hg0.txt" u 2:(1+\$2-floor(\$2)) w p pt 1 ps ptsize title "FERS board 0", \
#  "${QCPATH}/fers_run${RUN}_events_board0_hg0.txt" u 2:(2+\$2-floor(\$2)) w p pt 1 ps ptsize title "FERS board 1"
plot "${QCPATH}/digi_run${RUN}_events_ch${DIGI_CHANNELS_TO_ANALYZE[0]}.txt"  u (\$6*0.000000008+${OFFSET}):(0.000000008*\$6+${OFFSET}-floor(0.000000008*\$6+${OFFSET})) w d title "digitizer", \
 "${QCPATH}/fers_run${RUN}_events_board0_hg0.txt" u 2:(1+\$2-floor(\$2)) w d title "FERS board 0", \
 "${QCPATH}/fers_run${RUN}_events_board0_hg0.txt" u 2:(2+\$2-floor(\$2)) w d title "FERS board 1"
EOF
else
    echo "Exists: ${QCPATH}/time_series_run${RUN}.png . skipping"
fi

 
