#!/bin/bash
if ! [[ "$1" =~ ^[0-9]+$ ]]; then
    printf "Not valid parameter (1st position).\nUsage $0 <DIGITIZER_TS(AKA run number)> <FERS_RUN_NR>\n"
    exit 1
fi
if ! [[ "$2" =~ ^[0-9]+$ ]]; then
    printf "Not valid parameter (2nd position).\nUsage $0 <DIGITIZER_TS(AKA run number)> <FERS_RUN_NR>\n"
    exit 1
fi


TS=$1
RUN=$2
OUT_DIR=/eos/home-j/jkvas/DRD6/newtile/run${RUN}
RAWDIR=/eos/experiment/newtile/beamtests/26_05_t10/digi_raw/${TS}

#convert digitizer events
if [ ! -f run${RUN}/digi-events_ch0.txt ] ; then
    echo "converting digitizer TS ${TS} RUN ${RUN}"
    bash process_digi.sh ${TS} ${RUN}
fi

#convert FERS raw
if [ ! -f run${RUN}/dt5202_events_board0.txt ] ; then
    echo "generating histogram and event for FERS run ${RUN}"
    bash make_fers_events.sh ${RUN}
fi


#waveform sample
echo "generating waveform wample"
gnuplot <<EOF
set terminal pngcairo size 800,600;
set output "run${RUN}/waveform_sample${TS}.png";
set xlabel "digitizer time [bins]";
set ylabel "raw digi ADC [bin]";
set title "Digitizer run ${TS}"
set grid;
set xtics 100;
set xrange [0:1030];
plot for [i=0:2] '<grep -h -E "^[0-9]" $RAWDIR/wave0000000001_0_'.i.'.txt' every ::::200000 u (int(\$0)%1030==1029?NaN:int(\$0)%1030):1 w l lw 0.1 title "channel ".i;
EOF

echo "merging event for run ${RUN}"
bash merge_time.sh ${RUN}
