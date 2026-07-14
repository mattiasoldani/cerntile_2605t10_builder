#!/bin/bash
if ! [[ "$1" =~ ^[0-9]+$ ]]; then
  printf "Not a valid parameter.\nUsage $0 <DIGITIZER_TS(AKA run number)>\n"
  exit 1
fi
TS=$1
#OUT_DIR=/eos/home-j/jkvas/DRD6/newtile/tools/dt5202_raw_utils/raw2txt/dt5202txt

RAWDIR=/eos/experiment/newtile/beamtests/26_05_t10/digi_raw/${TS}
#pwd
gnuplot <<EOF
set terminal pngcairo size 800,600;
set output "waveform_sample${TS}.png";
set xlabel "digitizer time [bins]";
set ylabel "raw digi ADC [bin]";
set title "Digitizer run ${TS}"
set grid;
set xtics 100;
set xrange [0:1030];
plot for [i=0:2] '<grep -h -E "^[0-9]" $RAWDIR/wave0000000001_0_'.i.'.txt' every ::::200000 u (int(\$0)%1030==1029?NaN:int(\$0)%1030):1 w l lw 0.1 title "channel ".i;
EOF

#display waveform_sample.png
