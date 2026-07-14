#!/bin/bash
if ! [[ "$1" =~ ^[0-9]+$ ]]; then
  printf "Not a valid parameter.\nUsage $0 <RUN_NUMBER>\n"
  exit 1
fi

RUN=$1

for OFFSET in `seq -20 20` ; do
    python3 ./tools/event_time_merge.py          \
	    --file1 run${RUN}/digi-events_ch0.txt --col1 6            \
	    --file2 run${RUN}/dt5202_events_board0.txt --col2 2       \
	    --scale1 8.0            \
	    --scale2 1000000000 \
            --jitter 500            \
	    --max-drift 1 \
	    --event-offset ${OFFSET}             \
	    >run${RUN}/merged_time_o${OFFSET}.txt

    gnuplot <<EOF
run=${RUN};
set fit quiet;
set terminal png enhanced size 1600, 600 ;
set output "run${RUN}/run${RUN}_sync_o${OFFSET}.png" ;
set multiplot layout 1,2
a=0.000000008;
b=-1;
set grid
f(x)=a*x+b;
fit f(x) 'run'.run.'/merged_time_o${OFFSET}.txt' u 6:16 via a,b
set xlabel "run time [s]";
set ylabel "time difference (digi - FERS) [ns]" 
set title "run ".run.", time residuals" ;
plot 'run'.run.'/merged_time_o${OFFSET}.txt' u 16:((f(\$6)-\$16)*1000000000) title "offset ${OFFSET}"
set xlabel "Digi ADC, dummy pedestal subtracted [bin]";
set ylabel "FERS HG ADC [bin]"
set title "run ".run.", ADC correlation";
set cblabel "run time [s]"
plot 'run'.run.'/merged_time_o${OFFSET}.txt' u 12:21:16 palette ti "offset ${OFFSET}"

EOF

    # if [ ! ${OFFSET} -eq 0 ] ; then
    #    rm run${RUN}/merged_time_o${OFFSET}.txt
    # fi

done

