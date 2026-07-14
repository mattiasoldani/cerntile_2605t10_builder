#!/bin/bash
if ! [[ "$1" =~ ^[0-9]+$ ]]; then
  printf "Not a valid parameter.\nUsage $0 <RUN_NUMBER>\n"
  exit 1
fi


FERS_RAW=/eos/experiment/newtile/beamtests/26_05_t10/fers_daq/bin/DataFiles
RUN=$1
DT5202CONVERTER=/eos/home-j/jkvas/DRD6/newtile/tools/dt5202_raw_utils/raw2txt/dt5202txt

mkdir -p "run${RUN}" 
mapfile -t files < <(printf '%s\n' "${FERS_RAW}"/Run"${RUN}".*_list.dat | sort -V)
echo "processing Run $1 with files: ${files[@]}"

for BOARD in 0 1 ; do
    if [ ! -f run${RUN}/hist${RUN}_board${BOARD}.txt ] ; then
	rm -f Run${RUN}.dat.fifo
	#echo "sending data to FIFO  Run${RUN}.dat.fifo"
	mkfifo Run${RUN}.dat.fifo
	cat  "${files[@]}" >Run${RUN}.dat.fifo &
	echo "generating histogram for run ${RUN} board ${BOARD}"
	${DT5202CONVERTER} -w  Run${RUN}.dat.fifo \
			   --print_hg_histogram \
			   --print_lg_histogram \
			   --max_bin=8192 \
			   --board_id=${BOARD} \
			   >run${RUN}/hist${RUN}_board${BOARD}.txt
	#echo "clearing residual fifo using dd"    
	dd if=Run${RUN}.dat.fifo of=/dev/null iflag=nonblock status=none
	wait
	#echo "processing from FIFO is finished"
	rm -f Run${RUN}.dat.fifo
    fi

    if [ ! -f run${RUN}/dt5202_events_board${BOARD}.txt ] ; then
	rm -f Run${RUN}.dat.fifo
	# echo "sending data to FIFO  Run${RUN}.dat.fifo" 
	mkfifo Run${RUN}.dat.fifo
	cat  "${files[@]}" >Run${RUN}.dat.fifo &
	echo "generating event list run ${RUN} board ${BOARD}"
	${DT5202CONVERTER} -w  Run${RUN}.dat.fifo \
			   --channel=0 \
			   --board_id=${BOARD} \
			   --custom_expr="hg0" \
			   --print_txt_events | grep "#" \
						     >run${RUN}/dt5202_events_board${BOARD}.txt
	# echo "clearing residual fifo using dd"    
	dd if=Run${RUN}.dat.fifo of=/dev/null iflag=nonblock status=none
	wait
	# echo "processing from FIFO is finished"
	rm -f Run${RUN}.dat.fifo
    fi
done
