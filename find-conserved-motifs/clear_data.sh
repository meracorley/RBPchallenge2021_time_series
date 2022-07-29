#!/bin/bash

PATH_TO_INPUT_DATA_SET="./input_data_set/hnrnpa2b1_binding_sites_fshape"
EXPECTED_MOTIF_LENGTH=13

show_help() {                                 # Function: Print a help message.
  echo "Usage: $0 [ -i PATH_TO_INPUT_DATA_SET ] [ -l EXPECTED_MOTIF_LENGTH ]" 1>&2 
  exit 1
}

while getopts ":i:l:" options; do
  case "${options}" in
    i)                
      PATH_TO_INPUT_DATA_SET=${OPTARG}  
      ;;    
    l)                
      EXPECTED_MOTIF_LENGTH=${OPTARG}  
      ;;
    :)
      echo "Error: -${OPTARG} requires an argument."
      show_help
      ;;
    *)         
      show_help
      ;;
  esac
done

if [[ ! -d $PATH_TO_INPUT_DATA_SET ]]; then
	echo "${PATH_TO_INPUT_DATA_SET} is not directory!"
    exit 2
elif [[ $EXPECTED_MOTIF_LENGTH > 5 ]] && [[ $EXPECTED_MOTIF_LENGTH < 16 ]]; then
	echo "${EXPECTED_MOTIF_LENGTH} should be integer in range <6;15>!"
    exit 2
fi
cp remove_na.awk $PATH_TO_INPUT_DATA_SET
for f in ${PATH_TO_INPUT_DATA_SET}/*.csv; do 
	awk -v expected_motif_length=${EXPECTED_MOTIF_LENGTH} -f remove_na.awk < "$f" > "${f%.csv}.csvv" 
done
rm $PATH_TO_INPUT_DATA_SET/remove_na.awk
for f in ${PATH_TO_INPUT_DATA_SET}/*.csvv; do 
	mv -- "$f" "${f%.csvv}.csv"
done
for f in ${PATH_TO_INPUT_DATA_SET}/*.csv; do 
	lines_count=$(cat ${f} | wc -l)
	if [[ $lines_count == 1 ]]; then
		rm $f
	fi
done

