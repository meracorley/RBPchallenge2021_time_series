#!/bin/bash

PATH_TO_INPUT_DATA_SET="./input_data_set/hnrnpa2b1_binding_sites_fshape"
PATH_TO_RESULTS_SET="./results/hnrnpa2b1_binding_sites_fshape"
EXPECTED_MOTIF_LENGTH=13

show_help() {                                 # Function: Print a help message.
  echo "Usage: $0 [ -i PATH_TO_INPUT_DATA_SET ] [ -r PATH_TO_RESULTS_SET ] [ -l EXPECTED_MOTIF_LENGTH ]" 1>&2 
  exit 1
}

while getopts ":i:r:l:" options; do
  case "${options}" in
    i)                
      PATH_TO_INPUT_DATA_SET=${OPTARG}  
      ;;    
	r)                
      PATH_TO_RESULTS_SET=${OPTARG}  
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
elif [[ ! -d $PATH_TO_RESULTS_SET ]]; then
	echo "${PATH_TO_RESULTS_SET} is not directory!"
    exit 2
elif [[ $EXPECTED_MOTIF_LENGTH > 5 ]] && [[ $EXPECTED_MOTIF_LENGTH < 16 ]]; then
	echo "${EXPECTED_MOTIF_LENGTH} should be integer in range <6;15>!"
    exit 2
fi
mkdir $PATH_TO_RESULTS_SET/orig
for f in ${PATH_TO_INPUT_DATA_SET}/*.txt; do 
	file_name=$(basename "${f}")
	awk 'BEGIN {print "Reactivity,Sequence";} {printf("%s,%s\n",$1,$2);}' < "$f" > "${PATH_TO_RESULTS_SET}/orig/${file_name%.txt}.csv"
done
mkdir $PATH_TO_RESULTS_SET/$EXPECTED_MOTIF_LENGTH
python silence.py -i $PATH_TO_RESULTS_SET/orig -r $PATH_TO_RESULTS_SET/$EXPECTED_MOTIF_LENGTH -l $EXPECTED_MOTIF_LENGTH
./clear_data.sh -i $PATH_TO_RESULTS_SET/$EXPECTED_MOTIF_LENGTH -l $EXPECTED_MOTIF_LENGTH
exit 0
