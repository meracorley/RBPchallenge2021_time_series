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
fi
if [[ ! -d results ]]; then
	mkdir results
fi
INPUT_DATA_SET_PATH=$(realpath $PATH_TO_INPUT_DATA_SET)
RESULTS_PATH=$(realpath "results")
INPUT_NAME=$(basename $PATH_TO_INPUT_DATA_SET)
if [[ -d $RESULTS_PATH/$INPUT_NAME ]]; then
	rm -r $RESULTS_PATH/$INPUT_NAME
fi
mkdir $RESULTS_PATH/$INPUT_NAME
echo "Transform and clean data..."
./transform.sh -i $INPUT_DATA_SET_PATH -r $RESULTS_PATH/$INPUT_NAME -l $EXPECTED_MOTIF_LENGTH
echo "Done."
echo "Find conserved motifs..."
python find-conserved-motifs.py -i $RESULTS_PATH/$INPUT_NAME/$EXPECTED_MOTIF_LENGTH -r $RESULTS_PATH/$INPUT_NAME -l 13
echo "Done."
exit 0