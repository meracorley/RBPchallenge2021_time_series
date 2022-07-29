#!/usr/bin/python

import sys
import pandas as pd
import numpy as np
import getopt
import os

def read_config(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:r:l:",["input_data_path=","results_path=","expected_motif_length="])
    except getopt.GetoptError:
        print('silence.py -i <input_data_path> (required) -r <results_path> (required) -l <expected_motif_length>  (required)')
        sys.exit(1)
    results_path = "./results/hnrnpa2b1_binding_sites_fshape"
    expected_motif_length = 13
    required_arguments_count = 0
    for opt, arg in opts:
        if opt == '-h':
            print('silence.py -i <input_data_path> (required) -r <results_path> (required) -l <expected_motif_length>  (required)')
            sys.exit(0)
        elif opt in ("-i", "--input_data_path"):
            input_data_path = arg
            required_arguments_count = required_arguments_count + 1
        elif opt in ("-r", "--results_path"):
            results_path = arg
            required_arguments_count = required_arguments_count + 1
        elif opt in ("-l", "--expected_motif_length"):
            expected_motif_length = int(arg)
            required_arguments_count = required_arguments_count + 1
    if required_arguments_count != 3:
        print('silence.py -i <input_data_path> (required) -r <results_path> (required) -l <expected_motif_length>  (required)')
        sys.exit(1)
    return input_data_path, results_path, expected_motif_length
    
def process(input_file_path, output_file_path, expected_motif_length):
    results_df = pd.read_csv(input_file_path)
    values = results_df["Reactivity"].iloc[:].values
    mask = [0] * len(values)

    for i in range(len(values)):
        if abs(values[i]) > 1.0:
            for j in range(2 * expected_motif_length - 1):
                window_idx = i - expected_motif_length + 1 + j
                if window_idx < len(values) and mask[window_idx] == 0:
                    mask[window_idx] = 1;
    count = 0;
    for i in range(len(values)):
        if not np.isnan(values[i]) and mask[i] == 0:
            results_df["Reactivity"].iat[i] = float("nan")
            count = count + 1;

    results_df.to_csv(output_file_path,index=False, na_rep='NA')
    
def filter_useless_data_points(argv):
    input_data_path, results_path, expected_motif_length = read_config(argv)
    print('Input data path: {}'.format(input_data_path))
    print('Expected motif length: {}'.format(expected_motif_length))
    print('Results path: {}'.format(results_path))
    
    for filename in os.listdir(input_data_path):
        current_file = os.path.join(input_data_path,filename)
        if current_file.endswith(".csv") and os.path.isfile(current_file):
            process(current_file, os.path.join(results_path,filename), expected_motif_length)

def main(argv):
    filter_useless_data_points(argv)

if __name__ == "__main__":
    main(sys.argv[1:])