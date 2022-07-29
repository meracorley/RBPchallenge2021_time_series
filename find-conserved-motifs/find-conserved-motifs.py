#!/usr/bin/python

import stumpy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import cycle, combinations
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.special import comb
import glob
import sys
import os
import getopt

def find(data_path, results_path, m):
    data_files_pattern = '{}/*.csv'.format(data_path)
    data_files = glob.glob(data_files_pattern)

    Ts = [None] * len(data_files)
    seq = [None] * len(data_files)

    mmin = 1000.0
    mmax = -1000.0
    for i, s in enumerate(data_files):
      sample_df = pd.read_csv(s)
      Ts[i] = sample_df["Reactivity"].iloc[:].values
      seq[i] = sample_df["Sequence"].iloc[:].values
      current_min = np.nanargmin(Ts[i])
      if (Ts[i][current_min] < mmin):
        mmin = Ts[i][current_min]
      current_max = np.nanargmax(Ts[i])
      if (Ts[i][current_max] > mmax):
        mmax = Ts[i][current_max]
    mmax = mmax + 0.2
    mmin = mmin - 0.2

    conserved_motifs_list = []
    
    radius, Ts_idx, subseq_idx = stumpy.ostinato(Ts, m)
    
    tseq = '';
    for a in range(m):
        tseq = tseq + seq[Ts_idx][subseq_idx + a] 
    
    best_sample = data_files[Ts_idx].replace(os.path.join(data_path, ''),"").replace(".csv","")
    conserved_motifs_list.append(f'Lowest radius ({np.round(radius, 2)}) found in location {subseq_idx+1}-{subseq_idx+m+1} of data file {best_sample} (seed motif sequence: {tseq}).')
    
    seed_motif = Ts[Ts_idx][subseq_idx : subseq_idx + m]

    save_conserved_motif(seed_motif, results_path, m)
    nn = plot_motifs_alignment(plt, Ts, seq, Ts_idx, subseq_idx, seed_motif, data_files, m, data_path, results_path, conserved_motifs_list)
    save_conserved_motifs_list(results_path, conserved_motifs_list, m)
    plot_clustering_dendrogram(plt, Ts, Ts_idx, data_files, m, data_path, results_path, nn, seed_motif)
    plot_independent_motifs_matched_to_conserved_one(plt, Ts, data_path, data_files, m, results_path, nn, mmin, mmax)

def save_conserved_motif(seed_motif, results_path, m):
    np.savetxt(os.path.join(results_path, f'conserved-motif-{m}.csv'), np.asarray(seed_motif), delimiter=",")

def plot_motifs_alignment(plt, Ts, seq, Ts_idx, subseq_idx, seed_motif, data_files, m, data_path, results_path, conserved_motifs_list):
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[0] = 20
    fig_size[1] = 6
    plt.rcParams["figure.figsize"] = fig_size
    plt.rcParams['xtick.direction'] = 'out'

    plt.figure()
    x = np.linspace(0,1,m)
    nn = np.zeros(len(Ts), dtype=np.int64)
    nn[Ts_idx] = subseq_idx
    for i, e in enumerate(Ts):
        if i != Ts_idx:
            nn[i] = np.argmin(stumpy.core.mass(seed_motif, e))
            lw = 1
            label = None
            oseq = '';
            for a in range(m):
                val = seq[i][nn[i] + a]
                if isinstance(val, float) and np.isnan(val):
                    oseq = oseq + 'N'
                else:
                    oseq = oseq + val
            current_sample = data_files[i].replace(os.path.join(data_path, ''),"").replace(".csv","")            
            conserved_motifs_list.append(f'{current_sample}: {oseq} {nn[i]+1}-{nn[i]+m+1}')
        else:
            lw = 4
            label = 'Seed Motif'
        plt.plot(x, e[nn[i]:nn[i]+m], lw=lw, label=label)
    plt.title('All Motifs Alignment (fSHAPE)')
    plt.xlabel(f'Motif sequence (w = {m} nts)')
    plt.ylabel('Reactivity')
    plt.legend()
    plt.savefig(os.path.join(results_path, f'all-motifs-alignment-{m}.png'))
    plt.close()
    return nn

def save_conserved_motifs_list(results_path, conserved_motifs_list, m):
    with open(os.path.join(results_path, f'all-motifs-list-{m}.txt'), 'w') as fp:
        fp.write('\n'.join(conserved_motifs_list))

def get_map_of_data_points(data_files, data_path):
    data = {}
    for i, s in enumerate(data_files):
        data_file = s.replace(os.path.join(data_path, ''),"").replace(".csv","")    
        data[data_file] = i
    return data
    
def plot_clustering_dendrogram(plt, Ts, Ts_idx, data_files, m, data_path, results_path, nn, seed_motif):
    data = get_map_of_data_points(data_files, data_path)
    
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[1] = 4 * 10
    plt.rcParams["figure.figsize"] = fig_size

    plt.figure()
    fig, ax = plt.subplots()
    seed_name = data_files[Ts_idx].replace(os.path.join(data_path, ''),"").replace(".csv","")
    dp = np.zeros(int(comb(len(data), 2)))
    for i, a_c in enumerate(combinations(data.keys(), 2)):
        if seed_name == a_c[0]:
            other_idx = data[a_c[1]]
            other_motif = Ts[other_idx][nn[other_idx] : nn[other_idx] + m]
            dp[i] = stumpy.core.mass(seed_motif, other_motif)
        elif seed_name == a_c[1]:
            other_idx = data[a_c[0]]
            other_motif = Ts[other_idx][nn[other_idx] : nn[other_idx] + m]
            dp[i] = stumpy.core.mass(other_motif, seed_motif)
        else:
            first_idx = data[a_c[0]]
            first_motif = Ts[first_idx][nn[first_idx] : nn[first_idx] + m]
            second_idx = data[a_c[1]]
            second_motif = Ts[second_idx][nn[second_idx] : nn[second_idx] + m]
            dp[i] = stumpy.core.mass(first_motif, second_motif)
    Z = linkage(dp, optimal_ordering=True)
    dendrogram(Z, labels=[k for k in data.keys()], ax=ax)
    plt.ylabel('Z-Normalized Euclidean Distance')
    plt.title(f'Clustering (w = {m} nts)')
    plt.savefig(os.path.join(results_path, f'aligned-motifs-clustering-dendrogram-{m}.png'))
    plt.close()

def plot_independent_motifs_matched_to_conserved_one(plt, Ts, data_path, data_files, m, results_path, nn, mmin, mmax):
    fig_size = plt.rcParams["figure.figsize"]
    fig_size[1] = 10 * 20
    plt.rcParams["figure.figsize"] = fig_size

    plt.figure()
    ax = plot_vertical_signals(plt, Ts, data_path, data_files, m, mmin, mmax)
    for i in range(len(Ts)):
        y = ax[i].get_ylim()
        r = Rectangle((nn[i] / m, y[0]), 1, y[1]-y[0], alpha=0.3)
        ax[i].add_patch(r)
    plt.suptitle('fSHAPE', fontsize=14)
    plt.savefig(os.path.join(results_path, f'all-motifs-presented-independently-{m}.png'))
    plt.close()

def plot_vertical_signals(plt, Ts, data_path, data_files, m, mmin, mmax):
    fig, ax = plt.subplots(len(Ts), sharex=True, sharey=True)
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = cycle(prop_cycle.by_key()['color'])
    for i, e in enumerate(Ts):
        ax[i].plot(np.arange(0, len(e)) / m, e, color=next(colors), label=data_files[i].replace(os.path.join(data_path, ''),"").replace(".csv",""))
        ax[i].set_ylim((mmin, mmax))
        ax[i].legend()
        ax[i].set_ylabel('Reactivity')
    plt.subplots_adjust(hspace=0)
    plt.xlabel(f'Sequence split by motif length used (w = {m} nts)')
    return ax

def read_config(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:r:l:",["input_data_path=","results_path=","expected_motif_length="])
    except getopt.GetoptError:
        print('find-conserved-motifs.py -i <input_data_path> (required) -r <results_path> (required) -l <expected_motif_length>  (required)')
        sys.exit(1)
    results_path = "./results/hnrnpa2b1_binding_sites_fshape/13"
    expected_motif_length = 13
    required_arguments_count = 0
    for opt, arg in opts:
        if opt == '-h':
            print('find-conserved-motifs.py -i <input_data_path> (required) -r <results_path> (required) -l <expected_motif_length>  (required)')
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
        print('find-conserved-motifs.py -i <input_data_path> (required) -r <results_path> (required) -l <expected_motif_length>  (required)')
        sys.exit(1)
    return input_data_path, results_path, expected_motif_length

def find_conserved_motifs(argv):
    input_data_path, results_path, expected_motif_length = read_config(argv)
    print('Input data path: {}'.format(input_data_path))
    print('Expected motif length: {}'.format(expected_motif_length))
    print('Results path: {}'.format(results_path))
    find(input_data_path, results_path, expected_motif_length)

def main(argv):
    find_conserved_motifs(argv)

if __name__ == "__main__":
    main(sys.argv[1:])