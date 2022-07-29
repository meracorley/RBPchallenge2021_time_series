# find-conserved-motifs

A tool represented by a set of scripts looks for conserved motif(s) in a set of SHAPE/fSHAPE data that can be treated as time series data.

# How it works

Generally, we first preprocess the input data and next, searching for conserved motifs in the considered data set.

## Data preparation

All data files were transformed into CSV format to simplify analysis. Next, in every considered data file, we identify and silence all data points that we are sure cannot be part of the searching RNA-protein binding site. Next, we remove all undefined data points at both ends and shortened continuous sequences of undefined values into a single representative.

## Searching for conserved motifs

We identify the most central motif and its closest matches in remained data files. At the output, we provide: 
* a conserved motif identified for the particular length with the lowest radius,
*  hierarchical clustering dendrogram of aligned motifs,
* alignment of aligned motifs,
* a list of all motifs found for the particular data set.

Every motif is described by its DNA sequence, and its localization. 

# Methods

The [STUMPY library](https://github.com/TDAmeritrade/stumpy) provides algorithms for finding: (a) the most central motif among the subsequences with the best radius [**ostinato**] and (b) the closest match to the central one based on the distance profile [**mass**]. The central motif identified in a single file is the one with the smallest mean distance to nearest neighbors in remained time series files. In our case, we run the ostinato algorithm multiple times with different window sizes at the input. It's worth noticing that the mass algorithm always finds the solution, even if it is distant from the query. We also used other libraries like numpy and matplotlib to simplify data transformation and visualization. The crucial advantage of the proposed solutions is processing efficiency. The packages we used are easy to use, well-documented, open-source projects published on GitHub.

# Input

All input data files meet the same format in which every line consists of two required columns:

```
fSHAPE          Sequence
```

For example:

```
4.22243609882	T
2.33772397905	T
4.17355287021	G
...
```

# Requirements

* Linux, e.g., Ubuntu 22.04 LTS,
* Anaconda,
* Python 3.8.

# Installation with Anaconda

```
conda create --name rbp --file requirements.txt python=3.8
conda activate rbp
cd RBPchallenge2021_time_series/find-conserved-motifs 
./run.sh [ -i PATH_TO_INPUT_DATA_SET ] [ -l EXPECTED_MOTIF_LENGTH ]
```

# Available scripts manual

```
>./run.sh -h
Usage: ./run.sh [ -i PATH_TO_INPUT_DATA_SET ] [ -l EXPECTED_MOTIF_LENGTH ]
```

```
>./transform.sh -h
Usage: ./transform.sh [ -i PATH_TO_INPUT_DATA_SET ] [ -r PATH_TO_RESULTS_SET ] [ -l EXPECTED_MOTIF_LENGTH ]
```

```
>silence.py 
silence.py -i <input_data_path> (required) -r <results_path> (required) -l <expected_motif_length>  (required)
```

```
>find-conserved-motifs.py 
find-conserved-motifs.py -i <input_data_path> (required) -r <results_path> (required) -l <expected_motif_length>  (required)
```

# Example usage scenario

Input data set including fSHAPE values on known A2B1 binding sites is stored in the following directory **input_data_set/hnrnpa2b1_binding_sites_fshape**.

```
cd RBPchallenge2021_time_series/find-conserved-motifs
./run.sh -i ./input_data_set/hnrnpa2b1_binding_sites_fshape -l 13
```

Output:

```
Transform and clean data...
Input data path: RBPchallenge2021_time_series/find-conserved-motifs/results/hnrnpa2b1_binding_sites_fshape/orig
Expected motif length: 13
Results path: RBPchallenge2021_time_series/find-conserved-motifs/results/hnrnpa2b1_binding_sites_fshape/13
Done.
Find conserved motifs...
Input data path: RBPchallenge2021_time_series/find-conserved-motifs/results/hnrnpa2b1_binding_sites_fshape/13
Expected motif length: 13
Results path: RBPchallenge2021_time_series/find-conserved-motifs/results/hnrnpa2b1_binding_sites_fshape
Done.
```

Obtained results for the particular data set are stored in the following directory **results/hnrnpa2b1_binding_sites_fshape**.
