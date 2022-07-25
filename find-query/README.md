# find-query

This script looks for a given SHAPE or fSHAPE matrix profile in a set of input files.

# Input

The query and all other input files have the same file format in which every line has three columns:

```
fSHAPE          Sequence  SHAPE
```

For example:

```
0.531188539549  T         0.822363
0.279392225779  T         0.714092
1.41084255243   G         1.46544
...
```

The `Sequence` and `SHAPE` columns are optional. If absent, their default values will be `N` and `NaN` respectively.

# Usage

```
./find-query --query PATH-TO-QUERY-FILE PATH-TO-INPUT-1 PATH-TO-INPUT-2 ...
```

Example:

```
./find-query --query fshape-true-pattern.txt fshape-inputs/*
```

# How it works

1. Compute [matrix profile](https://pypi.org/project/matrixprofile/) for every input with the window size equal to the length of the query
2. Discover up to 10 motifs in every input file
3. Filter out any motifs with `NaN` (Not-a-Number) value inside
4. Sort the motifs according to the Z-normalized Euclidean distance

# Output

## Table

The script creates two CSV files with the following columns:

- *Sample*: name of the input file the motif comes from
- *Range*: indicies of nucleotides in the input file where the motif was found
- *Sequence*: sequence of the motif
- *Z-normalized*: the Z-normalized Euclidean distance of the motif to the query
- *Distance*: a regular euclidean distance of the motif to the query
- *Sequence-Score*: an integer representing similarity (higher is better) of the motif sequence to query sequence
- *fSHAPE-n*: fSHAPE value for every nucleotide in the motif
- *SHAPE-n*: SHAPE value for every nucleotide in the motif

The `output.csv` contains all results, while `output-filtered.csv` may have some of the motifs filtered out. In the latter case, motif is discarded if there is a strong discrepancy for any of its nucleotides (i.e. if the query requests it to be a strong signal, but the motif does not show it). 

The result usually requires manual curation. The table is sorted according to the *Z-normalized* column, but it is just a starting point and user should take into account *Distance* and *Sequence-Score* as well when selecting the matches.

## Plots

The script creates two plots:

- `motifs-highlighted.pdf` is a plot of the 10 best input files with the motif highlighted
- `motifs-only.pdf` is a plot of the 10 best motifs only (the remaining context is not shown)