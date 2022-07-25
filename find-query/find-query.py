#! /usr/bin/env python
import argparse
import csv
import math
import multiprocessing
import os
import random
import sys
from typing import List, Dict

import matrixprofile as mp
import numpy as np
from matplotlib import pyplot as plt


class Nucleotide:

    def __init__(self, fshape: float, base: str = 'N', shape: float = math.nan):
        self.fshape = fshape
        self.base = base
        self.shape = shape

    def __repr__(self):
        return str(self.__dict__)


class Input:

    @staticmethod
    def from_file(path: str):
        with open(path) as f:
            lines = f.readlines()

        lines = map(str.strip, lines)
        lines = map(lambda x: x.replace('NA', 'nan'), lines)
        lines = map(str.split, lines)

        nucleotides = []
        for line in lines:
            if len(line) == 1:
                nucleotides.append(Nucleotide(float(line[0])))
            elif len(line) == 2:
                nucleotides.append(Nucleotide(float(line[0]), line[1]))
            elif len(line) == 3:
                nucleotides.append(Nucleotide(float(line[0]), line[1], float(line[2])))
            else:
                raise RuntimeError(f'Invalid line: {line}')
        return Input(os.path.basename(path), nucleotides)

    def __init__(self, name: str, nucleotides: List[Nucleotide]):
        self.name = name
        self.nucleotides = nucleotides
        self.fshapes = []
        self.shapes = []
        self.bases = []
        self.profile: Dict = dict()
        self.motifs: List[Dict] = []

        for nt in nucleotides:
            self.fshapes.append(nt.fshape)
            self.shapes.append(nt.shape)
            self.bases.append(nt.base)

    def __repr__(self):
        return str(self.__dict__)

    def __len__(self):
        return len(self.fshapes)

    def compute_profile(self, query):
        if np.sum(np.isfinite(input.fshapes)) < len(query):
            return

        self.profile = mp.compute(self.fshapes,
                                  len(query),
                                  query.fshapes,
                                  n_jobs=multiprocessing.cpu_count(),
                                  preprocessing_kwargs={'window': len(query)})
        self.motifs = mp.discover.motifs(self.profile, int(len(query) / 2), 10)['motifs']

    def copy(self):
        obj = Input(self.name, self.nucleotides)
        obj.profile = self.profile
        obj.motifs = self.motifs
        return obj


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query', help='path to pattern data')
    parser.add_argument('--scramble',
                        help='shuffle each input file in a random manner to test robustness',
                        action='store_true')
    parser.add_argument('inputs', help='path to fSHAPE or SHAPE files', nargs='+')
    args = parser.parse_args()
    if not args.query:
        parser.print_help()
        sys.exit(1)
    return args


def separate_motifs(inputs: List[Input]):
    result = []
    for input in inputs:
        if len(input.motifs) == 1:
            result.append(input)
        else:
            for motif in input.motifs:
                copy = input.copy()
                copy.motifs = [motif]
                result.append(copy)
    return result


def filter_motifs_with_nans(inputs: List[Input], query):
    result = []
    for input in inputs:
        index = input.motifs[0]['motifs'][1]
        found = input.fshapes[index:index + len(query)]
        if np.isnan(found).any():
            continue
        result.append(input)
    return result


def filter_negative_motifs(inputs: List[Input], query):
    result = []
    for input in inputs:
        index = input.motifs[0]['motifs'][1]
        found = np.array(input.fshapes[index:index + len(query)])

        for i in range(len(query)):
            if query.fshapes[i] > 1.0 and found[i] <= 1.0:
                break
        else:
            result.append(input)
    return result


def euclidean(input: Input, query: Input):
    index = input.motifs[0]['motifs'][1]
    data = np.array(input.fshapes[index:index + len(query)])
    return np.linalg.norm(data - np.array(query.fshapes))


def znorm_euclidean(input: Input):
    return input.profile['mp'][input.motifs[0]['motifs'][1]]


def sequence_score(input: Input, query: Input):
    score = 0
    for i in range(len(query)):
        nt1 = query.bases[i]
        if nt1 == 'N':
            continue
        index = input.motifs[0]['motifs'][1]
        nt2 = input.bases[index + i]
        if nt1 == nt2:
            score += 2
        elif nt1 in 'CTU' and nt2 in 'CTU':
            score += 1
        elif nt1 in 'AG' and nt2 in 'AG':
            score += 1
    return score


def export_csv(inputs: List[Input], query: Input, path: str):
    objects = []
    for input in inputs:
        sample = os.path.splitext(input.name)[0]
        index = input.motifs[0]['motifs'][1]
        range_ = f'{index}-{index + len(query)}'
        sequence = ''.join(input.bases[index:index + len(query)])
        object = {
            'Sample': sample,
            'Range': range_,
            'Sequence': sequence,
            'Distance': euclidean(input, query),
            'Z-normalized': znorm_euclidean(input),
            'Sequence-Score': sequence_score(input, query)
        }
        for i in range(len(query)):
            object.update({f'fSHAPE-{i + 1}': input.fshapes[index + i]})
            object.update({f'SHAPE-{i + 1}': input.shapes[index + i]})
        objects.append(object)

    header = ['Sample', 'Range', 'Sequence', 'Z-normalized', 'Distance', 'Sequence-Score']
    for i in range(len(query)):
        header.append(f'fSHAPE-{i + 1}')
    for i in range(len(query)):
        header.append(f'SHAPE-{i + 1}')

    with open(path, 'w') as f:
        writer = csv.DictWriter(f, header)
        writer.writeheader()
        writer.writerows(objects)


def draw_everything_highlight_motifs(query, inputs):
    fig, axes = plt.subplots(len(inputs) + 1, 1, figsize=(6, 1.5 * len(inputs)))

    # draw query plot
    axes[0].plot(np.arange(len(query)), query.fshapes)
    axes[0].set_title('Query')

    # draw input plots with motifs marked on top
    for i, input in enumerate(inputs):
        axes[i + 1].set_title(input.name)

        data = input.fshapes
        length = len(data)
        axes[i + 1].plot(np.arange(length), data)
        assert len(input.motifs) == 1, "separate_motifs() function should have made input.motifs a 1-element array"

        mask = np.ones(length)
        motif = input.motifs[0]['motifs'][1]
        mask[motif:motif + len(query)] = 0

        xs = np.arange(length)
        ys = np.ma.masked_array(data, mask)
        axes[i + 1].plot(xs, ys, label=f'[{motif}:{motif + len(query)}]')

        mask = np.ones(length)
        for neighbor in input.motifs[0]['neighbors']:
            mask[neighbor:neighbor + len(query)] = 0
        ys = np.ma.masked_array(data, mask)
        axes[i + 1].plot(xs, ys, label=f'neighbours [{motif}:{motif + len(query)}]')

        axes[i + 1].legend()

    plt.tight_layout()
    plt.savefig('motifs-highlighted.pdf')


def draw_just_motifs(query, inputs):
    fig, axes = plt.subplots(len(inputs) + 1, 1, figsize=(6, 1.5 * len(inputs)))

    # draw query plot
    xs = np.arange(len(query))
    axes[0].set_title('Query')
    axes[0].plot(xs, query.fshapes)
    axes[0].set_xticks(xs)
    axes[0].set_xticklabels(query.bases)

    # draw plots with motifs
    for i, input in enumerate(inputs):
        index = input.motifs[0]['motifs'][1]
        xticklabels = []
        for j in range(len(query)):
            xticklabels.append(f'{index + j}\n{input.bases[index + j]}')

        axes[i + 1].set_title(input.name)
        axes[i + 1].plot(xs, input.fshapes[index:index + len(query)])
        axes[i + 1].set_xticks(xs)
        axes[i + 1].set_xticklabels(xticklabels)

    plt.tight_layout()
    plt.savefig('motifs-only.pdf')


if __name__ == '__main__':
    args = parse_args()

    query = Input.from_file(args.query)
    inputs = [Input.from_file(path) for path in args.inputs]

    if args.scramble:
        for input in inputs:
            random.shuffle(input.nucleotides)

    for input in inputs:
        input.compute_profile(query)

    inputs = separate_motifs(inputs)
    inputs = filter_motifs_with_nans(inputs, query)
    inputs.sort(key=znorm_euclidean)

    export_csv(inputs, query, 'output.csv')
    inputs = filter_negative_motifs(inputs, query)
    export_csv(inputs, query, 'output-filtered.csv')

    # use the best 10 profiles
    inputs = inputs[:10]

    draw_everything_highlight_motifs(query, inputs)
    draw_just_motifs(query, inputs)
