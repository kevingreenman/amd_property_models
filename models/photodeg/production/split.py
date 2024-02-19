# !/usr/bin/python3
from argparse import ArgumentParser
import csv
import os
import random
from collections import defaultdict
import numpy as np


def main():
    parser = ArgumentParser()
    add_args(parser)
    args = parser.parse_args()

    loaded_smiles = []
    loaded_targets = defaultdict(list)
    for path in args.train_data:
        new_smiles, new_targets = load_data(path, args.data_dir)
        loaded_smiles.extend(new_smiles)
        for i in range(len(new_smiles)):
            loaded_targets[new_smiles[i]].append(new_targets[i])
    # list_len=len(loaded_smiles)
    crossval_smiles = set(loaded_smiles)
    crossval_smiles = list(crossval_smiles)
    # set_len=len(crossval_smiles)
    # assert(list_len==set_len)
    crossval_targets = [np.mean(loaded_targets[s], axis=0) for s in crossval_smiles]

    os.makedirs(args.splits_dir, exist_ok=True)

    temp = list(zip(crossval_smiles, crossval_targets))
    random.seed(0)
    random.shuffle(temp)
    crossval_smiles, crossval_targets = zip(*temp)
    crossval_smiles, crossval_targets = list(crossval_smiles), list(crossval_targets)

    chunk_size = len(crossval_smiles) // args.num_splits
    remainder = len(crossval_smiles) % args.num_splits

    splits_smiles = []
    splits_targets = []
    for i in range(args.num_splits):
        splits_smiles.append(crossval_smiles[i*chunk_size:(i+1)*chunk_size])
        if i < remainder:
            splits_smiles[i].append(crossval_smiles[-i-1])
        splits_targets.append(crossval_targets[i*chunk_size:(i+1)*chunk_size])
        if i < remainder:
            splits_targets[i].append(crossval_targets[-i-1])

    for i in range(args.num_splits):
        split_path = os.path.join(args.splits_dir, f"split_{i+1}")
        os.makedirs(split_path, exist_ok=True)
        split_train_smiles = []
        for j in range(args.num_splits):
            if j != i:
                split_train_smiles.extend(splits_smiles[j])
        split_train_targets = []
        for j in range(args.num_splits):
            if j != i:
                split_train_targets.extend(splits_targets[j])
        write_data(
            smiles=split_train_smiles,
            targets=split_train_targets,
            path=os.path.join(split_path, "trainval.csv")
        )
        write_data(
            smiles=splits_smiles[i],
            targets=splits_targets[i],
            path=os.path.join(split_path, "eval.csv")
        )


def add_args(parser: ArgumentParser):
    parser.add_argument("--train_data", nargs="+", type=str)
    parser.add_argument("--splits_dir", type=str)
    parser.add_argument("--num_splits", type=int)
    parser.add_argument("--data_dir", type=str, default="/home/gridsan/kgreenman/amd_property_models/datasets/photodeg/")


def load_data(path, data_dir):
    with open(data_dir + path) as f:
        reader = csv.reader(f)
        next(reader)
        smiles = []
        targets = []
        for line in reader:
            targs = []
            if len(line) == 2:
                targs.extend([line[1]] + [''] * 5)
            else:
                targs.extend([
                    line[2], line[1], line[3], line[4], line[5], line[6]
                ])
            targs = [float(t) if t != '' else np.nan for t in targs]
            smiles.append(line[0])
            targets.append(targs)
    return smiles, targets


def write_data(smiles, targets, path):
    with open(path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles", "Air", "N2", "IPA", "NaN3", "Imidazole", "Bzq"])
        for i in range(len(smiles)):
            row = []
            row.append(smiles[i])
            row.extend(targets[i])
            writer.writerow(row)


if __name__ == "__main__":
    main()
