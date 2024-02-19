#!/usr/bin/python3
from argparse import ArgumentParser
import csv
import os
import random


def main():
    parser = ArgumentParser()
    add_args(parser)
    args = parser.parse_args()

    train_smiles = []
    train_solvents = []
    train_targets = []
    crossval_smiles = []
    crossval_solvents = []
    crossval_targets = []
    test_smiles = []
    test_solvents = []
    test_targets = []
    for path in args.train_only_data:
        new_smiles, new_solvents, new_targets = load_data(path)
        train_smiles.extend(new_smiles)
        train_solvents.extend(new_solvents)
        train_targets.extend(new_targets)
    if args.train_eval_data != ["None"]:
        for path in args.train_eval_data:
            new_smiles, new_solvents, new_targets = load_data(path)
            crossval_smiles.extend(new_smiles)
            crossval_solvents.extend(new_solvents)
            crossval_targets.extend(new_targets)
    if args.eval_only_data != ["None"]:
        for path in args.eval_only_data:
            new_smiles, new_solvents, new_targets = load_data(path)
            test_smiles.extend(new_smiles)
            test_solvents.extend(new_solvents)
            test_targets.extend(new_targets)

    os.makedirs(args.splits_dir, exist_ok=True)

    if args.train_eval_data != ["None"]:
        temp = list(zip(crossval_smiles, crossval_solvents, crossval_targets))
        random.seed(0)
        random.shuffle(temp)
        crossval_smiles, crossval_solvents, crossval_targets = zip(*temp)
        crossval_smiles, crossval_solvents, crossval_targets = list(crossval_smiles), list(crossval_solvents), list(crossval_targets)

        chunk_size = len(crossval_smiles) // args.num_splits
        remainder = len(crossval_smiles) % args.num_splits

        splits_smiles = []
        splits_solvents = []
        splits_targets = []
        for i in range(args.num_splits):
            splits_smiles.append(crossval_smiles[i*chunk_size:(i+1)*chunk_size])                
            splits_solvents.append(crossval_solvents[i*chunk_size:(i+1)*chunk_size])
            splits_targets.append(crossval_targets[i*chunk_size:(i+1)*chunk_size])
            if i < remainder:
                splits_smiles[i].append(crossval_smiles[-i-1])
                splits_solvents[i].append(crossval_solvents[-i-1])
                splits_targets[i].append(crossval_targets[-i-1])

        for i in range(args.num_splits):
            split_path = os.path.join(args.splits_dir, f"split_{i}")
            os.makedirs(split_path, exist_ok=True)
            split_train_smiles = []
            split_train_solvents = []
            split_train_targets = []
            split_train_smiles.extend(train_smiles)
            split_train_solvents.extend(train_solvents)
            split_train_targets.extend(train_targets)
            for j in range(args.num_splits):
                if j != i:
                    split_train_smiles.extend(splits_smiles[j])
                    split_train_solvents.extend(splits_solvents[j])
                    split_train_targets.extend(splits_targets[j])
            write_data(
                smiles=split_train_smiles,
                solvents=split_train_solvents,
                targets=split_train_targets,
                path=os.path.join(split_path, "trainval.csv")
            )
            write_data(
                smiles=splits_smiles[i],
                solvents=splits_solvents[i],
                targets=splits_targets[i],
                path=os.path.join(split_path, "eval.csv")
            )

    eval_path = os.path.join(args.splits_dir, "eval")
    os.makedirs(eval_path, exist_ok=True)
    write_data(
        smiles=train_smiles + crossval_smiles,
        solvents=train_solvents + crossval_solvents,
        targets=train_targets + crossval_targets,
        path=os.path.join(eval_path, "trainval.csv")
    )
    write_data(
        smiles=test_smiles,
        solvents=test_solvents,
        targets=test_targets,
        path=os.path.join(eval_path, "eval.csv")
    )


def add_args(parser: ArgumentParser):
    parser.add_argument("--train_only_data", nargs="+", type=str)
    parser.add_argument("--train_eval_data", nargs="+", type=str)
    parser.add_argument("--eval_only_data", nargs="+", type=str)
    parser.add_argument("--splits_dir", type=str)
    parser.add_argument("--num_splits", type=int)


def load_data(path):
    with open(path) as f:
        reader = csv.reader(f)
        next(reader)
        smiles = []
        solvents = []
        targets = []
        for line in reader:
            smiles.append(line[0])
            solvents.append(line[1])
            targets.append(line[2])
    return smiles, solvents, targets


def write_data(smiles, solvents, targets, path):
    with open(path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles", "solvent", "peakwavs_max"])
        for i in range(len(smiles)):
            writer.writerow([smiles[i], solvents[i], targets[i]])


if __name__ == "__main__":
    main()
