#!/usr/bin/python3
from argparse import ArgumentParser
import csv
import os
import random

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen


def main():
    parser = ArgumentParser()
    add_args(parser)
    args = parser.parse_args()

    train_smiles = []
    train_targets = []
    crossval_smiles = []
    crossval_targets = []
    test_smiles = []
    test_targets = []
    for path in args.train_only_data:
        new_smiles, new_targets = load_data(path)
        train_smiles.extend(new_smiles)
        train_targets.extend(new_targets)
    if args.train_eval_data != ["None"]:
        for path in args.train_eval_data:
            new_smiles, new_targets = load_data(path)
            crossval_smiles.extend(new_smiles)
            crossval_targets.extend(new_targets)
    if args.eval_only_data != ["None"]:
        for path in args.eval_only_data:
            new_smiles, new_targets = load_data(path)
            test_smiles.extend(new_smiles)
            test_targets.extend(new_targets)

    os.makedirs(args.splits_dir, exist_ok=True)

    if args.train_eval_data != ["None"]:
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
            splits_smiles.append(crossval_smiles[i * chunk_size:(i + 1) * chunk_size])
            if i < remainder:
                splits_smiles[i].append(crossval_smiles[-i - 1])
            splits_targets.append(crossval_targets[i * chunk_size:(i + 1) * chunk_size])
            if i < remainder:
                splits_targets[i].append(crossval_targets[-i - 1])

        for i in range(args.num_splits):
            split_path = os.path.join(args.splits_dir, f"split_{i}")
            os.makedirs(split_path, exist_ok=True)
            split_train_smiles = []
            split_train_smiles.extend(train_smiles)
            for j in range(args.num_splits):
                if j != i:
                    split_train_smiles.extend(splits_smiles[j])
            split_train_targets = []
            split_train_targets.extend(train_targets)
            for j in range(args.num_splits):
                if j != i:
                    split_train_targets.extend(splits_targets[j])
            write_data(
                smiles=split_train_smiles,
                targets=split_train_targets,
                path=os.path.join(split_path, "trainval.csv")
            )
            write_features(
                smiles=split_train_smiles,
                path=os.path.join(split_path, "trainval_features.csv")
            )
            write_data(
                smiles=splits_smiles[i],
                targets=splits_targets[i],
                path=os.path.join(split_path, "eval.csv")
            )
            write_features(
                smiles=splits_smiles[i],
                path=os.path.join(split_path, "eval_features.csv")
            )

    eval_path = os.path.join(args.splits_dir, "eval")
    os.makedirs(eval_path, exist_ok=True)
    write_data(
        smiles=train_smiles + crossval_smiles,
        targets=train_targets + crossval_targets,
        path=os.path.join(eval_path, "trainval.csv")
    )
    write_features(
        smiles=train_smiles + crossval_smiles,
        path=os.path.join(eval_path, "trainval_features.csv")
    )
    write_data(
        smiles=test_smiles,
        targets=test_targets,
        path=os.path.join(eval_path, "eval.csv")
    )
    write_features(
        smiles=test_smiles,
        path=os.path.join(eval_path, "eval_features.csv")
    )


def add_args(parser: ArgumentParser):
    parser.add_argument("--train_only_data", nargs="+", type=str)
    parser.add_argument("--train_eval_data", nargs="+", type=str)
    parser.add_argument("--eval_only_data", nargs="+", type=str)
    parser.add_argument("--splits_dir", type=str)
    parser.add_argument("--num_splits", type=int)
    parser.add_argument("--data_dir", type=str, default="/home/gridsan/kgreenman/amd_property_models/datasets/logp/")


def load_data(path, data_dir):
    with open(data_dir + path) as f:
        reader = csv.reader(f)
        next(reader)
        smiles = []
        targets = []
        for line in reader:
            smiles.append(line[0])
            targets.append(line[1])
    return smiles, targets


def write_data(smiles, targets, path):
    with open(path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles", "LogP"])
        for i in range(len(smiles)):
            writer.writerow([smiles[i], targets[i]])


def write_features(smiles, path):
    features = []
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        tpa = AllChem.CalcTPSA(mol) / 100
        mol_r = Crippen.MolMR(mol) / 100
        features.append([tpa, mol_r])

    with open(path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['topological polar surface area', 'molecular radius'])
        writer.writerows(features)


if __name__ == "__main__":
    main()
