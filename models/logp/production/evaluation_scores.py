#!/usr/bin/python3
from argparse import ArgumentParser
import csv
import os
import random

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen
import numpy as np

def main():
    parser = ArgumentParser()
    add_args(parser)
    args = parser.parse_args()

    smiles = []
    targets = []

    eval_smiles = []
    eval_preds = []
    eval_var = []

    for path in args.train_only_data:
        new_smiles, new_targets = load_data(path)
        smiles.extend(new_smiles)
        targets.extend(new_targets)
    if args.train_eval_data != ["None"]:
        for path in args.train_eval_data:
            new_smiles, new_targets = load_data(path)
            smiles.extend(new_smiles)
            targets.extend(new_targets)
    if args.eval_only_data != ["None"]:
        for path in args.eval_only_data:
            new_smiles, new_targets = load_data(path)
            smiles.extend(new_smiles)
            targets.extend(new_targets)

    target_dict = {}
    for i in range(len(smiles)):
        target_dict[smiles[i]] = targets[i]
    
    if args.train_eval_data != ["None"]:
        for i in range(args.num_splits):
            path = os.path.join(args.splits_dir, f"split_{i}", "eval_preds.csv")
            new_smiles, new_values, new_var = load_data(path, include_variance=True)
            eval_smiles.extend(new_smiles)
            eval_preds.extend(new_values)
            eval_var.extend(new_var)

    if args.eval_only_data != ["None"]:
        path = os.path.join(args.splits_dir, "eval", "eval_preds.csv")
        new_smiles, new_values, new_var = load_data(path, include_variance=True)
        eval_smiles.extend(new_smiles)
        eval_preds.extend(new_values)
        eval_var.extend(new_var)
    
    eval_targets = []
    for s in eval_smiles:
        eval_targets.append(target_dict[s])
    path = os.path.join(args.splits_dir, "eval_preds.csv")
    with open(path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(["smiles", "LogP", "ensemble_variance", "target"])
        for i in range(len(eval_smiles)):
            writer.writerow([eval_smiles[i], eval_preds[i], eval_var[i], eval_targets[i]])
    targets_array = np.array(eval_targets, dtype=float)
    preds_array = np.array(eval_preds, dtype=float)
    error = preds_array - targets_array
    rmse = np.sqrt(np.mean(np.square(error)))
    path = os.path.join(args.splits_dir, "eval_score.csv")
    with open(path, "w") as f:
        f.write(f"rmse, {rmse}")
    


def add_args(parser: ArgumentParser):
    parser.add_argument("--train_only_data", nargs="+", type=str)
    parser.add_argument("--train_eval_data", nargs="+", type=str)
    parser.add_argument("--eval_only_data", nargs="+", type=str)
    parser.add_argument("--splits_dir", type=str)
    parser.add_argument("--num_splits", type=int)


def load_data(path, include_variance=False):
    with open(path) as f:
        reader = csv.reader(f)
        next(reader)
        smiles = []
        targets = []
        var = []
        for line in reader:
            smiles.append(line[0])
            targets.append(line[1])
            if include_variance:
                var.append(line[2])
    if include_variance:
        return smiles, targets, var
    else:
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
        mol_r = Crippen.MolMR(mol) /100
        features.append([tpa, mol_r])

    with open(path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['topological polar surface area', 'molecular radius'])
        writer.writerows(features)

if __name__ == "__main__":
    main()