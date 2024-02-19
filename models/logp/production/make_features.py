#!/usr/bin/python3
from argparse import ArgumentParser
import csv
import os
import random
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen

smiles_path = sys.argv[1]
features_path = sys.argv[2]

def main(smiles_path, features_path):
    smiles = load_data(smiles_path)
    features_dir = os.path.dirname(features_path)
    os.makedirs(features_dir, exist_ok=True)
    write_features(smiles=smiles, path=features_path)


def load_data(path):
    with open(path) as f:
        reader = csv.reader(f)
        next(reader)
        smiles = []
        targets = []
        for line in reader:
            smiles.append(line[0])
    return smiles


def write_features(smiles, path):
    features = []
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        tpa = AllChem.CalcTPSA(mol) / 100
        mol_r = Crippen.MolMR(mol) /100
        features.append([tpa, mol_r])

    with open(path, 'w', newline = '') as f:
        writer = csv.writer(f)
        writer.writerow(['topological polar surface area', 'molecular radius'])
        writer.writerows(features)


if __name__ == "__main__":
    main(smiles_path=smiles_path, features_path=features_path)
