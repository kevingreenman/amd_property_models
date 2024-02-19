#!/usr/bin/python3
import sys
import csv
import os
from argparse import ArgumentParser

def add_args(parser: ArgumentParser):
    parser.add_argument("--train_data", nargs="+", type=str)
    parser.add_argument("--splits_dir", type=str)
    parser.add_argument("--num_splits", type=int)
parser = ArgumentParser()
add_args(parser)
args = parser.parse_args()
splits_dir=args.splits_dir
num_splits=args.num_splits
train_data=args.train_data

smiles=[]
preds=[]
unc=[]

for i in range(1,(num_splits+1)):
    path=os.path.join(splits_dir, f'split_{i}', 'eval_preds.csv')
    with open(path) as f:
        reader=csv.reader(f)
        next(reader)
        for line in reader:
            smiles.append(line[0])
            preds.append(line[1])
            unc.append(line[7])

target_dict={}
for train_path in train_data:
    with open(train_path) as f:
        reader=csv.reader(f)
        next(reader)
        for line in reader:
            if len(line)==2:
                target_dict[line[0]]=line[1]
            else:
                target_dict[line[0]]=line[2]

path=os.path.join(splits_dir, 'eval_preds.csv')
with open(path, 'w') as f:
    writer=csv.writer(f)
    writer.writerow(['smiles', 'Air', 'Air_ensemble_variance', "target"])
    for i in range(len(smiles)):
        writer.writerow([smiles[i], preds[i], unc[i], target_dict[smiles[i]]])


