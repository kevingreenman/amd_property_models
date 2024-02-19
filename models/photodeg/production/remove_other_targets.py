#!/usr/bin/python3

import sys
import csv
import numpy as np
import os

input_path=sys.argv[1]
output_path=sys.argv[2]

lines=[]
with open(input_path, 'r') as f:
    reader=csv.reader(f)
    next(reader)
    for line in reader:
        lines.append(line)

os.makedirs(os.path.dirname(output_path), exist_ok=True)

with open(output_path, 'w') as f:
    writer=csv.writer(f)
    writer.writerow(['smiles', 'log10_Air_rate', 'ensemble_variance'])
    for line in lines:
        writer.writerow([line[0], line[1], line[7]])