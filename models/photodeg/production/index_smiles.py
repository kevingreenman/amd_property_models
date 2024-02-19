#!/usr/bin/python3

import sys
import csv

input_path = sys.argv[1]
output_path = sys.argv[2]

smiles=[]

with open(input_path,'r') as f:
    reader=csv.reader(f)
    next(reader)
    for line in reader:
        smiles.append(line[0])

with open(output_path,'w') as f:
    writer=csv.writer(f)
    writer.writerow(['index', 'smiles'])
    for i, s in enumerate(smiles):
        writer.writerow([i, s])
