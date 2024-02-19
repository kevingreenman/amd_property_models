#!/usr/bin/python3

import csv
import sys

input_path = sys.argv[1]
output_path = sys.argv[2]

data = []

with open(input_path, "r") as f:
    reader = csv.reader(f)
    for line in reader:
        data.append(line)

with open(output_path, "w") as f:
    writer = csv.writer(f)
    for line in data:
        writer.writerow(line[1:])