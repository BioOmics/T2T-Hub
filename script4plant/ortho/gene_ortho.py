# -*- coding: utf-8 -*-

import csv
import os
import sys

csv.field_size_limit(sys.maxsize)

if len(sys.argv) != 3:
    print("Usage: python gene_ortho.py <Orthogroups.tsv> <output_dir>")
    sys.exit(1)

input_file = sys.argv[1]
output_dir = sys.argv[2]

if not os.path.isfile(input_file):
    print(f"Error: input file not found: {input_file}")
    sys.exit(1)

os.makedirs(output_dir, exist_ok=True)

column_dicts = {}

with open(input_file, 'r', encoding='utf-8') as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    header = next(reader)

    for col_name in header[1:]:
        column_dicts[col_name] = {}

    for row in reader:
        if not row:
            continue
        orthogroup = row[0]
        for col_name, value in zip(header[1:], row[1:]):
            if value:
                values = value.split(', ')
                for val in values:
                    column_dicts[col_name][val] = orthogroup

for col_name, col_dict in column_dicts.items():
    output_file = os.path.join(output_dir, f"{col_name}.gene_ortho.txt")
    with open(output_file, 'w', encoding='utf-8') as outfile:
        for key, val in col_dict.items():
            outfile.write(f"{key}\t{val}\n")

print(f"Done. Results written to: {output_dir}")