# -*- coding: utf-8 -*-

import os
import argparse


def load_species_ids(species_file):
    speciesid_dict = {}
    with open(species_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            key, value = line.split(':', 1)
            key = key.strip()
            value = value.strip()
            if value.endswith('.pep'):
                value = value[:-4]
            speciesid_dict[key] = value
    return speciesid_dict


def main():
    parser = argparse.ArgumentParser(
        description="Generate tree_species_dict.txt from SequenceIDs.txt and SpeciesIDs.txt"
    )
    parser.add_argument(
        "--working-directory",
        required=True,
        help="Path like .../reference_pep/OrthoFinder/Results_xxx/WorkingDirectory"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output tree_species_dict.txt path"
    )
    args = parser.parse_args()

    working_directory = os.path.abspath(args.working_directory)
    output_file_path = os.path.abspath(args.output)

    input_file1_path = os.path.join(working_directory, 'SequenceIDs.txt')
    input_file2_path = os.path.join(working_directory, 'SpeciesIDs.txt')

    if not os.path.isfile(input_file1_path):
        raise FileNotFoundError(f"SequenceIDs.txt not found: {input_file1_path}")
    if not os.path.isfile(input_file2_path):
        raise FileNotFoundError(f"SpeciesIDs.txt not found: {input_file2_path}")

    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)

    speciesid_dict = load_species_ids(input_file2_path)

    tree_species_dict = {}
    with open(input_file1_path, 'r', encoding='utf-8') as file1:
        for line1 in file1:
            line1 = line1.strip()
            if not line1:
                continue
            index, transid = line1.split(':', 1)
            index = index.strip()
            transid = transid.strip()

            key = index.split('_')[0]
            if key in speciesid_dict:
                tree_species_dict[transid] = speciesid_dict[key]

    with open(output_file_path, 'w', encoding='utf-8') as result_file:
        for key, value in tree_species_dict.items():
            result_file.write(f"{key}: {value}\n")

    print(f"Done. Output written to: {output_file_path}")


if __name__ == "__main__":
    main()