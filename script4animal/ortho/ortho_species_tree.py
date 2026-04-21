# -*- coding: utf-8 -*-

import os
import re
import argparse
from tqdm import tqdm
from collections import OrderedDict


def is_float(value: str) -> bool:
    if "_" in value:
        return False
    try:
        float(value)
        return True
    except ValueError:
        return False


def load_key_value_file(filepath: str, sep: str, maxsplit: int = 1) -> dict:
    data = {}
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            key, value = line.split(sep, maxsplit)
            data[key.strip()] = value.strip()
    return data


def treeorder(tree: str, sequence_id_mapping: dict) -> list:
    regex = re.compile(r'([a-zA-Z0-9_]+[^,:)(]+)(?=:)')
    tree_order = []
    matches = regex.findall(tree)
    for item in matches:
        if is_float(item):
            continue
        if item in sequence_id_mapping:
            tree_order.append(sequence_id_mapping[item])
    return tree_order


def convert_tree(tree: str, sequence_id_mapping: dict) -> str:
    regex = re.compile(r'([a-zA-Z0-9_]+[^,:)(]+)(?=:)')
    matches = regex.findall(tree)
    for old in matches:
        if is_float(old):
            continue
        if old in sequence_id_mapping:
            new = sequence_id_mapping[old]
            tree = tree.replace(old, new, 1)
    return tree


def main():
    parser = argparse.ArgumentParser(
        description="Generate ortho_species_tree.txt from OrthoFinder latest Results_* directory."
    )
    parser.add_argument(
        "--user-id",
        required=True,
        help="UserUpload ID, e.g. 83356"
    )
    parser.add_argument(
        "--result-dir",
        required=True,
        help="Latest OrthoFinder Results_* directory, e.g. .../reference_pep/OrthoFinder/Results_Apr06"
    )
    parser.add_argument(
        "--tree-species-dict",
        required=True,
        help="tree_species_dict.txt path"
    )
    args = parser.parse_args()

    user_id = args.user_id
    result_dir = args.result_dir
    tree_species_dict_path = args.tree_species_dict

    working_directory = os.path.join(result_dir, "WorkingDirectory")
    gene_tree_directory = os.path.join(working_directory, "Trees_ids")
    sequence_id_file_path = os.path.join(working_directory, "SequenceIDs.txt")

    output_file_path = (
        f"/public/workspace/biobigdata/project/t2thub/UserUpload/"
        f"{user_id}/reference_pep/ortho_species_tree.txt"
    )

    if not os.path.isdir(result_dir):
        raise FileNotFoundError(f"Results directory not found: {result_dir}")
    if not os.path.isdir(working_directory):
        raise FileNotFoundError(f"WorkingDirectory not found: {working_directory}")
    if not os.path.isdir(gene_tree_directory):
        raise FileNotFoundError(f"Trees_ids directory not found: {gene_tree_directory}")
    if not os.path.isfile(sequence_id_file_path):
        raise FileNotFoundError(f"SequenceIDs.txt not found: {sequence_id_file_path}")
    if not os.path.isfile(tree_species_dict_path):
        raise FileNotFoundError(f"tree_species_dict.txt not found: {tree_species_dict_path}")

    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)

    sequence_id_mapping = load_key_value_file(sequence_id_file_path, ":", 1)
    gene_species_dict = load_key_value_file(tree_species_dict_path, ": ", 1)

    merged_lines = OrderedDict()

    tree_files = sorted(
        [
            f for f in os.listdir(gene_tree_directory)
            if f.endswith(".txt")
        ]
    )

    print(f"Found {len(tree_files)} tree files in {gene_tree_directory}")

    for tree_filename in tqdm(tree_files, desc="Processing orthogroups"):
        tree_file_path = os.path.join(gene_tree_directory, tree_filename)
        orthogroup = tree_filename.replace(".txt", "")

        with open(tree_file_path, "r", encoding="utf-8") as tree_file:
            tree = tree_file.read().strip()

        if not tree:
            continue

        tree_new = convert_tree(tree, sequence_id_mapping)
        tree_order = treeorder(tree, sequence_id_mapping)

        species = []
        for item in tree_order:
            if item in gene_species_dict:
                species.append(gene_species_dict[item])

        merged_lines[orthogroup] = {
            "species": species,
            "tree": tree_new
        }

    with open(output_file_path, "w", encoding="utf-8") as output_file:
        for orthogroup, data in merged_lines.items():
            species_str = ";".join(data["species"])
            output_file.write(f"{orthogroup}\t{species_str}\t{data['tree']}\n")

    print(f"Done. Output written to: {output_file_path}")


if __name__ == "__main__":
    main()