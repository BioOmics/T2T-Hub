# -*- coding: utf-8 -*-

import os
import gzip
import argparse
from collections import defaultdict
from tqdm import tqdm


def get_reference_pep_dir(working_directory: str) -> str:
    """
    working_directory:
      .../reference_pep/OrthoFinder/Results_xxx/WorkingDirectory
    return:
      .../reference_pep
    """
    return os.path.dirname(os.path.dirname(os.path.dirname(working_directory)))


def load_sequence_ids(sequence_file: str) -> dict:
    """
    read SequenceIDs.txt
    format:
      0_0: gene_name ...
    """
    mapping = {}
    with open(sequence_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            key, value = line.split(":", 1)
            mapping[key.strip()] = value.strip()
    return mapping


def load_species_ids(species_file: str) -> dict:
    """
    read SpeciesIDs.txt
    format:
      0: 07369.pep
      1: Arabidopsis_thaliana_Col-CEN.pep
    Processing:
      Remove the trailing .pep
    Returns:
      {"0": "07369", "1": "Arabidopsis_thaliana_Col-CEN", ...}
    """
    mapping = {}
    with open(species_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            key, value = line.split(":", 1)
            key = key.strip()
            value = value.strip()
            if value.endswith(".pep"):
                value = value[:-4]
            mapping[key] = value
    return mapping


def load_cluster_mappings(cluster_file: str):
    """
    read clusters_OrthoFinder_I1.2.txt_id_pairs.txt

    Preprocessing required:
    - Remove all content before 'begin'
    - Remove the last standalone ')'

    Finally, keep only lines like:
      0 1_9259 1_9262 ... $
      1 ...
    """
    with open(cluster_file, "r", encoding="utf-8") as f:
        content = f.read()

    if "begin" not in content:
        raise ValueError(f"'begin' not found in cluster file: {cluster_file}")

    content = content.split("begin", 1)[1].strip()

    if content.endswith(")"):
        content = content[:-1].strip()

    clusterid_mapping = {}
    gene_orthogroup_mapping = defaultdict(list)

    for block in content.split("$"):
        parts = block.strip().split()
        if not parts:
            continue

        og_id = parts[0]
        genes = set(parts[1:])

        clusterid_mapping[og_id] = genes
        for gene in genes:
            gene_orthogroup_mapping[gene].append(og_id)

    return clusterid_mapping, gene_orthogroup_mapping


def infer_ortho_type(first_column, second_column, matching_orthogroups,
                     clusterid_mapping, gene_orthogroup_mapping):
    """
    Determine the ortho type based on the original logic.
    Returns:
      ortho_type, selected_orthogroups
    """
    species_prefix = second_column.split("_")[0] + "_"

    if len(matching_orthogroups) == 1:
        og = matching_orthogroups[0]

        if og not in gene_orthogroup_mapping.get(second_column, []):
            return None, []

        num_elements = 0
        ortho = None

        for value in clusterid_mapping[og]:
            if value.startswith(species_prefix):
                num_elements += 1
                if num_elements > 1:
                    ortho = "1-M"
                    break

        if num_elements == 1:
            ortho = "1-1"

        if ortho is None:
            return None, []

        return ortho, [og]

    counts = 0
    num_elements = 0

    for og in matching_orthogroups:
        if og in gene_orthogroup_mapping.get(second_column, []):
            counts += 1
            for value in clusterid_mapping[og]:
                if counts == 1 and value.startswith(species_prefix):
                    num_elements += 1

    if num_elements == 1 and counts == 1:
        ortho = "1-1"
    elif num_elements == 1 and counts > 1:
        ortho = "M-1"
    elif num_elements > 1 and counts == 1:
        ortho = "1-M"
    elif num_elements > 1 and counts > 1:
        ortho = "M-M"
    else:
        return None, []

    if ortho in ("M-M", "M-1"):
        return ortho, matching_orthogroups

    selected = []
    for og in matching_orthogroups:
        if second_column in clusterid_mapping[og]:
            selected.append(og)

    return ortho, selected


def generate_combine_records(aa, working_directory, clusterid_mapping,
                             gene_orthogroup_mapping, sequence_id_mapping):
    """
    Scan Blast{aa}_*.txt.gz to generate combine records.
    Each record:
      [search_gene, trans_gene, ortho, orthogroup, similarity, score]
    """
    records = []

    for filename in tqdm(sorted(os.listdir(working_directory)), desc=f"Blast{aa}"):
        if not filename.startswith(f"Blast{aa}_"):
            continue
        if not filename.endswith(".gz"):
            continue

        file_path = os.path.join(working_directory, filename)

        with gzip.open(file_path, "rt", encoding="utf-8") as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 3:
                    continue

                first_column = fields[0]
                second_column = fields[1]

                search_gene = sequence_id_mapping.get(first_column)
                trans_gene = sequence_id_mapping.get(second_column)

                if search_gene is None or trans_gene is None:
                    continue

                matching_orthogroups = gene_orthogroup_mapping.get(first_column, [])
                if not matching_orthogroups:
                    continue

                ortho, selected_ogs = infer_ortho_type(
                    first_column,
                    second_column,
                    matching_orthogroups,
                    clusterid_mapping,
                    gene_orthogroup_mapping
                )

                if ortho is None:
                    continue

                similarity = fields[2]
                score = fields[-1]

                for og in selected_ogs:
                    orthogroup = f"OG{int(og):07d}"
                    records.append([
                        search_gene,
                        trans_gene,
                        ortho,
                        orthogroup,
                        similarity,
                        score
                    ])

    return records


def write_combine_file(records, output_file):
    with open(output_file, "w", encoding="utf-8") as f:
        for row in records:
            f.write("\t".join(row) + "\n")


def aggregate_final(records):
    """
    Aggregate by orthogroup + search_gene
    Output format:
      search_gene  trans_genes  ortho_types  similarities  scores  orthogroup
    """
    merged = defaultdict(lambda: defaultdict(lambda: {
        "trans_genes": [],
        "ortho": [],
        "similarity": [],
        "score": []
    }))

    for search_gene, trans_gene, ortho, orthogroup, similarity, score in records:
        merged[orthogroup][search_gene]["trans_genes"].append(trans_gene)
        merged[orthogroup][search_gene]["ortho"].append(ortho)
        merged[orthogroup][search_gene]["similarity"].append(similarity)
        merged[orthogroup][search_gene]["score"].append(score)

    final_lines = []
    for orthogroup, search_dict in merged.items():
        for search_gene, data in search_dict.items():
            final_lines.append([
                search_gene,
                ",".join(data["trans_genes"]),
                ",".join(data["ortho"]),
                ",".join(data["similarity"]),
                ",".join(data["score"]),
                orthogroup
            ])
    return final_lines


def write_final_file(lines, output_file):
    with open(output_file, "w", encoding="utf-8") as f:
        for row in lines:
            f.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate combine and final gene blast files from OrthoFinder WorkingDirectory."
    )
    parser.add_argument(
        "--working-directory",
        required=True,
        help="Path like .../reference_pep/OrthoFinder/Results_xxx/WorkingDirectory"
    )
    parser.add_argument(
        "--aa",
        required=True,
        help="Species numeric ID, e.g. 0"
    )
    args = parser.parse_args()

    working_directory = os.path.abspath(args.working_directory)
    aa = str(args.aa)

    if not os.path.isdir(working_directory):
        raise FileNotFoundError(f"working-directory not found: {working_directory}")

    reference_pep_dir = get_reference_pep_dir(working_directory)
    combine_dir = os.path.join(reference_pep_dir, "combine")
    final_dir = os.path.join(reference_pep_dir, "final")

    os.makedirs(combine_dir, exist_ok=True)
    os.makedirs(final_dir, exist_ok=True)

    cluster_file = os.path.join(working_directory, "clusters_OrthoFinder_I1.2.txt_id_pairs.txt")
    sequence_file = os.path.join(working_directory, "SequenceIDs.txt")
    species_file = os.path.join(working_directory, "SpeciesIDs.txt")

    if not os.path.isfile(cluster_file):
        raise FileNotFoundError(f"cluster file not found: {cluster_file}")
    if not os.path.isfile(sequence_file):
        raise FileNotFoundError(f"SequenceIDs.txt not found: {sequence_file}")
    if not os.path.isfile(species_file):
        raise FileNotFoundError(f"SpeciesIDs.txt not found: {species_file}")

    species_mapping = load_species_ids(species_file)
    if aa not in species_mapping:
        raise KeyError(f"aa={aa} not found in SpeciesIDs.txt")

    species_name = species_mapping[aa]

    combine_file = os.path.join(combine_dir, f"combine{aa}.txt")
    final_file = os.path.join(final_dir, f"{species_name}_geneblast.txt")

    clusterid_mapping, gene_orthogroup_mapping = load_cluster_mappings(cluster_file)
    sequence_id_mapping = load_sequence_ids(sequence_file)

    records = generate_combine_records(
        aa=aa,
        working_directory=working_directory,
        clusterid_mapping=clusterid_mapping,
        gene_orthogroup_mapping=gene_orthogroup_mapping,
        sequence_id_mapping=sequence_id_mapping
    )

    write_combine_file(records, combine_file)
    print(f"Combine file written to: {combine_file}")

    final_lines = aggregate_final(records)
    write_final_file(final_lines, final_file)
    print(f"Final file written to: {final_file}")

    print("Done.")


if __name__ == "__main__":
    main()