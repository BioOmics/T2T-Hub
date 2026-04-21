#!/bin/bash

namaList="$1"

if [ ! -f "$namaList" ] || [ ! -s "$namaList" ]; then
    echo "Error: namelist file not found or empty: $namaList"
    exit 1
fi

# 取第一列 ID (get the first column ID)
ID=$(awk 'NR==1{print $1}' "$namaList")
if [ -z "${ID:-}" ]; then
    echo "Error: failed to parse ID from $namaList"
    exit 1
fi

cd "/pathway/UserUpload/${ID}"

WORKDIR=$(pwd)
REFDIR="${WORKDIR}/reference_pep"
USERPEP="${WORKDIR}/genome.re.pep"

echo "WORKDIR: $WORKDIR"
echo "REFDIR: $REFDIR"
echo "USERPEP: $USERPEP"

# 检查 reference_pep 目录是否存在且包含参考蛋白文件 (check if reference_pep directory exists and contains reference protein files)
if [ ! -d "$REFDIR" ]; then
    echo "Error: reference_pep directory not found: $REFDIR"
    exit 1
fi

if ! find "$REFDIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.faa" -o -name "*.pep" -o -name "*.pep.fa" \) | grep -q .; then
    echo "Error: reference_pep directory exists but contains no protein files: $REFDIR"
    exit 1
fi

# 检查用户自己的蛋白文件 (check if user protein file exists)
if [ ! -f "$USERPEP" ]; then
    echo "Error: user protein file not found: $USERPEP"
    exit 1
fi

# 复制当前用户蛋白到 reference_pep 目录 (copy current user protein to reference_pep directory)
BASENAME=$(basename "$WORKDIR")
TARGETPEP="${REFDIR}/${BASENAME}.pep"
cp -f "$USERPEP" "$TARGETPEP"
echo "Copied user pep to: $TARGETPEP"


# 对目录下所有蛋白文件去掉序列行中的 '.' ( remove '.' from sequence lines in all protein files in the directory)
find "$REFDIR" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.faa" -o -name "*.pep" -o -name "*.pep.fa" \) -print0 | \
while IFS= read -r -d '' f; do
    echo "Cleaning dots in: $f"
    sed -i '/^>/! s/\.//g' "$f"

    tmp="${f}.tmp"
    awk '/^>/ {print $1; next} {print}' "$f" > "$tmp" && mv "$tmp" "$f"
done


OLD_ENV="${CONDA_DEFAULT_ENV:-}"

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
else
    echo "Error: conda not found."
    exit 1
fi

if ! command -v orthofinder >/dev/null 2>&1; then
    echo "Error: orthofinder not found in orthofinderV3 environment."
    exit 1
fi

echo "Running OrthoFinder in env: ${CONDA_DEFAULT_ENV:-unknown}"
cd "$WORKDIR"
orthofinder -f reference_pep -t 20 -a 10

echo "OrthoFinder finished."

# 查找最新的 Results_* 目录 (find the latest Results_* directory)
ORTHO_BASE="${WORKDIR}/reference_pep/OrthoFinder"
LATEST_RESULT=$(find "$ORTHO_BASE" -maxdepth 1 -type d -name "Results_*" | sort | tail -n 1)

if [ -z "${LATEST_RESULT:-}" ]; then
    echo "Error: no OrthoFinder result directory found under $ORTHO_BASE"
    exit 1
fi

echo "LATEST_RESULT: $LATEST_RESULT"

ORTHO_TSV="${LATEST_RESULT}/Orthogroups/Orthogroups.tsv"
GENE_ORTHO_OUT="${WORKDIR}/reference_pep/gene_ortho"
GENE_ORTHO_PY="/pathway/script4plant/ortho/gene_ortho.py"

ORTHO_SPECIES_TREE_PY="/pathway/script4plant/ortho/ortho_species_tree.py"
TREE_SPECIES_PY="/pathway/script4plant/ortho/tree_species.py"

COMBINE_FINAL_PY="/pathway/script4plant/ortho/gene_blast_combine_final.py"
WORKING_DIR="${LATEST_RESULT}/WorkingDirectory"
TREE_SPECIES_DICT="${WORKDIR}/reference_pep/tree_species_dict.txt"
SPECIES_IDS="${WORKING_DIR}/SpeciesIDs.txt"

echo "ORTHO_TSV: $ORTHO_TSV"
echo "GENE_ORTHO_OUT: $GENE_ORTHO_OUT"
echo "GENE_ORTHO_PY: $GENE_ORTHO_PY"
echo "ORTHO_SPECIES_TREE_PY: $ORTHO_SPECIES_TREE_PY"
echo "TREE_SPECIES_PY: $TREE_SPECIES_PY"
echo "TREE_SPECIES_DICT: $TREE_SPECIES_DICT"
echo "COMBINE_FINAL_PY: $COMBINE_FINAL_PY"
echo "WORKING_DIR: $WORKING_DIR"
echo "SPECIES_IDS: $SPECIES_IDS"

if [ ! -f "$ORTHO_TSV" ]; then
    echo "Error: Orthogroups.tsv not found: $ORTHO_TSV"
    exit 1
fi

if [ ! -f "$GENE_ORTHO_PY" ]; then
    echo "Error: gene_ortho.py not found: $GENE_ORTHO_PY"
    exit 1
fi

if [ ! -f "$ORTHO_SPECIES_TREE_PY" ]; then
    echo "Error: ortho_species_tree.py not found: $ORTHO_SPECIES_TREE_PY"
    exit 1
fi

if [ ! -f "$TREE_SPECIES_PY" ]; then
    echo "Error: tree_species.py not found: $TREE_SPECIES_PY"
    exit 1
fi

if [ ! -f "$COMBINE_FINAL_PY" ]; then
    echo "Error: combine_final.py not found: $COMBINE_FINAL_PY"
    exit 1
fi

if [ ! -f "$SPECIES_IDS" ]; then
    echo "Error: SpeciesIDs.txt not found: $SPECIES_IDS"
    exit 1
fi

mkdir -p "$GENE_ORTHO_OUT"

python "$GENE_ORTHO_PY" "$ORTHO_TSV" "$GENE_ORTHO_OUT"
echo "gene_ortho.py finished."

python "$TREE_SPECIES_PY" \
    --working-directory "$WORKING_DIR" \
    --output "$TREE_SPECIES_DICT"
echo "tree_species.py finished."

python "$ORTHO_SPECIES_TREE_PY" \
    --user-id "$ID" \
    --result-dir "$LATEST_RESULT" \
    --tree-species-dict "$TREE_SPECIES_DICT"
echo "ortho_species_tree.py finished."

AA=$(awk -F': ' -v id="${ID}.pep" '$2 == id {print $1}' "$SPECIES_IDS" | head -n 1)

if [ -z "${AA:-}" ]; then
    echo "Error: failed to find aa for ${ID}.pep in $SPECIES_IDS"
    exit 1
fi

echo "Detected aa for current user: $AA"

python "$COMBINE_FINAL_PY" \
    --working-directory "$WORKING_DIR" \
    --aa "$AA"
echo "gene_blast_combine_final.py finished."

# 切回原环境
if [ -n "${OLD_ENV:-}" ]; then
    conda activate "$OLD_ENV"
else
    conda deactivate || true
fi

echo "Back to env: ${CONDA_DEFAULT_ENV:-base}"
echo "All done."