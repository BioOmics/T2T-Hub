# HMM 数据库路径 (HMM database path)
HMM="/pathway/script4animal/animalTF/all_TF.hmm"
if [[ ! -s "${HMM}.h3f" ]]; then
    echo "Preparing HMM database..."
    hmmpress "$HMM"
fi

if [[ ! -s "./genome.re.pep_output/tf_classification.txt" ]]; then
    pep="genome.re.pep"
    CPU=30
    E=1e-5

    # TF 分类 (TF classification)
    hmmsearch --cpu $CPU --tblout tmp_tbl.txt -E $E "$HMM" "$pep" > /dev/null

    awk '
        BEGIN{OFS="\t"}
        !/^#/ {
            query=$1; tf=$3; e=$5;
            if(!(query in best) || e<best_e[query]) {
                best[query]=$1"\t"tf"\tTF"
                best_e[query]=e
            }
        }
        END{for(i in best) print best[i]}
    ' tmp_tbl.txt > tf_classification.txt
    rm -f tmp_tbl.txt

    # TR 分类 (TR classification)
    human_fasta="/pathway/script4animal/animalTF/Homo_sapiens_Cof_protein.fa"
    human_annot="/pathway/script4animal/animalTF/Homo_sapiens_Cof.txt"

    cut -d " " -f1 "$pep" | sed '/^[^>]/ s/\.//g' > tmp_TR.pep

    if [[ ! -s "$(dirname "$human_fasta")/human_cof_db.dmnd" ]]; then
        diamond makedb --in "$human_fasta" -d "$(dirname "$human_fasta")/human_cof_db"
    fi

    diamond blastp \
        -q tmp_TR.pep \
        -d "$(dirname "$human_fasta")/human_cof_db.dmnd" \
        -o tmp_raw.tsv \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --evalue 1e-5 \
        --id 50 \
        --query-cover 50 \
        --max-target-seqs 1 --max-hsps 1 \
        --threads 30 --quiet

    awk -vFS="\t" -vOFS="\t" 'NR>1{print $2"\t"$4}' "$human_annot" > symbol2family.txt
    awk -vFS="\t" -vOFS="\t" 'NR==FNR{fam[$1]=$2;next} {print $1 "\t" fam[$2] "\tTR"}' symbol2family.txt tmp_raw.tsv > tr_classification.txt
    rm -f tmp_raw.tsv symbol2family.txt tmp_TR.pep

    # 合并并对第一列去重（保留第一条记录） (Merge and deduplicate the first column (keep the first record))
    mkdir -p genome.re.pep_output
    cat tf_classification.txt tr_classification.txt | \
    awk '!a[$1]++' > genome.re.pep_output/tf_classification.txt

    echo "$(pwd)"
fi
    cat tf_classification.txt tr_classification.txt | \
    awk '!a[$1]++' > genome.re.pep_output/tf_classification.txt
