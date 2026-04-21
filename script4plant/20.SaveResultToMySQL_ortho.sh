# Usage: bash 20.SaveResultToMySQL.sh namaList

namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read i;do
python /pathway/script4plant/20.getAnotateSeq_ortho.py ${i}
python /pathway/script4plant/20.getGeneFunction.py ${i}
python /pathway/script4plant/20.getGeneTable.py ${i}
bash /pathway/script4plant/20.spe_gene_sql_ortho.sh ${i}
done