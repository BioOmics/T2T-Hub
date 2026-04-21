# Usage: bash 14.getKEGG.sh namaList

namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read i;do
species=$i
echo ${species}

cp /pathway/script4animal/11.ko00001.filter4Animal.txt ./
cp /pathway/script4animal/14.GOterm.txt ./

awk -vFS="\t" -vOFS="\t" '{print $5, $3}' 11.ko00001.filter4Animal.txt > level3.txt
awk -vFS="\t" -vOFS="\t" '{print $5, $2}' 11.ko00001.filter4Animal.txt > level2.txt
awk -vFS="\t" -vOFS="\t" '{print $5, $1}' 11.ko00001.filter4Animal.txt > level1.txt

awk -vFS="\t" -vOFS="\t" '$2 != "" {print $2, $1}' ./${species}.pep2ko > level.gene2ko.txt

join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 level3.txt) <(sort -t $'\t' -k1,1 level.gene2ko.txt) > level3.merged.txt
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 level2.txt) <(sort -t $'\t' -k1,1 level.gene2ko.txt) > level2.merged.txt
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 level1.txt) <(sort -t $'\t' -k1,1 level.gene2ko.txt) > level1.merged.txt

awk -F'\t' 'BEGIN {OFS="\t"} {print "P", $2, $3}' level3.merged.txt > level3.final.txt
awk -F'\t' 'BEGIN {OFS="\t"} {print "C", $2, $3}' level2.merged.txt > level2.final.txt
awk -F'\t' 'BEGIN {OFS="\t"} {print "F", $2, $3}' level1.merged.txt > level1.final.txt

cat level*.final.txt | awk '!a[$1$2$3]' > kegg_term.txt

python3 /pathway/script4animal/14.getKEGGnpy.py
mv kegg_dict.npy ${species}_kegg.npy
rm kegg_term.txt level*.txt 11.ko00001.filter4Animal.txt

# =====================================================

awk -vFS="\t" '{if(NF==9){print $0}}' ./${species}.interpro.gff3 | grep "InterPro" | grep "GO:" | awk -vFS=";" -vORS="" '{for(i=1;i<=NF;i++){if($i~/Target=/){split($i, a, " "); print a[1]"\t"}; if($i~/Ontology_term=/){print $i}}; print "\n"}' | sed 's/Target=//; s/Ontology_term=//; s/"//g' | sed "s/,/\t/g" | awk -vFS="\t" -vOFS="\t" '{for(i=2;i<=NF;i++){print $i,$1}}' | cut -f1,3,4 -d"." | sort -u > go.txt

awk -vFS="\t" -vOFS="\t" 'BEGIN{while(getline<"14.GOterm.txt"){a[$1]=a[$1]?a[$1]";"$2"\t"$3:$2"\t"$3}}{if(a[$1]){split(a[$1],rec,";");for(i in rec){print rec[i],$2}}}' go.txt | sort -k1,1 -k2,2n > go_term.txt

python3 /pathway/script4animal/14.getGOnpy.py
mv go_dict.npy ${species}_go.npy
rm go.txt go_term.txt 14.GOterm.txt

done
