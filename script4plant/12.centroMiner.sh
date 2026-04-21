# Usage: bash 12.centroMiner.sh namaList

# 检查是否有染色体的长度大于 100M（100,000,000）(check if any chromosome length exceeds 100M)
over_100M=$(awk '{if ($2 > 100000000) {print "yes"; exit}}' genome.re.fa.fai)

# 根据判断结果执行不同的 quartet.py 命令 (execute different quartet.py commands based on the result)
if [ "$over_100M" == "yes" ]; then
    # 如果有染色体大于 100M，执行带 --noplot 参数的命令 (execute the command with --noplot if any chromosome exceeds 100M)
	echo "Executing CentroMiner with no plot generation (--noplot) as some chromosomes exceeds 100M in size."
    quartet.py CentroMiner -i genome.re.fa --gene genome.renamed.gff -t 30 -r 100 -l 10000 --noplot
else
    # 如果没有染色体大于 100M，执行不带 --noplot 参数的命令 (execute the command without --noplot if no chromosome exceeds 100M)
    quartet.py CentroMiner -i genome.re.fa --gene genome.renamed.gff -t 30 -r 100 -l 10000
fi

echo -e "Chrid\tStart\tEnd\tLength\tTRlength\tTRcoverage" > genome.centromere.info
length_file="genome.telo.info"
declare -A chr_lengths
while read -r line; do
  chr=$(echo "$line" | cut -f1)
  length=$(echo "$line" | cut -f2)
  chr_lengths[$chr]=$length
done < <(grep -v "#" "$length_file")
chr_length_str=$(for chr in "${!chr_lengths[@]}"; do echo "$chr:${chr_lengths[$chr]}"; done | tr '\n' ',')

centromere_file="./Candidates/*.candidate"
cat $centromere_file | grep -v -e "#" -e "@" | sort -k1,1V -k5,5nr | \
awk -vOFS="\t" -v lengths="$chr_length_str" 'BEGIN {
  split(lengths, chr_lengths, ",");
  for (i in chr_lengths) {
    split(chr_lengths[i], arr, ":");
    length_map[arr[1]] = arr[2];
  }
} {
  if ($2 > 100000 && $3 < length_map[$1] - 100000) {
    print $0;
  }
}' |  awk '!a[$1]++' | awk -vOFS="\t" '{print $1,$2,$3,$4,$5,$6}' >> genome.centromere.info

column -t genome.centromere.info > genome.centromere.col.info
mv genome.centromere.col.info genome.centromere.info

tar -zcf genome.Candidates.tar.gz Candidates

