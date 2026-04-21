# Usage: bash 00.renameGenome.sh namaList

# Extract Genome File and GFF3 File paths from metadata.txt

genome_file=$(cat metadata.txt | grep "Genome File:" | awk -F': ' '{print $2}')
gff3_file=$(cat metadata.txt | grep "GFF3 File:" | awk -F': ' '{print $2}')
repeat_file=$(cat metadata.txt | grep "Repeat File:" | awk -F': ' '{print $2}')
ncRNA_file=$(cat metadata.txt | grep "ncRNA File:" | awk -F': ' '{print $2}')

# 解压文件 (decompress files)
decompress_file() {
    local file=$1
    case $file in
        *.gz)
            gunzip -c "$file" > $(basename "$file" .gz) # 解压并输出为去掉 .gz 后缀的文件名 (decompress and output as the filename without the .gz suffix)
            echo "$(basename "$file" .gz)"
            ;;
        *.zip)
            unzip -q "$file" -d temp_dir
            real_file=$(ls temp_dir | head -n 1)  # 获取解压后的第一个文件 (get the first file after decompression)
            mv "temp_dir/$real_file" "$real_file"
            echo "$real_file"
            ;;
        *.tar)
            tar -xf "$file" -C temp_dir
            real_file=$(ls temp_dir | head -n 1)  # 获取解压后的第一个文件(get the first file after decompression)
            mv "temp_dir/$real_file" "$real_file"
            echo "$real_file"
            ;;
        *.tar.gz)
            tar -xzf "$file" -C temp_dir
            real_file=$(ls temp_dir | head -n 1)  # 获取解压后的第一个文件 (get the first file after decompression)
            mv "temp_dir/$real_file" "$real_file"
            echo "$real_file"
            ;;
        *.bz2)
            bunzip2 -c "$file" > $(basename "$file" .bz2) # 解压并输出为去掉 .bz2 后缀的文件名 (decompress and output as the filename without the .bz2 suffix)
            echo "$(basename "$file" .bz2)"
            ;;
        *)
            echo "$file" # 非压缩文件直接返回 (return the non-compressed file directly)
            ;;
    esac
}

# 解压 genome 文件 (decompress genome file)
genome_real_file=$(decompress_file "$genome_file")
mv "$genome_real_file" genome.fa

# 解压 gff3 文件 (decompress gff3 file)
gff3_real_file=$(decompress_file "$gff3_file")
mv "$gff3_real_file" genome.gff

# 解压 repeat 文件 (decompress repeat file)
if [ -f "$repeat_file" ]; then
	repeat_real_file=$(decompress_file "$repeat_file")
	mv "$repeat_real_file" repeat.gff3
	cp repeat.gff3 repeat.raw.gff3
fi

# 解压 ncRNA 文件 (decompress ncRNA file)
if [ -f "$ncRNA_file" ]; then
	ncRNA_real_file=$(decompress_file "$ncRNA_file")
	mv "$ncRNA_real_file" ncRNA.gff3
	cp ncRNA.gff3 ncRNA.raw.gff3
fi

# ===================================
# 检查 genome.gff 和 genome.fa 是否匹配, 以及gff3文件是否符合规范 (check if genome.gff and genome.fa match, and if the gff3 file is compliant with the specifications)
# 1.必须要有gene的feature，2.必须要有ID和Parent属性，3.必须明确正反链 (1. feature must have gene, 2. must have ID and Parent attributes, 3. must specify strand)
# 4.feature的起始位置必须小于终止位置，5.gff3中长度必须小于等于fasta文件中的长度，6. 序列名必须一致 (4. start position must be less than end position, 5. length in gff3 must be less than or equal to length in fasta file, 6. sequence names must be consistent)
python3 ./00.validate_genome.py genome.gff genome.fa
if [ $? -ne 0 ]; then
    exit 1
fi

# ===================================
# gap number
assembly-stats genome.fa > genome.stats && cat genome.stats
gapNm=$(cat genome.stats | grep -E "^Gaps" | awk -F '[=, ]' '{print $4}')

n90=$(cat genome.stats | grep -E "^N90" | awk -F '[=, ]' '{print $4}')
if [[ $n90 -gt 10000000 ]]; then
    minLength=10000000
else
    minLength=1000000
fi
echo "Selected minLength: $minLength"

python3 - <<EOF
from Bio import SeqIO

def calculate_fa_size(fasta_file, output_file=None):
    """
    Calculate sequence lengths for each chromosome in a FASTA file.
    Args:
        fasta_file (str): Path to the input FASTA file.
        output_file (str, optional): Path to save the genome size results. Defaults to None.
    Returns:
        None
    """
    genome_size = []
    
    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        chromosome = record.id
        length = len(record.seq)
        genome_size.append((chromosome, length))
    
    # Print results
    # print("Chromosome\tLength")
    # for chrom, size in genome_size:
    #    print(f"{chrom}\t{size}")
    
    # Save to a file if output path is provided
    if output_file:
        with open(output_file, "w") as f:
            for chrom, size in genome_size:
                f.write(f"{chrom}\t{size}\n")

# Example usage
calculate_fa_size("genome.fa", "genome.size")
EOF

cat genome.size

chrNm=$(awk -v minLength=${minLength} '{if ($2 >= minLength){print $1}}' genome.size | wc -l)

# 如果gapNm大于chrNm*2，或者大于50个，就认为这个基因组不符合PlanT2T的要求，退出 (If gapNm is greater than chrNm*2, or greater than 50, it is considered that this genome does not meet the requirements of PlanT2T, exit)
if [[ $gapNm -gt $((chrNm*2)) ]] || [[ $gapNm -gt 50 ]]; then
    echo "Too many gaps in the genome, PlanT2T requires that the number of gaps is less than 2 times the number of chromosomes or less than 50."
    exit 1
fi
# ===================================
# 检查蛋白质文件 (Check the protein file)
gffread -y check.pep -g genome.fa genome.gff
if [ ! -f "check.pep" ] || [ ! -s "check.pep" ]; then
    echo "Cannot generate protein sequences from the GFF3 file, please check your GFF3 file."
    exit 1
fi
grep ">" check.pep | cut -d " " -f1 | awk '{ if (a[$1]++) { print "Found duplicate:", $1; exit 1 } }'
if [ $? -ne 0 ]; then
    echo "Exiting due to duplicate sequence names when generating protein sequences"
    exit 1
fi
MidStopAndNoStopFreq=$(python3 ./16.pepstatic.py check.pep | awk '{print ($4+$5)/$2*100}')
if (( $(echo "$MidStopAndNoStopFreq > 10" | bc -l) )); then
    echo "Too many (${MidStopAndNoStopFreq}% > 10%) sequences with internal stop codons or no stop codons in the protein sequences, please check your GFF3 or Genome file."
    exit 1
fi

# ===================================
# 备份用户提供的原始文件 (Backup the original files provided by the user)
cp genome.fa genome.raw.fa
cp genome.gff genome.raw.gff

# ===================================
BeforeFilter=$(wc -l genome.size | cut -d " " -f1)
echo "Chromosome number before filtering: $BeforeFilter"

# 过滤长度至少为minLength的序列，并保存它们的ID (Filter for sequences that are at least minLength long and save their IDs)
awk -v minLength=${minLength} '{if ($2 >= minLength){print $1}}' genome.size | sort -k1,1V > ids.txt && cat ids.txt

# Check the number of sequences after filtering
AfterFilter=$(wc -l ids.txt | cut -d " " -f1)
if [ $AfterFilter -eq 0 ]; then
    echo "Error: All sequences are shorter than the minimum length of $minLength bp."
    exit 1
fi
echo "Chromosome number after filtering: $AfterFilter"

# 提取具有过滤后的ID的序列并保存到新文件中 (Extract sequences with the filtered IDs from the genome file and save to a new file)
seqkit grep -f ids.txt genome.fa > genome.fa.new

# 过滤GFF3文件以仅包含与过滤后的序列对应的条目，并保存到新文件中 (Filter the GFF3 file to include only entries corresponding to the filtered sequences and save to a new file)
awk 'BEGIN{OFS="\t"; FS="\t"} NR==FNR{a[$1]; next} $1 in a {print $0}' ids.txt genome.gff > genome.gff.new

# ===================================
# 查看染色体名字是否以chr开头（不区分大小写），如果不是，则重命名为Chr1,Chr2... (Check if chromosome names start with "chr" (case-insensitive), if not, rename them to "Chr1", "Chr2", etc.)
ChrCount=$(grep -c -E "^[Cc][Hh][Rr]" ids.txt)
if [ $ChrCount -ne $AfterFilter ]; then
    echo "Chromosome names do not start with 'Chr', renaming them to 'Chr1', 'Chr2', etc."
    awk -vFS="\t" -vOFS="\t" '{print $1,"Chr"NR}' ids.txt > ids.txt.new
    mv ids.txt.new ids.txt
    # 对基因组文件中的序列进行重命名 (Rename the sequences in the genome file)
    awk -vFS="\t" -vOFS="\t" 'NR==FNR{a[$1]=$2; next} /^>/{print ">"a[substr($1,2)],$2; next} {print}' ids.txt genome.fa.new > genome.fa
    # 对GFF3文件中的序列进行重命名 (Rename the sequences in the GFF3 file)
    awk -vFS="\t" -vOFS="\t" 'NR==FNR{a[$1]=$2; next} $1 in a {$1=a[$1]; print}' ids.txt genome.gff.new > genome.gff
    rm genome.fa.new genome.gff.new
	
	# 对repeat文件中的序列进行重命名 (Rename the sequences in the repeat file)
	if [ -f "$repeat_file" ]; then
		awk -vFS="\t" -vOFS="\t" 'NR==FNR{a[$1]=$2; next} $1 in a {$1=a[$1]; print}' ids.txt repeat.gff3 > repeat.gff3.new
		mv repeat.gff3.new repeat.gff3
	fi
	
	# 对ncNRA文件中的序列进行重命名 (Rename the sequences in the ncNRA file) 
	if [ -f "$ncNRA_file" ]; then
		awk -vFS="\t" -vOFS="\t" 'NR==FNR{a[$1]=$2; next} $1 in a {$1=a[$1]; print}' ids.txt ncNRA.gff3 > ncNRA.gff3.new
		mv ncNRA.gff3.new ncNRA.gff3
	fi
else
    echo "Chromosome names start with 'Chr', no need to rename."
    mv genome.fa.new genome.fa
    mv genome.gff.new genome.gff
fi
# ===================================

python3 - <<EOF
from Bio import SeqIO

def calculate_fa_size(fasta_file, output_file=None):
    """
    Calculate sequence lengths for each chromosome in a FASTA file.
    Args:
        fasta_file (str): Path to the input FASTA file.
        output_file (str, optional): Path to save the genome size results. Defaults to None.
    Returns:
        None
    """
    genome_size = []
    
    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        chromosome = record.id
        length = len(record.seq)
        genome_size.append((chromosome, length))
    
    # Print results
    # print("Chromosome\tLength")
    # for chrom, size in genome_size:
    #    print(f"{chrom}\t{size}")
    
    # Save to a file if output path is provided
    if output_file:
        with open(output_file, "w") as f:
            for chrom, size in genome_size:
                f.write(f"{chrom}\t{size}\n")

# Example usage
calculate_fa_size("genome.fa", "genome.size")
EOF

cat genome.size

awk -vFS="\t" -vOFS="\t" 'BEGIN{while(getline<"genome.size"){s+=$2}}{print $1,$1,$2,$2/s,"chromosome"}' genome.size > genome.list
rm genome.size ids.txt* check.pep *.fai
