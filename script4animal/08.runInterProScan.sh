# Usage: bash 08.runInterProScan.sh namaList

# Get the absolute path of the input file (获取输入文件的绝对路径)
namaList=$(realpath $1)

# Process each line of the input file (逐行处理输入文件)
cat ${namaList} | cut -f1 | while read i; do
    # Extract protein sequences, remove extensions, and save to a temporary file (提取蛋白质序列，去除扩展名，并保存到临时文件)
    cut -d " " -f1 genome.re.pep | sed '/^[^>]/ s/\.//g' > ${i}.pep

    # Count the number of protein sequences in the pep file (统计pep文件中的蛋白质序列数量)
    seq_count=$(grep -c ">" ${i}.pep)

    # If the number of sequences is greater than 50,000, split the file and process in batches (如果序列数量大于50,000，则分割文件并批量处理)
    if [ ${seq_count} -gt 50000 ]; then
        # Calculate the number of split files (each file contains a maximum of 50,000 sequences) (计算分割文件的数量，每个文件包含最多50,000个序列)
        split_count=$((seq_count / 50000 + 1))

        # Split the pep file into smaller files based on sequences using awk (根据序列使用awk将pep文件分割成更小的文件)
        awk -v max_seqs=50000 '
        BEGIN {seq_num=0; file_num=1; out_file="";} 
        /^>/ {
            seq_num++; 
            if (seq_num > max_seqs) {
                seq_num=1;
                file_num++;
            }
            if (seq_num == 1) {
                close(out_file);
                out_file = "temp_" file_num ".pep";
            }
            print > out_file;
        }
        /^[^>]/ { print > out_file; }
        ' ${i}.pep

        # Process each split file with InterProScan (使用InterProScan处理每个分割文件)
        for part_file in temp_*.pep; do
            # Run InterProScan analysis for each split file (对每个分割文件运行InterProScan分析)
            /pathway/interproscan-5.72-103.0/interproscan.sh -i ${part_file} -b ${part_file}.interpro -f GFF3 -cpu 30 --iprlookup --goterms -pa -dp
            rm -rf ${part_file} # Remove the temporary split files after processing (处理后删除临时分割文件)
        done

        # Merge the analysis results from the split files
        cat temp_*.interpro.gff3 > ${i}.interpro.gff3
        rm -rf temp_*.interpro.gff3 # Remove the individual results from the split files (删除分割文件的单独结果)
    else
        # If the number of sequences is less than or equal to 50,000, directly run InterProScan analysis (如果序列数量小于或等于50,000，则直接运行InterProScan分析)
        /pathway/interproscan-5.72-103.0/interproscan.sh -i ${i}.pep -b ${i}.interpro -f GFF3 -cpu 30 --iprlookup --goterms -pa -dp
    fi

    # Clean up temporary files (清理临时文件)
    rm -rf ${i}.pep temp
done
