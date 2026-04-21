# Usage: bash 08.runInterProScan.sh namaList

# Get the absolute path of the input file
namaList=$(realpath $1)

# Process each line of the input file
cat ${namaList} | cut -f1 | while read i; do
    # Extract protein sequences, remove extensions, and save to a temporary file
    cut -d " " -f1 genome.re.pep | sed '/^[^>]/ s/\.//g' > ${i}.pep

    # Count the number of protein sequences in the pep file
    seq_count=$(grep -c ">" ${i}.pep)

    # If the number of sequences is greater than 50000, split the file
    if [ ${seq_count} -gt 50000 ]; then
        # Calculate the number of split files (each file contains a maximum of 50000 sequences)
        split_count=$((seq_count / 50000 + 1))

        # Split the pep file into smaller files based on sequences using awk
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

        # Process each split file
        for part_file in temp_*.pep; do
            # Run InterProScan analysis for each split file
            /pathway/interproscan-5.72-103.0/interproscan.sh -i ${part_file} -b ${part_file}.interpro -f GFF3 -cpu 30 --iprlookup --goterms -pa -dp
            rm -rf ${part_file} # Remove the temporary split files after processing
        done

        # Merge the analysis results from the split files
        cat temp_*.interpro.gff3 > ${i}.interpro.gff3
        rm -rf temp_*.interpro.gff3 # Remove the individual results from the split files
    else
        # If the number of sequences is less than or equal to 50000, directly run InterProScan
        /pathway/interproscan-5.72-103.0/interproscan.sh -i ${i}.pep -b ${i}.interpro -f GFF3 -cpu 30 --iprlookup --goterms -pa -dp
    fi

    # Clean up temporary files
    rm -rf ${i}.pep temp
done
