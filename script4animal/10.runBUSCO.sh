# Usage: bash 10.runBUSCO.sh namaList

/pathway/busco-5.8.2/bin/busco -i genome.re.fa -o busco_result -m genome --offline -c 30 -f \
-l /pathway/busco-5.8.2/buscodb/metazoa_odb10