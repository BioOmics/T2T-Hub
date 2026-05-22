# Usage:
# bash 09.runKoFamScan.sh namaList

namaList=$(realpath $1)

cat ${namaList} | cut -f1 | while read i
do

    if [ ! -s "${i}.pep2ko" ]; then

        ########################################
        # 1. 准备 pep 文件
        ########################################
        cp genome.re.pep ${i}.pep

        ########################################
        # 2. 检查是否存在 >100k 序列
        ########################################
        awk '
        BEGIN{
            RS=">";
            FS="\n"
        }
        NR>1{

            len=0

            for(i=2;i<=NF;i++){
                len += length($i)
            }

            split($1,a," ")
            id=a[1]

            if(len > 100000){
                print id
            }
        }' ${i}.pep > ${i}.tooLong.id

        ########################################
        # 3. 根据情况决定输入文件
        ########################################
        if [ -s ${i}.tooLong.id ]; then

            echo "[INFO] Found sequences >100k in ${i}"

            seqkit grep -v -f ${i}.tooLong.id \
                ${i}.pep > ${i}.filtered.pep

            inputPep=${i}.filtered.pep

        else

            inputPep=${i}.pep

        fi

        ########################################
        # 4. 跑 KoFamScan
        ########################################
        exec_annotation \
            -o ${i}.pep2ko.raw.txt \
            -p /public/workspace/biobigdata/project/Plant2t/software/kofamscan/profiles \
            -k /public/workspace/biobigdata/project/Plant2t/software/kofamscan/ko_list \
            --cpu 20 \
            -E 1e-5 \
            -f mapper \
            --report-unannotated \
            ${inputPep} \
            --tmp-dir /tmp/kotmp_${i}

        ########################################
        # 5. 整理正常注释结果
        ########################################
        sed 's/\.[0-9]*//g' ${i}.pep2ko.raw.txt \
            | awk '!a[$1]++' \
            > ${i}.pep2ko.tmp

        ########################################
        # 6. 若存在超长序列，则补空注释
        ########################################
        if [ -s ${i}.tooLong.id ]; then

            sed 's/\.[0-9]*//g' ${i}.tooLong.id \
                | awk '{print $1"\t"}' \
                > ${i}.tooLong.empty

            cat ${i}.pep2ko.tmp ${i}.tooLong.empty \
                | awk '!a[$1]++' \
                > ${i}.pep2ko

        else

            mv ${i}.pep2ko.tmp ${i}.pep2ko

        fi

        ########################################
        # 7. 清理临时文件
        ########################################
        rm -rf \
            ${i}.pep \
            ${i}.filtered.pep \
            ${i}.pep2ko.raw.txt \
            ${i}.pep2ko.tmp \
            ${i}.tooLong.id \
            ${i}.tooLong.empty \
            /tmp/kotmp_${i}

    fi

    ########################################
    # 8. 检查结果
    ########################################
    if [ ! -s "${i}.pep2ko" ]; then

        echo "File empty: ${i}.pep2ko"
        exit 1

    else

        echo "Success: ${i}.pep2ko"

    fi

done
