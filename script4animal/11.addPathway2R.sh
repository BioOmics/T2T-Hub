#!/bin/bash

read -r -d '' CODE_TO_ADD << EOM
    ## Load pathway data 'ko00001.T2THub.txt' as 'pathway4Animal'
    data_file <- system.file("extdata", "ko00001.T2THub.txt", package = pkgname)
    if (file.exists(data_file)) {
        pathway4Animal <- read.table(data_file, header = TRUE, sep = "takeplacename", stringsAsFactors = F, quote = "", check.names = F)
        assign("pathway4Animal", pathway4Animal, envir = .GlobalEnv)
        packageStartupMessage("Loaded 'pathway4Animal' from ko00001.T2THub.txt for KEGG annotation")
		packageStartupMessage(" ")
		packageStartupMessage("See more information in T2T-Hub:")
		packageStartupMessage("NJU site: https://biobigdata.nju.edu.cn/t2thub/")
		packageStartupMessage("ZJU site: https://bis.zju.edu.cn/t2thub/")
    } else {
        packageStartupMessage("Pathway data 'ko00001.T2THub.txt' not found.")
		packageStartupMessage("If you have any questions, please open an issue on GitHub")
		packageStartupMessage("GitHub site: https://github.com/BioOmics/T2T-Hub")
    }
EOM

OrgDB=${1}
ZZZ_R_FILE="${OrgDB}/R/zzz.R"

if [[ ! -f "$ZZZ_R_FILE" ]]; then
    echo "Error: $ZZZ_R_FILE does not exist."
    exit 1
fi

onLoadStartLine=$(grep -n '.onLoad <- function' "$ZZZ_R_FILE" | cut -d ':' -f 1)
onLoadBraceLine=$(awk "/^ *{/{print NR; exit}" "$ZZZ_R_FILE")

if [[ -z "$onLoadStartLine" ]] || [[ -z "$onLoadBraceLine" ]]; then
    echo "Error: .onLoad function or its opening brace not found in $ZZZ_R_FILE."
    exit 1
fi

if grep -q "Load pathway data 'ko00001.T2THub.txt'" "$ZZZ_R_FILE"; then
    echo "The code block already exists in $ZZZ_R_FILE. Skipping insertion."
else
    awk -v code="$CODE_TO_ADD" -v braceLine="$onLoadBraceLine" '
    NR == braceLine { print; print code; next }
    { print }
    ' "$ZZZ_R_FILE" > "${ZZZ_R_FILE}.tmp" && mv "${ZZZ_R_FILE}.tmp" "$ZZZ_R_FILE"
    
    echo "Code block inserted into $ZZZ_R_FILE."

    sed -i 's/takeplacename/\\t/g' "$ZZZ_R_FILE"
fi
