#!/bin/sh -e

echo "$MAPPINGFILE"
echo "$1"

FUNC_MAPPING_MODE="${FUNC_MAPPING_MODE:-0}"

if [ "$FUNC_MAPPING_MODE" = "0" ]; then
  # Mode 0: query GO (2-column mapping)
  awk -v db_name="$1" '
  BEGIN {
    outindex=db_name"_func.index";
    outfile=db_name"_func";
    outdb=db_name"_func.dbtype";
    offset=0;

    printf("") > outfile;
    printf("") > outindex;
  }
  FNR==NR {
    id[$1][$2]=1; next
  }

  $2 in id {
    size=0;
    for (key in id[$2]) {
      if (size >0) {
        printf("\n") >> outfile;
        size = size + 1;
      }
      printf("%s", key) >> outfile;
      size=size + length(key);
    }
    printf("\0\n") >> outfile;
    size=size + 2;

    print $1"\t"offset"\t"size >> outindex;
    offset = offset+size;
  }

  END {
    # Gene Ontology datatype
    printf("%c%c%c%c",21,0,0,0) > outdb
  }
  ' "$MAPPINGFILE" "$1.lookup"

elif [ "$FUNC_MAPPING_MODE" = "1" ]; then
  # Mode 1: query GO qualifier namespace evidenceCode evidence (5-6 column mapping,
  # evidence is optional)
  # Mapping file columns: $1=query $2=GO $3=qualifier $4=namespace $5=evidenceCode $6=evidence (optional)
  # Stored entry format:  "GO qualifier namespace evidenceCode evidence" (space-separated,
  # evidence may be empty)
  awk -v db_name="$1" '
  BEGIN {
    outindex=db_name"_func.index";
    outfile=db_name"_func";
    outdb=db_name"_func.dbtype";
    offset=0;

    printf("") > outfile;
    printf("") > outindex;
  }
  FNR==NR {
    evidence=(NF>=6)?$6:"";
    id[$1][$2]=$3" "$4" "$5" "evidence; next
  }

  $2 in id {
    size=0;
    for (key in id[$2]) {
      if (size > 0) {
        printf("\n") >> outfile;
        size = size + 1;
      }
      entry=key" "id[$2][key];
      printf("%s", entry) >> outfile;
      size=size + length(entry);
    }
    printf("\0\n") >> outfile;
    size=size + 2;

    print $1"\t"offset"\t"size >> outindex;
    offset = offset+size;
  }

  END {
    # Gene Ontology datatype
    printf("%c%c%c%c",21,0,0,0) > outdb
  }
  ' "$MAPPINGFILE" "$1.lookup"
fi


# Create Gene Ontology graph file (_func_gog)
if [ -n "$GOOBO" ] && [ -f "$GOOBO" ]; then
    awk '
    BEGIN {
        RS = ""
        FS = "\n"
        OFS = ","
    }
    {
        if ($1 != "[Term]") next

        id        = ""
        name      = ""
        namespace = ""
        is_a_ids  = ""
        obsolete  = 0

        for (i = 2; i <= NF; i++) {
            line = $i

            if (line ~ /^id: /)             { id        = substr(line, 5) }
            else if (line ~ /^name: /)      { name      = substr(line, 7) }
            else if (line ~ /^namespace: /) { namespace = substr(line, 12) }
            else if (line ~ /^is_obsolete: true/) { obsolete = 1 }
            else if (line ~ /^is_a: /) {
                split(line, parts, " ")
                go_id = parts[2]
                if (is_a_ids == "")
                    is_a_ids = go_id
                else
                    is_a_ids = is_a_ids ";" go_id
            }
        }

        if (obsolete) next

        print id, name, is_a_ids, namespace
    }
    ' "$GOOBO" | sort -t',' -k1,1 > "$1_func_gog"
fi
# fail() {
#     echo "Error: $1"
#     exit 1
# }

# notExists() {
# 	  [ ! -e "$1" ]
# }

# hasCommand () {
#     command -v "$1" >/dev/null 2>&1
# }

# notExists "$1" && echo "$1 not found!" && exit 1;

# hasCommand awk
# hasCommand gunzip
# hasCommand touch
# hasCommand tar

# TAXDBNAME="$1"
# TMP_PATH="$2"

# STRATEGY=""
# if hasCommand aria2c; then STRATEGY="$STRATEGY ARIA"; fi
# if hasCommand curl;   then STRATEGY="$STRATEGY CURL"; fi
# if hasCommand wget;   then STRATEGY="$STRATEGY WGET"; fi
# if [ "$STRATEGY" = "" ]; then
#     fail "No download tool found in PATH. Please install aria2c, curl or wget."
# fi

# downloadFile() {
#     URL="$1"
#     OUTPUT="$2"
#     set +e
#     for i in $STRATEGY; do
#         case "$i" in
#         ARIA)
#             FILENAME=$(basename "${OUTPUT}")
#             DIR=$(dirname "${OUTPUT}")
#             if aria2c -c --max-connection-per-server="$ARIA_NUM_CONN" --allow-overwrite=true -o "${FILENAME}.aria2" -d "$DIR" "$URL"; then
#                 mv -f -- "${OUTPUT}.aria2" "${OUTPUT}"
#                 return 0
#             fi
#             ;;
#         CURL)
#             if curl -L -C - -o "${OUTPUT}.curl" "$URL"; then
#                 mv -f -- "${OUTPUT}.curl" "${OUTPUT}"
#                 return 0
#             fi
#             ;;
#         WGET)
#             if wget -O "${OUTPUT}.wget" -c "$URL"; then
#                 mv -f -- "${OUTPUT}.wget" "${OUTPUT}"
#                 return 0
#             fi
#             ;;
#         esac
#     done
#     set -e
#     fail "Could not download $URL to $OUTPUT"
# }

# if { [ "${DBMODE}" = "1" ] && notExists "${TAXDBNAME}_taxonomy"; } || { [ "${DBMODE}" = "0" ] && { notExists "${TAXDBNAME}_names.dmp" || notExists "${TAXDBNAME}_nodes.dmp" || notExists "${TAXDBNAME}_merged.dmp"; }; }; then
#     if [ "$DOWNLOAD_NCBITAXDUMP" -eq "1" ]; then
#         # Download NCBI taxon information
#         if notExists "${TMP_PATH}/ncbi_download.complete"; then
#             echo "Download taxdump.tar.gz"
#             downloadFile "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" "${TMP_PATH}/taxdump.tar.gz"
#             tar -C "${TMP_PATH}" -xzf "${TMP_PATH}/taxdump.tar.gz" names.dmp nodes.dmp merged.dmp delnodes.dmp
#             touch "${TMP_PATH}/ncbi_download.complete"
#             rm -f "${TMP_PATH}/taxdump.tar.gz"
#         fi
#         NCBITAXINFO="${TMP_PATH}"
#     fi
#     if [ "${DBMODE}" = "1" ]; then
#         # shellcheck disable=SC2086
#         "${MMSEQS}" createbintaxonomy "${NCBITAXINFO}/names.dmp" "${NCBITAXINFO}/nodes.dmp" "${NCBITAXINFO}/merged.dmp" "${TAXDBNAME}_taxonomy" ${VERBOSITY_PAR} \
#             || fail "createbintaxonomy failed"
#     else
#         cp -f "${NCBITAXINFO}/names.dmp"    "${TAXDBNAME}_names.dmp"
#         cp -f "${NCBITAXINFO}/nodes.dmp"    "${TAXDBNAME}_nodes.dmp"
#         cp -f "${NCBITAXINFO}/merged.dmp"   "${TAXDBNAME}_merged.dmp"
#         cp -f "${NCBITAXINFO}/delnodes.dmp" "${TAXDBNAME}_delnodes.dmp"
#     fi
# fi

# if notExists "${TAXDBNAME}_mapping"; then
#     if [ "$DOWNLOAD_MAPPING" -eq "1" ]; then
#         # Download the latest UniProt ID mapping to extract taxon identifiers
#         if notExists "${TMP_PATH}/mapping_download.complete"; then
#             echo "Download idmapping.dat.gz"
#             URL="https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"
#             downloadFile "$URL" "${TMP_PATH}/idmapping.dat.gz"
#             gunzip -c "${TMP_PATH}/idmapping.dat.gz" | awk '$2 == "NCBI_TaxID" {print $1"\t"$3 }' > "${TMP_PATH}/taxidmapping"
#             touch "${TMP_PATH}/mapping_download.complete"
#             rm -f "${TMP_PATH}/idmapping.dat.gz"
#         fi
#         MAPPINGFILE="${TMP_PATH}/taxidmapping"
#     fi

#     if [ "$MAPPINGMODE" = "0" ]; then
#         awk 'NR == FNR { f[$1] = $2; next } $2 in f { print $1"\t"f[$2] }' \
#             "$MAPPINGFILE" "${TAXDBNAME}.lookup" > "${TAXDBNAME}_mapping"
#     else
#         awk 'FNR == 1 { fidx++; } fidx == 1 { tax[$1] = $2; next; } fidx == 2 { source[$1] = tax[$2]; next; } fidx == 3 { print $1"\t"source[$3]; next; }' \
#             "$MAPPINGFILE" "${TAXDBNAME}.source" "${TAXDBNAME}.lookup" > "${TAXDBNAME}_mapping"
#     fi
# fi

# if [ -n "$REMOVE_TMP" ]; then
#     if [ "$DOWNLOAD_NCBITAXDUMP" -eq "1" ]; then
#         rm -f "${TMP_PATH}/names.dmp" "${TMP_PATH}/nodes.dmp" "${TMP_PATH}/merged.dmp" "${TMP_PATH}/delnodes.dmp"
#         rm -f "${TMP_PATH}/ncbi_download.complete"
#     fi
#     if [ "$DOWNLOAD_MAPPING" -eq "1" ]; then
#         rm -f "${TMP_PATH}/taxidmapping" "${TMP_PATH}/mapping_download.complete"
#     fi
#     rm -f "${TMP_PATH}/createtaxdb.sh"
# fi
