#!/bin/bash

# Usage: ./createfuncgog.sh <obo_file>
# Output: ${DB_NAME}_func_gog

if [ -z "$1" ]; then
    echo "Usage: $0 <obo_file>"
    echo "       DB_NAME environment variable must be set"
    exit 1
fi

if [ -z "$DB_NAME" ]; then
    echo "Error: DB_NAME environment variable is not set"
    echo "  export DB_NAME=your_db_name"
    exit 1
fi

OBO="$1"
OUT="${DB_NAME}_func_gog"

if [ ! -f "$OBO" ]; then
    echo "Error: File not found: $OBO"
    exit 1
fi

echo "Input : $OBO"
echo "Output: $OUT"

awk '
BEGIN {
    RS = ""          # paragraph mode (blank line = record separator)
    FS = "\n"
    OFS = ","
}
{
    # Only process [Term] blocks
    if ($1 != "[Term]") next

    id        = ""
    name      = ""
    namespace = ""
    is_a_ids  = ""
    obsolete  = 0

    for (i = 2; i <= NF; i++) {
        line = $i

        if (line ~ /^id: /)        { id        = substr(line, 5) }
        else if (line ~ /^name: /) { name      = substr(line, 7) }
        else if (line ~ /^namespace: /) { namespace = substr(line, 12) }
        else if (line ~ /^is_obsolete: true/) { obsolete = 1 }
        else if (line ~ /^is_a: /) {
            # is_a: GO:XXXXXXX ! description  →  take only GO:XXXXXXX
            split(line, parts, " ")
            go_id = parts[2]
            if (is_a_ids == "")
                is_a_ids = go_id
            else
                is_a_ids = is_a_ids ";" go_id
        }
    }

    # Skip obsolete terms
    if (obsolete) next

    # Skip terms with no is_a (root terms like GO:0003674, GO:0005575, GO:0008150)
    # Uncomment below if you want to skip them:
    # if (is_a_ids == "") next

    print id, name, is_a_ids, namespace
}
' "$OBO" | sort -t',' -k1,1 > "$OUT"

COUNT=$(wc -l < "$OUT")
echo "Done : $COUNT terms written to $OUT"