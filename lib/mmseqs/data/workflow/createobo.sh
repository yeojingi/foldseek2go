#!/bin/sh -e

echo "$GOOBO"

awk -v db_name="$1" '
BEGIN {
  outfile = db_name "_func_gog";
  goid = "";
  goname = "";
  goparents = "";
  ns = "";
  printf("") > outfile;
}

/^\[Term\]/ {
  # 이전 항목 출력
  if (goid != "" && is_obsolete != 1) {
    print goid","goname","goparents","ns >> outfile;
  }
  # 초기화
  goid = "";
  goname = "";
  goparents = "";
  ns = "";
  is_obsolete = 0;
}

/^id: GO:/ {
  goid = $2;
}

/^name:/ {
  $1 = ""; sub(/^ /, ""); goname = $0;
}

/^namespace:/ {
  ns = $2;
}

/^(is_a|relationship: part_of)/ {
  match($0, /GO:[0-9]+/, arr);
  if (arr[0] != "") {
    if (goparents != "") {
      goparents = goparents ";" arr[0];
    } else {
      goparents = arr[0];
    }
  }
}

/^is_obsolete: true/ {
  is_obsolete = 1;
}

END {
  if (goid != "" && is_obsolete != 1) {
    print goid","goname","goparents","ns >> outfile;
  }
}
' "$GOOBO"