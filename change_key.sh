for file in bib_refs/*.bib*; do
    awk '
    BEGIN { RS="@" ; FS="\n"; ORS=""; }
    NR > 1 {
        last_name = ""; year = "";  # Reset variables for each entry
        for (i = 1; i <= NF; i++) {
            if (tolower($i) ~ /^author\s*=/) {
                split($i, a, "=");
                gsub(/[\{\}]/, "", a[2]);
                split(a[2], authors, ",");
                split(authors[1], name_parts, " ");
                last_name = name_parts[length(name_parts)];
            }
            if (tolower($i) ~ /^year\s*=/) {
                split($i, y, "=");
                gsub(/[\{\}]/, "", y[2]);
                year = y[2];
            }
        }
        if (!last_name) last_name = "UnknownAuthor";
        if (!year) year = "UnknownYear";
        sub(/^(\w+\{)[^,]+/, "\\1" last_name year, $0);
        print "@" $0 "\n";
    }' "$file" > "${file%.bib}_updated.bib"
done
