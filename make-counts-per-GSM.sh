#! /bin/bash -

infile=geo-rnaseq-immune-cells.csv
outDir="../data/rnaseq-counts"
# This script attempts to generate a count
# vector per GSM number in $infile.

if [[ $# -eq 0 ]]; then
    # start one thread per line of $infile
    total=$(($(wc -l "$infile" | cut -f1 -d' ') - 1))
    parallel --will-cite --skip-first-line -a "$infile" "${0};0" | \
        pv -ls $total > /dev/null
    exit
fi

IFS=';' read -ra fields <<< "$@"
id="${fields[3]}"
file="${fields[2]%.anno.tsv}"
cell_type="${fields[7]}"

fl_path="$outDir/$cell_type"
mkdir -p "$fl_path"
cd "$fl_path"

# do the extration

echo "completed $id at $time"
