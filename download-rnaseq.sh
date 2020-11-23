#! /bin/bash -

infile=geo-rnaseq-immune-cells.csv
outDir="../data/rnaseq-raw"
# The script takes a line of $infile
# and downloads the respective expression file into
# ../data/rnaseq-raw/<cell.type>/
# If no arguments is given, it starts one thtead of
# itselfe with the lines of $infile as arguments.


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
curl --data-urlencode "file=$file" --data-urlencode "acc=$id" \
    --data-urlencode "format=file" https://www.ncbi.nlm.nih.gov/geo/download/ \
    --output "${file%.gz}" -s -S --compressed

echo "completed $id at $time"
