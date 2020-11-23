#! /bin/bash -

infile=geo-rnaseq-immune-cells.csv
outDir="../data/rnaseq-raw-studies-series"
error_log=rnaseq-studies-series_download.err
# The script takes a gse number
# and downloads the respective expression files into
# $outDir/<studie>/
# If no arguments is given, it starts one thread of
# itselfe for each gse number in the second column
# of the csv file $infile as arguments.

if [[ $# -eq 0 ]]; then
    # start one thread per line of $infile
    gse_file=$(mktemp -p /dev/shm/)
    trap "rm $gse_file" EXIT
    trap "exit" INT TERM
    export LC_ALL=C
    cut -f2 -d';' "$infile" | tail -n+2 | sort | uniq > $gse_file
    total=$(wc -l "$gse_file" | cut -f1 -d' ')
    parallel --will-cite -a "$gse_file" "$0" 2> "$error_log" | \
        pv -ls $total > /dev/null
    exit
fi

gse="$*"
fl_path="$outDir/$gse"
mkdir -p "$fl_path"
cd "$fl_path"

wget -q "ftp://ftp.ncbi.nlm.nih.gov/geo/series/${gse:0:-3}nnn/$gse/suppl/*"

echo "completed $gse at $time"
