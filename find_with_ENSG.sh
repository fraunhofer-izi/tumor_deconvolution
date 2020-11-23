#! /bin/bash -

# Finds files with enseble gene ids.

error_log="$(pwd -P)/find_with_ENSG.err"
inDir="../data/rnaseq-raw-studies"
outFile="$(pwd -P)/with_ENSG.out"

cd "$inDir"

engs_files=$(mktemp -p /dev/shm/)
trap "rm $engs_files" EXIT
trap "exit" INT TERM
export LC_ALL=C
find . -type f ! -name "*annot*" | sort > $engs_files
total=$(wc -l "$engs_files" | cut -f1 -d' ')
parallel --will-cite -a "$engs_files" \
    "head {} | grep -q '^ENSG[0-9]' && echo '{}' || echo none" \
    2> "$error_log" | pv -ls $total | grep -v "^none$" > "$outFile"
