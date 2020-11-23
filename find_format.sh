#! /bin/bash -

# Finds all GSE numbers that seem to have the given format.

name="$1"
test_code="$2"
pre_filter="$3"

error_log="$(pwd -P)/find_$name.err"
inDir="../data/rnaseq-raw-studies"
outFile="$(pwd -P)/format_$name.out"
ensFile="$(pwd -P)/with_ENSG.out"
textFile="$(pwd -P)/text_files.out"
outGSE="$(pwd -P)/${name}_gse.out"

if [[ "$pre_filter" = "ENS" ]]; then
    files="$ensFile"
else
    files="$textFile"
fi

cd "$inDir"

total=$(wc -l "$files" | cut -f1 -d' ')
printf -v prog '%s  && echo "{}" || echo none' "$test_code"
parallel --will-cite -a "$files" "$prog" 2> "$error_log" | \
    pv -ls $total | grep -v "^none$" | tee "$outFile" | \
    grep -o "GSE[0-9]*"  | sort -u > "$outGSE"
