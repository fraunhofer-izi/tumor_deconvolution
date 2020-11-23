#! /bin/bash -

# Extract information of the Soft files into a text dictionary.

base_dir=../data/rnaseq-raw-studies-soft
error_log="$(pwd -P)/make-soft-dict.err"

cd "$base_dir"

extractThat() {
    cd "$(dirname "$*")"
    fl="$(basename "$*")"
    out_file="${fl%soft}dic"
    touch "$out_file"
    hg19_count=$(grep -oc "hg19" "$fl")
    printf 'hg19_hits=%d\n' $hg19_count >> "$out_file"
    hg38_count=$(grep -oc "hg38" "$fl")
    printf 'hg38_hits=%d\n' $hg38_count >> "$out_file"
    if (( $hg19_count > $hg38_count )); then
        printf 'genome=hg19\n' >> "$out_file"
    else
        printf 'genome=hg38\n' >> "$out_file"
    fi
    echo "completed $* at $(date)"
}
export -f extractThat

files=$(mktemp -p /dev/shm/)
trap "rm $files" EXIT
trap "exit" INT TERM
find . -type f -name "*soft" > $files
total=$(wc -l "$files" | cut -f1 -d' ')
parallel --will-cite -a "$files" extractThat 2> "$error_log" | \
    pv -ls $total > /dev/null
exit
