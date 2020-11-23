#! /bin/bash -

# Unpacks all .tar, .gz and .tar.gz in the given directory.
unpacker="$(readlink -f ${BASH_SOURCE[0]})"
error_log="$(dirname "$unpacker")/unpack_run.err"
inDir="../data/rnaseq-raw-studies"
errorDir="../data/rnaseq-raw-studies-errors"
export errorDir="$(dirname "$unpacker")/$errorDir"

if [[ $# -eq 0 ]]; then
    printf 'Unpacking %s...\n' "$(readlink -f "$inDir")" > "$error_log"
    printf 'Crawl 1:\n' >> "$error_log"
    cd "$inDir"
    run_number=1
else
    base_dir="$1"
    run_number="$2"
    cd "$base_dir"
    printf 'Crawl %i:\n' "$run_number" >> "$error_log"
fi

echo "Unpacking crawl number $run_number:"

unpackThat() {
    log=$(mktemp -p /dev/shm/)
    trap "rm $log" EXIT
    trap "exit" INT TERM
    exec 2> >(tee -a "$log" >&2)
    inDir="$(dirname "$*")"
    cd "$inDir"
    fl="$(basename "$*")"
    dump() {
        errorDump="$errorDir/$(basename "$inDir")"
        dumpLog="$errorDump/$fl.log"
        mkdir -p "$errorDump"
        mv "$fl" "$errorDump/$fl"
        cp "$log" "$dumpLog"
    }
    if [[ "$*" =~ .tar.gz$ ]] || [[ "$*" =~ .tgz$ ]]; then
        tar -xzf "$fl" && rm "$fl" || dump
    elif [[ "$*" =~ .gz$ ]]; then
        if [[ -z "$(du --threshold=2G "$fl")" ]]; then
            gzip -df "$fl" || dump
        else
            pigz -df "$fl" || dump
        fi
    elif [[ "$*" =~ .tar$ ]]; then
        tar -xf "$fl" && rm "$fl" || dump
    fi
    echo "completed $* at $(date)"
}
export -f unpackThat

packed=$(mktemp -p /dev/shm/)
trap "rm $packed" EXIT
trap "exit" INT TERM
find . -type f \( -name "*tar" -o -name "*gz" -o -name "*tgz" \) | sort -r > $packed
total=$(wc -l "$packed" | cut -f1 -d' ')
parallel --will-cite -a "$packed" unpackThat 2>> "$error_log" | \
    pv -ls $total  > /dev/null

if (( "$total" > 0 )); then
    $unpacker "$base_dir" $(( $run_number + 1 ))
fi
