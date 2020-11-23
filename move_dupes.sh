#! /bin/bash -

# renames duplicate <format> samples

format="$1"
error_log="$(pwd -P)/move_${format}_dupes.err"
inDir="../data/rnaseq-raw-studies"
inFile="$(pwd -P)/${format}.out"

cd "$inDir"

dupes=$(mktemp -p /dev/shm/)
trap "rm $dupes" EXIT
trap "exit" INT TERM
grep -o "GSM[0-9]*" "${format}.out" | sort | uniq -c | awk '$1>1 {print $2}' > $dupes
total=$(wc -l "$dupes" | cut -f1 -d' ')

moveThat(){
    path="$*"
    newpath="${path/GSM/dup-GSM}"
    if [[ ! -f "$newpath" ]] || [[ -f "$path" ]]; then
        mv "$path" "$newpath"
        echo "Moved $path"
    else
        echo "Unmoved $path"
    fi
}
export -f moveThat

parallel --will-cite -a "$txts" \
    "moveThat \"$(grep {} "$inFile" | head -n1)\"" 2> "$error_log" | \
    pv -ls $total > /dev/null
