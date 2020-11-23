#! /bin/bash -

##############################################################################
#
#    Copyright © 2019 Gesellschaft zur Förderung der angewandten Forschung e.V.
#    acting on behalf of its Fraunhofer Institut für Zelltherapie und Immunologie.
#    All rights reserved. Contact: dominik.otto@izi.fraunhofer.de
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License as published by the Free Software
#    Foundation; either version 3 of the License, or (at your option) any later
#    version.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT ANY
#    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#    PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along with
#    this program; if not, see <http://www.gnu.org/licenses/>.
#
##############################################################################

START=$(date +%s)
inDir="/input"
inFile="$inDir/input.csv"
[ "$MODE" = "debug" ] && CORES=16 || CORES=4
[ "$MODE" = "debug" ] && TIME=43200 || TIME=172800
export CORES
export TIME
export outDir="/deconvs"
export outFile="/output/predictions.csv"
export grain_file="grains"
export python_script="deconvolute.py"
export characteristics="chars.pkl"
export time_form="+%Y-%m-%d %H:%M:%S"
export sample_count_file="samplecount"
export END=$(( $START + $TIME ))

mkdir -p "$outDir"
mkdir -p "$(dirname "$outFile")"

printf 'File: %s\n' "$inFile"
cat "$inFile"
printf '\n'

declare -A colmap
export header="$(head -n1 "$inFile")"
colmap() {
    # because hashmaps cannot be exported we emulate it with a function
    i=0
    printf '%s,' "$header" | while read -d, colname; do
        [[ "${1^^}" == "${colname^^}" ]] && printf "$i" && break
        ((i++))
    done
}
export -f colmap

make_empty(){
    set_name="$1"
    sample_name="$2"
    out_file="$outDir/$set_name.$sample_name.deconv"
    cat "$grain_file" | xargs -n1 -I{} printf '%s,%s,%s,%f\n' \
        "$set_name" "$sample_name" "{}" "0.0" > "$out_file"
}
export -f make_empty

tail -n+2 "$inFile" | while read line; do
    IFS=',' read -ra info <<< "$line,"
    export set_name="${info[$(colmap 'dataset.name')]}"
    file="$inDir/${info[$(colmap 'hugo.expr.file')]}"
    printf '%s,' $(head -n1 "$file") | \
        tr ',' '\n' | tr '\r' '\n' | tail -n+2 | \
        xargs -n1 -I{} bash -c 'make_empty "$set_name" "{}"'
done
ls -1 "$outDir" | wc -l > "$sample_count_file"

reduce_samples(){
    (
        flock -x 200
        samples=$(<$sample_count_file)
        printf '%i' $(($samples-1)) > $sample_count_file
    ) 200<> $sample_count_file
}
export -f reduce_samples

update(){
    cp header "$outFile"
    chmod a+rw "$outFile"
    find "$outDir" -type f -name "*.deconv" -exec cat {} + >> "$outFile"
    if [ "$MODE" = "debug" ]; then
        printf '[%s] updated %s, %i remaining\n' \
            "$(date -u "$time_form")" "$outFile" "$(cat $sample_count_file)"
    fi
}
export -f update

auto_update(){
    trap 'update; exit' INT TERM HUP
    update
    while :; do
        inotifywait -q "$outDir" > /dev/null
        update
    done
}
auto_update &
export upadtePID=$!
trap '>&2 printf "Interrupted..."; exit 130' INT
trap '>&2 printf "Terminated ...\n"; exit 143' TERM
trap '>&2 printf "Connection to user lost ...\n"; exit 129' HUP
trap 'ec=$?; kill $upadtePID; exit $ec' EXIT

run_deconv(){
    set_name="$1"
    set_file="$2"
	scale="$3"
	norm="$4"
    cancer_type="$5"
    platform="$6"
    sample_name="$7"
    out_file="$outDir/$set_name.$sample_name.deconv"
    >&2 printf '[%s] Starting %s: %s %s...\n' "$(date -u "$time_form")" \
        "$SLOT" "$set_name" "$sample_name"
    printf -v cmd 'python -u "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s"' \
        "$python_script" "$characteristics" "$set_name" \
        "$set_file" "$sample_name" "$scale" "$norm" "$cancer_type" "$platform" \
        "$out_file" "$(dirname "$outFile")"
    if [ "$MODE" = "debug" ]; then
        printf '[Slot %s] Runnig: %s\n' "$SLOT" "$cmd"
    fi
    if [[ -z "$SLOT" ]]; then
        eval "$cmd" # outputs progress
    else
        2>&1 eval "$cmd" | \
            stdbuf -oL -eL sed "s/^/[Slot $SLOT: $set_name $sample_name] /g"
    fi
    reduce_samples
    >&2 printf '[%s] Finished %s: %s %s...\n' "$(date -u "$time_form")" \
        "$SLOT" "$set_name" "$sample_name"
}
export -f run_deconv

if (($CORES==1)); then
    slotvar=SINGLE
else
    slotvar=SLOT
fi
tail -n+2 "$inFile" | sed 's/,\s*,/,NA,/g' | sed 's/,\s*,/,NA,/g' | \
    while read line; do
    IFS=',' read -ra info <<< "$line,"
    set_name="${info[$(colmap 'dataset.name')]// /_}"
    c_type="${info[$(colmap 'cancer.type')]// /_}"
	scale="${info[$(colmap 'scale')]// /_}"
    norm="${info[$(colmap 'normalization')]// /_}"
    platform="${info[$(colmap 'platform')]// /_}"
    platform="${platform//[^a-zA-Z0-9_]/}"
    [ $FEATURES = "HUGO" ] \
        && file="$inDir/${info[$(colmap 'hugo.expr.file')]// /_}" \
        || file="$inDir/${info[$(colmap 'ensg.expr.file')]// /_}"
    printf '%s,' $(head -n1 "$file") | \
        tr ',' '\n' | tr '\r' '\n' | tail -n+2 | xargs -n1 \
        printf '"%s" "%s" "%s" "%s" "%s" "%s" "%s"\n' \
        "$set_name" "$file" "$scale" "$norm" "$c_type" "$platform"
done | xargs -n1 -I{} -P$CORES --process-slot-var=$slotvar bash -c 'run_deconv {}'

