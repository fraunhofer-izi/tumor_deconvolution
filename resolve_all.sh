#! /bin/bash -

source /etc/profile.d/modules.sh
module load R/3.5.2-0

export inPath=../data/rnaseq-raw-studies/
export outPath=../data/rnaseq-R/
error_log=resolve_all.err

# The script attempts to resolve all GSE with the R script.

mkdir -p "$outPath"

gse_file=$(mktemp -p /dev/shm/)
export jobtag=$(printf "gse extract %s" "$(date)" | sha256sum | head -c 10)
trap_task() {
    >&2 printf 'Cleaning up ...\n'
    rm $gse_file
    >&2 printf 'Killing all associated slurm jobs ...\n'
    sleep 5 &
    squeue -h -o "%i;%k" | grep ";jobtag $jobtag$" | cut -f1 -d';' | \
        xargs -P0 -n1 scancel 2> /dev/null
    # some jobs may still beeing submitted
    wait
    squeue -h -o "%i;%k" | grep ";jobtag $jobtag$" | cut -f1 -d';' | \
        xargs -n1 scancel 2> /dev/null
}
trap trap_task EXIT
trap "exit 1" INT TERM HUP
export LC_ALL=C
total=$(find "$inPath" -mindepth 1 -maxdepth 1 -type d | \
        sed 's|.*/||g' | sort | tee "$gse_file" | wc -l)

mem_gse(){
    # returns GSE specific memory requirements
    GSE="$1"
    large_gse=( GSE100501 GSE113675 GSE117156 GSE87254 GSE112658
        GSE119428 GSE107185 GSE88933 GSE116672 GSE23316 GSE110612 GSE76270 )
    for gse in "${large_gse[@]}"; do
        [[ "$gse" == "$GSE" ]] && printf '240G' && return
    done
    printf '35G'
}
export -f mem_gse

cores_gse(){
    # returns GSE specific CPU requirements
    GSE="$1"
    large_gse=( GSE100501 GSE113675 GSE117156 GSE87254 GSE112658
        GSE119428 GSE107185 GSE88933 GSE116672 GSE23316 GSE110612 GSE76270 )
    for gse in "${large_gse[@]}"; do
        [[ "$gse" == "$GSE" ]] && printf '16' && return
    done
    printf '1'
}
export -f cores_gse

bash_job(){
    GSE="$1"
    if R/resolve_general.R "$inPath/$GSE" $(cores_gse "$GSE") &> "$outPath/$GSE.log"
    then
        printf "=%b" "\x00"; >&2 printf "Success finish: %s\n" "$GSE"
    else
        printf "f%b" "\x00"; >&2 printf "Error finish: %s (%s)\n" "$GSE" "$outPath/$GSE.log"
    fi
}
export -f bash_job

slurm_job(){
    GSE="$1"
    srun --quiet --job-name="extract $GSE" --mem=$(mem_gse "$GSE") \
        --cpus-per-task=$(cores_gse "$GSE") --time=0-05:00:00 \
        --comment="jobtag $jobtag" bash -c "bash_job $GSE"
}
export -f slurm_job

cat $gse_file | xargs -n1 -I{} -P0 bash -c 'slurm_job {}' 2> "$error_log" | \
    pv -0ls $total > /dev/null

rm $gse_file
trap - EXIT

printf 'Good extractions: %s ' $(grep -a "^Success finish: " "$error_log" | wc -l)
printf 'Bad extractions: %s ' $(grep -a "^Error finish: " "$error_log" | wc -l)
printf 'Slurm canceled extractions: %s\n' $(grep -a "^slurmstepd.*CANCELLED" "$error_log" | wc -l)
printf 'Read %s for more information!\n' "$error_log"
