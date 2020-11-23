#! /bin/bash -

# Finds all none binary files.

inDir="../data/rnaseq-raw-studies"
outFile="$(pwd -P)/text_files.out"

cd "$inDir"

find . -type f -exec grep -Il . {} + | grep -v "annot" | sort > "$outFile"
