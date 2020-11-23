#! /bin/bash -

# Shows remaining files with ensemble gene ids and unidendified format.

inFile="$(pwd -P)/with_ENSG.out"

grep -vFf <(cat format_*.out) "$inFile" | grep -v "annot\.txt$"
