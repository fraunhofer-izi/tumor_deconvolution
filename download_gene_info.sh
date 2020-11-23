#! /bin/bash

outDir="../data"

cd "$outDir"
wget -c ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
gunzip Homo_sapiens.gene_info.gz

wget -O HUGO.tsv "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_prev_sym&col=gd_aliases&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_ensembl_id&col=gd_vega_ids&col=gd_enz_ids&col=gd_ccds_ids&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
