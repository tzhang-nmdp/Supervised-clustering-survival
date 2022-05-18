#!/bin/sh

/usr/local/bin/Rscript /Supervised-clustering-survival/R/SCCW_supervised_clustering.R -i /Supervised-clustering-survival/Example/genomic_data_vcf.RData -d  /Supervised-clustering-survival/Example -o test_gene -k 5 
