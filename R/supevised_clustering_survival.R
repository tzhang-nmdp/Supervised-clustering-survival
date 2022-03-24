################################################################################################################################################################################# 
# This program is for supervised clustering using two clustering algorithms Kmeans and Mclust, with weighted gene-variant mutation materix by statistical effects and GO similarity, under the cross-validation framework on multiple phenotypes ( mds sbutype, ipssr, hct-ci, overall survival)

# author: Tao Zhang tzhang@nmdp.org

# Example:
# Rscript supervised_clustering_gene_variant-docu-L1-inb-Copy2.R \
# -i  ${genomic_data}.RData \ # input matrix for genomic data
# -o v/g \ # running model option ( 'v' for variant level of common variant analysis, 'g' for gene level of rare variant analysis)
# -d ${outdir} \  # output directory
# -k 5 # k_fold setting for cross-validation

################################################################################################################################################################################
# loading the R packages

library(optparse)
library(cluster)
library(reghelper)
library(bios2mds)
library(reshape2)
library("RColorBrewer")
library(Hmisc)
library(stringr)
library(glmnet)
library(GOSemSim)
library(org.Hs.eg.db)
library('biomaRt')
library(pheatmap)
library(ggplot2)
library(survival)
library(survminer)
library(doParallel)
library(gridExtra)
library(mclust)
library(hash)

################################################################################################################################################################################
# loading the data and parameter settings

option_list<-list(make_option(c("-i", "--input_file"),type="character",help="input file", default=NA,metavar="filename"), make_option(c("-d", "--output_dir"),type="character",help="output dir", default=NA,metavar="filename"), make_option(c("-o", "--output_file"),type="character",help="output file", default=NA,metavar="filename"), make_option(c("-k", "--input_kfolds"),type="character",help="input kfolds", default=NA,metavar="filename"))
opt_parser<-OptionParser(option_list=option_list)
opt=parse_args(opt_parser)

# load genomic data
input<-as.character(opt$input_file)
load(input)

# outdir setting
outdir<-as.character(opt$output_dir)
system(paste('mkdir ', outdir,sep=''))
system(paste('sudo chmod 777 -R ', outdir,sep=''))

# clinical information data
clin_data<-read.table("clin_data.csv",sep="\t",header=T,stringsAsFactors = F,comment.char = "")

# gene/varaint dict file 
all_gene<-read.table("dbNSFP4.0_gene.complete_id",sep="\t",header=F,stringsAsFactors = F,comment.char = "")
variant_gene_id_dict<-read.table("variant_gene_dict",sep="\t",header=F,stringsAsFactors = F,comment.char = "")

# workflow control indicator
opc<-as.character(opt$output_file)
k_folds<-as.numeric(opt$input_kfolds)

# parameter setting
kk_x_list<-2:12
delta_list<-c( 0.01,seq(0,1,by=0.1),0.99)   
km_name<-'genomic_vcf_'

################################################################################################################################################################################
     
# supervised clustering workflow   
if (opc=="g") {                                
    km_name<-paste(km_name,opc,sep="")  
    save.image(file=paste('germ_somatic_vcf_gene_',opc,'.RData' ,sep="") )        
    germ_somatic_vcf_reg<-distance_L1_GO_regulation_g(germ_somatic_vcf, marker_index, -6)    
    germ_somatic_vcf_reg_cluster<-supervised_clustering(germ_somatic_vcf_reg, marker_index, k_folds, marker_cutoff_metrics,km_name)      
    germ_somatic_vcf_reg_cluster_qc<-score_qc_all(outdir,km_name,marker_cutoff_metrics, kk_x_list, delta_list,k_folds)
    save.image(file=paste(outdir,'/',km_name,opc,'.RData' ,sep="") )        
} else if (opc=="v") {                                
    km_name<-paste(km_name,opc,sep="")  
    save.image(file=paste('germ_somatic_vcf_gene_',opc,'.RData' ,sep="") )        
    germ_somatic_vcf_reg<-distance_L1_GO_regulation_v(germ_somatic_vcf, marker_index, -6)    
    germ_somatic_vcf_reg_cluster<-supervised_clustering(germ_somatic_vcf_reg, marker_index, k_folds, marker_cutoff_metrics,km_name)      
    germ_somatic_vcf_reg_cluster_qc<-score_qc_all(outdir,km_name,marker_cutoff_metrics, kk_x_list, delta_list,k_folds)
    save.image(file=paste(outdir,'/',km_name,opc,'.RData' ,sep="") )         
} 
