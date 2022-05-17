################################################################################################################################################################################# 
# This program is for supervised clustering using two clustering algorithms Kmeans and Mclust, with weighted gene-variant mutation materix by statistical effects and GO similarity, under the cross-validation framework on multiple metrics (overall survival /clustering)

# author: Tao Zhang tzhang@nmdp.org

################################################################################################################################################################################
# loading the R packages

library(optparse)
library(dplyr)
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
clin_data<-read.table("dat.csv.cr",sep="\t",header=T,stringsAsFactors = F,comment.char = "")
clin_data_kn<-read.table("dat.csv.kn_tp53_del5q_mono7",sep="\t",header=T,stringsAsFactors = F,comment.char = "")

# gene/varaint dict file 
all_gene<-read.table("dbNSFP4.0_gene.complete_id",sep="\t",header=F,stringsAsFactors = F,comment.char = "")
variant_gene_id_dict<-read.table("germ_somatic_variant_gene_dict.10ab_reg.cr_all",sep="\t",header=F,stringsAsFactors = F,comment.char = "")

# workflow control indicator
opc<-as.character(opt$output_file)
k_folds<-as.numeric(opt$input_kfolds)

# parameter setting
kk_x_list<-2:12
delta_list<-c(0.01,seq(0,1,by=0.1),0.99)   
marker_cutoff_metrics<-'2_2_3_1'                                                         
# variant count/percentage cutoff         
vc<-10
vpc<-0.9
marker_cutoff_metrics<-'2_2_3_1'

# km_name: 'gene' for gene-based genomic matrix / 'variant' for variant-based genomic matrix
km_name<-paste('germ_somatic_vcf_gene_',opc,sep="")

################################################################################################################################################################################ 
print("main workflow") 

# pre-regulation of genomic matrix 
genomic_matrix_reg<-distance_L1_GO_regulation_v(genomic_matrix, dim(genomic_matrix)[2]-1, -6)

# supervised-clustering of genomic matrix 
genomic_matrix_cluster<-supervised_clustering(genomic_matrix, dim(genomic_matrix_reg)[2]-1, k_folds, marker_cutoff_metrics, km_name)

# quality control of supervised clustering
genomic_matrix_cluster_qc<-score_qc_all(outdir,km_name,marker_cutoff_metrics, kk_x_list, delta_list,k_folds)

# save the data
save.image(file=paste(outdir,"/",'germ_somatic_vcf_gene_',opc,'.RData' ,sep=""))        
