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
library(randomForestSRC)
library(survXgboost)
library(purrr)
#library(memisc)
library(mlr)
pacman::p_load_gh("IyarLin/survXgboost")
pacman::p_load("survival")
pacman::p_load("xgboost")

# load the packages
sources_path <- c("/Supervised-clustering-survival/R/")
file.sources = list.files(sources_path,pattern="*.R")
file.sources =file.sources[which(file.sources!='SCCW_supervised_clustering.R')]
#sapply(file.sources,source,.GlobalEnv)
map(paste0(sources_path, file.sources), source) 

################################################################################################################################################################################
# loading the data and parameter settings

option_list<-list(make_option(c("-i", "--input_file"),type="character",help="input file", default=NA,metavar="filename"), make_option(c("-d", "--output_dir"),type="character",help="output dir", default=NA,metavar="filename"), make_option(c("-o", "--output_file"),type="character",help="output file", default=NA,metavar="filename"), make_option(c("-k", "--input_kfolds"),type="character",help="input kfolds", default=NA,metavar="filename"))
opt_parser<-OptionParser(option_list=option_list)
opt=parse_args(opt_parser)

# load genomic data and extra information
input<-as.character(opt$input_file)
load(input)

# outdir setting
outdir<-as.character(opt$output_dir)
system(paste('mkdir ', outdir,sep=''))
system(paste('chmod 777 -R ', outdir,sep=''))

# workflow control indicator
opc<-as.character(opt$output_file)
k_folds<-as.numeric(opt$input_kfolds)

# parameter setting
kk_x_list<-3:9
delta_list<-c(0.01,seq(0.1,0.9,by=0.1),0.99)
marker_cutoff_metrics<-'2_2_3_1'
num_inb<-10
num_oob<-5
# variant count/percentage cutoff
vc<-10
vpc<-0.9

# km_name: 'gene' for gene-based genomic matrix / 'variant' for variant-based genomic matrix
km_name<-paste('germ_somatic_',opc,sep="")

################################################################################################################################################################################
print("main workflow")

# pre-regulation of genomic matrix
if (opc=="test_gene")
    {
    germ_somatic_vcf_0.000001_gene_clin_kn_tmp$intxsurv<-0                                          
    germ_somatic_vcf_0.000001_gene_clin_kn_tmp$dead<-0         
    germ_somatic_vcf_0.000001_gene_clin_kn_tmp[,c('intxsurv','dead')]<-clin_data_kn[match(rownames(germ_somatic_vcf_0.000001_gene_clin_kn_tmp),clin_data_kn$formattedRID_LabCorpID),c('intxsurv','dead')] 
    # remove low variation variables
    germ_somatic_vcf_0.000001_gene_clin_kn_tmp<-germ_somatic_vcf_0.000001_gene_clin_kn_tmp[,c(which(colSums(germ_somatic_vcf_0.000001_gene_clin_kn_tmp[,1:(dim(germ_somatic_vcf_0.000001_gene_clin_kn_tmp)[2]-2)])<=dim(germ_somatic_vcf_0.000001_gene_clin_kn_tmp)[1]*vpc & colSums(germ_somatic_vcf_0.000001_gene_clin_kn_tmp[,1:(dim(germ_somatic_vcf_0.000001_gene_clin_kn_tmp)[2]-2)])>=vc),dim(germ_somatic_vcf_0.000001_gene_clin_kn_tmp)[2]-1,dim(germ_somatic_vcf_0.000001_gene_clin_kn_tmp)[2])]
    # subset for testing
    germ_somatic_vcf_0.000001_gene_clin_kn_tmp<-germ_somatic_vcf_0.000001_gene_clin_kn_tmp[,c(1:2000,dim(germ_somatic_vcf_0.000001_gene_clin_kn_tmp)[2]-1,dim(germ_somatic_vcf_0.000001_gene_clin_kn_tmp)[2])]
    genomic_matrix_reg<-distance_L1_GO_regulation_v(germ_somatic_vcf_0.000001_gene_clin_kn_tmp, dim(germ_somatic_vcf_0.000001_gene_clin_kn_tmp)[2]-1, -6)
    } else if ((!(grepl("test",opc))) && (!(grepl("_L1",opc)))) {
    genomic_matrix_reg<-distance_regulation(input_genomic_matrix, dim(input_genomic_matrix)[2]-1, -6) 
    } else if ((!(grepl("test",opc))) && (grepl( "gene",opc))) {
    genomic_matrix_reg<-distance_L1_GO_regulation_g(input_genomic_matrix, dim(input_genomic_matrix)[2]-1, -6)   
    } else if ((!(grepl("test",opc))) && (grepl("variant",opc))) {
    genomic_matrix_reg<-distance_L1_GO_regulation_v(input_genomic_matrix, dim(input_genomic_matrix)[2]-1, -6)    
    }

# supervised-clustering of genomic matrix
genomic_matrix_cluster<-supervised_clustering(genomic_matrix_reg, dim(genomic_matrix_reg)[2]-1, k_folds, marker_cutoff_metrics, km_name)

# quality control of supervised clustering
genomic_matrix_cluster_qc<-score_qc_all(outdir,km_name,marker_cutoff_metrics, kk_x_list, delta_list,k_folds)

# optimal setting grid_search
optimal_cluster<-optimal_clustering_setting(genomic_matrix_cluster, genomic_matrix_cluster_qc, input_matrix_r, km_name, marker_index, marker_cutoff_metrics, num_oob)
clin_data_tmp_all_com<-optimal_cluster[[1]]
kk_x_best<-optimal_cluster[[2]]
delta_best<-optimal_cluster[[3]]
save(genomic_matrix_reg,genomic_matrix_cluster,genomic_matrix_cluster_qc,optimal_cluster,file=paste(outdir,"/",km_name,"_tmp.RData",sep=""))

# survival model summary
survival_summary_optimal_clustering<-survival_summary_model(clin_data_tmp_all_com[,c('km_id','mdstype','ipssr','HCT.CI','intxsurv', 'dead')],kk_x_best,delta_best, km_name)
write.table(survival_summary_optimal_clustering,paste(outdir,"/",km_name,"_kk_x_",kk_x_best,"_delta_",delta_best,"_km_survival_summary_model.csv",sep=""),sep="\t",col.names=T,row.names=T,quote=F)

# save the data
save.image(file=paste(outdir,"/",km_name,".RData",sep=""))
