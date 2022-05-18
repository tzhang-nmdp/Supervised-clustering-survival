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

# load the packages
sources_path <- c("/Supervised-clustering-survival/R")
file.sources = list.files(sources_path,pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

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
clin_data<-read.table("/Supervised-clustering-survival/Example/dat.csv.cr",sep="\t",header=T,stringsAsFactors = F,comment.char = "")
clin_data_kn<-read.table("/Supervised-clustering-survival/Example/dat.csv.kn_tp53_del5q_mono7",sep="\t",header=T,stringsAsFactors = F,comment.char = "")

# gene/varaint dict file
all_gene<-read.table("/Supervised-clustering-survival/Example/dbNSFP4.0_gene.complete_id",sep="\t",header=F,stringsAsFactors = F,comment.char = "")
variant_gene_id_dict<-read.table("/Supervised-clustering-survival/Example/germ_somatic_variant_gene_dict.10ab_reg.cr_all.cr",sep="\t",header=F,stringsAsFactors = F,comment.char = "")

# workflow control indicator
opc<-as.character(opt$output_file)
k_folds<-as.numeric(opt$input_kfolds)

# parameter setting
kk_x_list<-3:9
delta_list<-c(0.01,seq(0,1,by=0.1))
marker_cutoff_metrics<-'2_2_3_1'
# variant count/percentage cutoff
vc<-10
vpc<-0.9

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

# optimal setting grid_search
optimal_set<-genomic_matrix_cluster_qc[genomic_matrix_cluster_qc$log_rank_mc==0 & genomic_matrix_cluster_qc$siho_score_min>0,]
optimal_set<-optimal_set[order(optimal_set$c_index_med, decreasing=T),]
kk_x_best<-unlist(str_split(row.names(optimal_set)[1],'_'))[1]
delta_best<-unlist(str_split(row.names(optimal_set)[1],'_'))[2]

# extract the clustsering information for all samples with optimal setting grid_search
optimal_genomic_matrix<-genomic_matrix_cluster[[kk_x_best-2]][[which(delta_list==delta_best)]]
load(paste(outdir,"/",km_name,"_delta_",delta_best,"_K_",kk_x_best,"_oob_",oob,"marker_cutoff",marker_cutoff_metrics,"_val_model.RData",sep="")) 
km_mc_cluster_id_score_matrix <-validation_predict(km_mc_model, input_matrix_r, stat_go_weight_vector, marker_index, marker_cutoff_metrics, kk_x, num_oob)   
km_cluster_id<-km_mc_cluster_id_score_matrix[[2]]
mc_cluster_id<-km_mc_cluster_id_score_matrix[[3]]
input_matrix_r$km_id<-km_cluster_id[match(rownames(input_matrix_r),names(km_cluster_id))]
input_matrix_r$mc_id<-mc_cluster_id[match(rownames(input_matrix_r),names(mc_cluster_id))]
clin_data_tmp_com<-clin_data[clin_data$formattedRID_LabCorpID %in% rownames(input_matrix_r),]
clin_data_tmp_com$km_id<-input_matrix_r[match(clin_data_tmp_com$formattedRID_LabCorpID,rownames(input_matrix_r)),'km_id']
clin_data_tmp1_com<-clin_data[clin_data$formattedRID_LabCorpID %in% names(cluster_inb_oob_rel),]
clin_data_tmp1_com$km_id<-cluster_inb_oob_rel[match(clin_data_tmp1_com$formattedRID_LabCorpID,names(cluster_inb_oob_rel))]   
clin_data_tmp_all_com<-rbind(clin_data_tmp_com,clin_data_tmp1_com)

# survival model summary
survival_summary_optimal_clustering<-survival_summary_model(clin_data_tmp_all_com, dim(optimal_genomic_matrix)[2]-1, km_name)

# save the data
save.image(file=paste(outdir,"/",'germ_somatic_vcf_gene_',opc,'.RData' ,sep=""))
