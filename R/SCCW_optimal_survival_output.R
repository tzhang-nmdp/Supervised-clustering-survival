# the function for optimal survival summary
optimal_clustering_setting<-function(genomic_matrix_cluster, genomic_matrix_cluster_qc, input_matrix_r, km_name, marker_index, marker_cutoff_metrics, num_oob)
    {
    # optimal setting grid_search
    optimal_set<-genomic_matrix_cluster_qc[genomic_matrix_cluster_qc$log_rank_mc==0 & genomic_matrix_cluster_qc$siho_score_min>0,]
    optimal_set<-optimal_set[order(optimal_set$c_index_med, decreasing=T),]
    kk_x_best<-unlist(str_split(row.names(optimal_set)[1],'_'))[1]
    delta_best<-unlist(str_split(row.names(optimal_set)[1],'_'))[2]
    print(paste(km_name,"has_delta_best:",delta_best,"and_K_best:",kk_x_best))

    # extract the clustsering information for all samples with optimal setting grid_search
    optimal_genomic_sc<-genomic_matrix_cluster[[as.numeric(kk_x_best)-2]]
    optimal_genomic_sc_cluster_id<-optimal_genomic_sc@km_cluster_inb_oob_rel[[which(delta_list==delta_best)]][[1]]
    clin_data_tmp1_com<-clin_data[clin_data$formattedRID_LabCorpID %in% names(optimal_genomic_sc_cluster_id),]
    clin_data_tmp1_com$km_id<-optimal_genomic_sc_cluster_id[match(clin_data_tmp1_com$formattedRID_LabCorpID,names(optimal_genomic_sc_cluster_id))]   
        
    # validation subset samples
    if (exists("input_matrix_r"))
        {print("input_matrix_r")
        load(paste(outdir,"/",km_name,"_delta_",delta_best,"_K_",kk_x_best,"_oob_1marker_cutoff",marker_cutoff_metrics,"_val_model.RData",sep="")) 
        km_mc_cluster_id_score_matrix <-validation_predict(km_mc_model, input_matrix_r, stat_go_weight_vector, marker_index, marker_cutoff_metrics, as.numeric(kk_x_best), num_oob)   
        km_cluster_id<-km_mc_cluster_id_score_matrix[[2]]
        mc_cluster_id<-km_mc_cluster_id_score_matrix[[3]]
        input_matrix_r$km_id<-km_cluster_id[match(rownames(input_matrix_r),names(km_cluster_id))]
        input_matrix_r$mc_id<-mc_cluster_id[match(rownames(input_matrix_r),names(mc_cluster_id))]
        clin_data_tmp_com<-clin_data[clin_data$formattedRID_LabCorpID %in% rownames(input_matrix_r),]
        clin_data_tmp_com$km_id<-input_matrix_r[match(clin_data_tmp_com$formattedRID_LabCorpID,rownames(input_matrix_r)),'km_id']

        # combine all samples    
        clin_data_tmp_all_com<-rbind(clin_data_tmp_com,clin_data_tmp1_com)
        } else {
        clin_data_tmp_all_com<-clin_data_tmp1_com
        }

    # plot heatmap for different clinical phenotype                    
    annotation_col<-as.data.frame(clin_data_tmp_all_com[,c('km_id','mdstype','ipssr','HCT.CI')])
    colnames(annotation_col)<-c('subgroup','mdstype','ipssr','HCT.CI')
    annotation_col<-annotation_col[order(annotation_col$subgroup),]
    options(repr.plot.width=12, repr.plot.height=6)
    p0<-pheatmap(t(annotation_col),cluster_row=F,cluster_col=F,fontsize=20)                  
    save_pheatmap_pdf(p0,paste(outdir,"/",km_name,"_kk_x_",kk_x_best,"_delta_",delta_best,"_oob_1_km_clin.pdf",sep=""))
 
    
    # coxph mutlivariate forest plot of Kmeans clustering
    clin_data_tmp_all_com$km_id<-as.factor( clin_data_tmp_all_com$km_id)
    clin_data_tmp_all_com$km_id<-relevel(clin_data_tmp_all_com$km_id,ref='1')
    cox <- coxph(Surv(clin_data_tmp_all_com$intxsurv, clin_data_tmp_all_com$dead) ~ ipssr + mdstype + HMA + CHEMO +km_id, data=clin_data_tmp_all_com)
    p1=ggforest(cox, data=clin_data_tmp_all_com,fontsize =0.8)
    options(repr.plot.width=8, repr.plot.height=6) 
    tiff(paste(outdir,"/",km_name,"_kk_x_",kk_x_best,"_delta_",delta_best,"_oob_1_km_clin_cox_val.tiff",sep=""),width = 2000, height =1500,res=300) 
    print(p1 +  theme_set(theme_grey(base_size = 15)) )
    dev.off()  
    
    # keplan-Meier survival curve of Kmeans clustering
    surv<- survminer::surv_fit(Surv(clin_data_tmp_all_com$intxsurv, clin_data_tmp_all_com$dead) ~km_id, data=clin_data_tmp_all_com)
    p2<-ggsurvplot(surv,data = clin_data_tmp_all_com, conf.int = FALSE,  pval = F,pval.method = TRUE, risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(), font.main = c(20, "bold", "black"),font.x = c(20, "bold.italic", "black"),font.y = c(20, "bold.italic", "black"),font.tickslab = c(20, "bold", "black"), font.table=c(20, "bold", "black"))
    options(repr.plot.width=8, repr.plot.height=7)  
    tiff(paste(outdir,"/",km_name,"_kk_x_",kk_x_best,"_delta_",delta_best,"_oob_1_km_clin_surv_val.tiff",sep=""),width = 2000, height =1500,res=300) 
    print(p2)
    dev.off()   
    
    return(list(clin_data_tmp_all_com, kk_x_best, delta_best))
    }
