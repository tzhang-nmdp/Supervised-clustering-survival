# these functions are under development (No available data for the fully testing)

# Example :

# supervised_clustering_val_all(germ_somatic_vcf_reg, marker_index, marker_cutoff_metrics, km_name) 
    
#########################################################################################################################################################################
# The validation/prediction function  
         
validation_predict <- function ( model, input_matrix, weight_vector, marker_index, marker_cutoff_metrics, kk_x, num_oob)  
    {
    # transform orignal gene/variant matrix into weighted matrix
    input_matrix_tf<-t(apply(input_matrix[,1:(marker_index-1)],1, function(x) x * weight_vector ))
                             
    # retreive training models                         
    km_model<-model[[1]]
    mc_model<-model[[2]]
                             
    # predicting cluster id of validation data                         
    km_cluster_id<-cluster_label(input_matrix_tf,km_model[["centers"]]) 
    mc_cluster_id<-predict(mc_model, input_matrix_tf)$classification  
    names(km_cluster_id)<-rownames(input_matrix_tf)
    names(mc_cluster_id)<-rownames(input_matrix_tf)              
    km_mc_cluster_id_matrix<-cbind(km_cluster_id,mc_cluster_id)                     
    colnames(km_mc_cluster_id_matrix)<-c('km_id','mc_id')
                             
    # calculate score_err of multiple phenotypes for each cluster of validation data                           
    km_mc_cluster_id_score <-score_err_all(km_mc_cluster_id_matrix, marker_cutoff_metrics, kk_x,num_oob)
                             
    return (list(km_mc_cluster_id_score,km_cluster_id,mc_cluster_id))
    }
    
################################################################################################################################################################################ 
# The main workflow function for validation of supervised clustering with delta parameter and weight transformation                                                             
supervised_clustering_val_all<- function(input_matrix, marker_index, marker_cutoff_metrics, km_name) 
    {
    flush.console()    
    print(paste("The validation starting at :", Sys.time(),sep=" ")) 
    
    # basic setting
    num_oob<-10
    oob<-1
    kk_x_list<-2:12
    delta_list<-c( 0.01,seq(0,1,by=0.1),0.99)
    
    # load training model settings for validation 
    # loop checking of cluster number k 
    for (kk_x in kk_x_list)     
        { 
                  
        # loop checking of weight factor delta
        for (delta in delta_list)          
            { 
            load(paste(outdir,"/",km_name_t,"_delta_",delta,"_K_",kk_x,"_oob_",oob,"marker_cutoff",marker_cutoff_metrics,"_val_model.RData",sep="")) 
            flush.console() 
            print(paste("# predict cluster id and calculate the score for multiple phenotypes", Sys.time(),sep=" "))  
            km_mc_cluster_id_score_matrix <-validation_predict(km_mc_model, input_matrix, stat_go_weight_vector, marker_index, marker_cutoff_metrics, kk_x, num_oob)   
            write.table(km_mc_cluster_id_score_matrix[[1]],paste(outdir,"/",km_name,"_delta_",delta,"_K_",kk_x,"_oob_",oob,"marker_cutoff",marker_cutoff_metrics,"_replication.sig_validation.csv",sep=""), sep="\t", quote=F,col.names=T,row.names=F) 
            #print(km_mc_cluster_id_score_matrix)

            km_cluster_id<-km_mc_cluster_id_score_matrix[[2]]
            mc_cluster_id<-km_mc_cluster_id_score_matrix[[3]]
            input_matrix$km_id<-km_cluster_id[match(rownames(input_matrix),names(km_cluster_id))]
            input_matrix$mc_id<-mc_cluster_id[match(rownames(input_matrix),names(mc_cluster_id))]
            save(km_mc_cluster_id_score_matrix,input_matrix,km_mc_model, stat_go_weight_vector, marker_cutoff_metrics, km_cluster_id, mc_cluster_id, km_name, kk_x, delta, file=paste(outdir,"/",km_name,"_delta_",delta,"_K_",kk_x,"_oob_",oob,"marker_cutoff",marker_cutoff_metrics,"_replication.sig_validation.RData",sep=""))

            # retrieve sample id and clin info for each cluster ( all clusters)     
            clin_data_tmp<<-clin_data[clin_data$formattedRID_LabCorpID %in% rownames(input_matrix),]
            clin_data_tmp$km_id<-input_matrix[match(clin_data_tmp$formattedRID_LabCorpID,rownames(input_matrix)),'km_id']
            clin_data_tmp$mc_id<-input_matrix[match(clin_data_tmp$formattedRID_LabCorpID,rownames(input_matrix)),'mc_id']

            # The validation table or plot heatmap for different clinical phenotype of Kmeans clustering 
            annotation_col<-as.data.frame(clin_data_tmp[,c('km_id','mdstype','ipssr','HCT.CI')])
            colnames(annotation_col)<-c('subgroup','mdstype','ipssr','HCT.CI')
            annotation_col<-annotation_col[order(annotation_col$subgroup),]
            options(repr.plot.width=12, repr.plot.height=6)
            p<-pheatmap(t(annotation_col),cluster_row=F,cluster_col=F,fontsize=20)                  
            save_pheatmap_pdf(p,paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob,"_km_clin_val.pdf",sep=""))

             if (length(unique(clin_data_tmp$km_id))>1)   
                {            
                # plot coxph / survival / logrank of Kmeans clustering
                 clin_data_tmp$km_id<-as.factor(clin_data_tmp$km_id)
                 cox <- coxph(Surv(clin_data_tmp$intxsurv, clin_data_tmp$dead) ~ ipssr + mdstype + HMA + CHEMO +km_id, data=clin_data_tmp)
                 p=ggforest(cox, data=clin_data_tmp,fontsize = 1)
                 options(repr.plot.width=8, repr.plot.height=8) 
                 tiff(paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob,"_km_clin_cox_val.tiff",sep=""),width = 2000, height =1500,res=300) 
                 print(p)
                 dev.off()   

                 # keplan-Meier survival curve of Kmeans clustering
                 surv<- survfit(Surv(clin_data_tmp$intxsurv, clin_data_tmp$dead) ~km_id, data=clin_data_tmp)
                 options(repr.plot.width=8, repr.plot.height=12)  
                 p1<-ggsurvplot(surv, conf.int = FALSE, data = clin_data_tmp, pval = T,pval.method = TRUE, risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(), 
                 font.main = c(20, "bold", "black"),font.x = c(20, "bold.italic", "black"),font.y = c(20, "bold.italic", "black"),font.tickslab = c(20, "bold", "black"), font.table=c(20, "bold", "black"))
                 tiff(paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob,"_km_clin_surv_val.tiff",sep=""),width = 2000, height =1500,res=300) 
                 print(p1)
                 dev.off()  
                 }

                 
             #The validation table or plot heatmap for different clinical phenotype of Mclust clustering
             annotation_col<-as.data.frame(clin_data_tmp[,c('mc_id','mdstype','ipssr','HCT.CI')])
             colnames(annotation_col)<-c('subgroup','mdstype','ipssr','HCT.CI')
             annotation_col<-annotation_col[order(annotation_col$subgroup),]
             options(repr.plot.width=12, repr.plot.height=6)
             p<-pheatmap(t(annotation_col),cluster_row=F,cluster_col=F,fontsize=20)                  
             save_pheatmap_pdf(p,paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob,"_mc_clin_val.pdf",sep=""))

             if (length(unique(clin_data_tmp$mc_id))>1)   
                {            
                 # plot coxph / survival / logrank of Mclust clustering
                 clin_data_tmp$mc_id<-as.factor(clin_data_tmp$mc_id)
                 cox <- coxph(Surv(clin_data_tmp$intxsurv, clin_data_tmp$dead) ~ ipssr + mdstype + HMA + CHEMO +mc_id, data=clin_data_tmp)
                 p=ggforest(cox, data=clin_data_tmp,fontsize = 1)
                 options(repr.plot.width=8, repr.plot.height=8) 
                 tiff(paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob,"_mc_clin_cox_val.tiff",sep=""),width = 2000, height =1500,res=300) 
                 print(p)
                 dev.off()   

                 # keplan-Meier survival curve of Mclust clustering
                 surv<- survfit(Surv(clin_data_tmp$intxsurv, clin_data_tmp$dead) ~mc_id, data=clin_data_tmp)
                 options(repr.plot.width=8, repr.plot.height=12)  
                 p1<-ggsurvplot(surv, conf.int = FALSE, data = clin_data_tmp, pval = T,pval.method = TRUE, risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
                 font.main = c(20, "bold", "black"),font.x = c(20, "bold.italic", "black"),font.y = c(20, "bold.italic", "black"),font.tickslab = c(20, "bold", "black"), font.table=c(20, "bold", "black"))
                 tiff(paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob,"_mc_clin_surv_val.tiff",sep=""),width = 2000, height =1500,res=300) 
                 print(p1)
                 dev.off()  
                }
             print(paste("The validation ending:", Sys.time(),sep=" "))   
            }
        }
     }                       
                  
