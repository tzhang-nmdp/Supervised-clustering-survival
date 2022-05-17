################################################################################################################################################################################           
# 2.3 clustering function for cluster fitting and scoring of in-bag/inb and out-of-bag/oob samples 
                             
clustering_cv_fit<-function(matrix_random_split, marker_index, marker_cutoff_metrics, kk_x, delta, oob, num_inb, num_oob,km_name)
    {  
    flush.console()
    print(paste(oob,"folds cross-validation starting:", Sys.time(),sep=" "))            
    options(warn = -1)
    matrix_random_split_inb<-matrix_random_split[matrix_random_split$io!=as.character(oob),]
    matrix_random_split_oob<-matrix_random_split[matrix_random_split$io==as.character(oob),] 
    
    # 2.3.1 transformation of weighted matrix by statistical effect
    print(paste("Stat data weighting:", Sys.time(),sep=" "))       
    matrix_weight_go_tmp<-distance_weight_normalization(matrix_random_split_inb,marker_index)
    matrix_weight_stats_random_split_inb<-matrix_weight_go_tmp[[1]]
    ds_wt_matrix<-matrix_weight_go_tmp[[2]]
    
    # 2.3.2 transformation of weighted matrix by GO similarity    
    print(paste("Go data weighting:", Sys.time(),sep=" ")) 
    if (grepl('gene', km_name))
        {
        matrix_weight_go_tmp<-go_weight_normalization_g(matrix_random_split_inb,marker_index)            
    } else if (grepl('variant',km_name)) {
        matrix_weight_go_tmp<-go_weight_normalization_v(matrix_random_split_inb,marker_index)          
    } else {print ("### km_name error!!!")}
    matrix_weight_go_random_split_inb<-matrix_weight_go_tmp[[1]]   
    go_wt_matrix<-matrix_weight_go_tmp[[2]] 
    
    # 2.3.3 calculate clustering optimal metrics of kmeans and Mclust clustering for cluster number "k"  and weight factor "delta"
    matrix_weight_random_split_inb<-matrix_weight_final(matrix_weight_stats_random_split_inb, matrix_weight_go_random_split_inb, delta)      
    K_optimizer<-K_delta_calculator(matrix_weight_random_split_inb, marker_index, oob, kk_x, delta)
    silho_score<-K_optimizer[1]
    BIC_score_list<-K_optimizer[2][[1]]
    BIC_score_list[is.na(BIC_score_list)]<- min(BIC_score_list[!is.na(BIC_score_list)])-10000
    BIC_model_name<-c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
    BIC_model<-BIC_model_name[which(BIC_score_list==max(BIC_score_list))]
    print(paste("Combined data weighting:", Sys.time(),sep=" "))           

    # 2.3.4 transform out-of-bag samples
    stat_go_weight_vector<-apply(matrix_weight_stats_random_split_inb[,1:(marker_index-1)],2, max)
    matrix_weight_random_split_oob<-matrix_random_split_oob
    matrix_weight_random_split_oob[,1:(marker_index-1)]<-t(apply(matrix_random_split_oob[,1:(marker_index-1)],1, function(x) x * stat_go_weight_vector ))
                                            
    # 2.3.5 run clustering for in-bag samples     
    set.seed(1)
    cluster_s3_inb_km <-kmeans(matrix_weight_random_split_inb[,1:(marker_index-1)], kk_x, iter.max = 20, nstart = 5)
    cluster_s3_inb_mc <- Mclust(matrix_weight_random_split_inb[,1:(marker_index-1)], G=kk_x,modelNames=BIC_model)
    qq<-1
    while (class(cluster_s3_inb_mc)=='NULL')
        {
        cluster_s3_inb_mc <- Mclust(matrix_weight_random_split_inb[,1:(marker_index-1)], G=kk_x-qq)
        qq<-qq+1
        }
    cluster_s3_inb_km_id <-cluster_s3_inb_km$cluster
    cluster_s3_inb_mc_id <-cluster_s3_inb_mc$classification    
    matrix_weight_random_split_inb$km_id<-cluster_s3_inb_km_id
    matrix_weight_random_split_inb$mc_id<-cluster_s3_inb_mc_id    
    print(paste(oob,"folds cross-validation kmeans clustering done:", Sys.time(),sep=" "))     
  
    # 2.3.6 store training model parameter and weight confusing matrix
    km_mc_model<-list()
    km_mc_model[[1]]<-cluster_s3_inb_km
    km_mc_model[[2]]<-cluster_s3_inb_mc
                                                                 
    # 2.3.7 calculate and identify the best clusterscore_err for supervision 
    score_err_inb_tmp <-score_err_all(matrix_weight_random_split_inb, marker_cutoff_metrics, kk_x,num_inb)
    flush.console()
    print(paste(oob,"folds cross-validation inb scoring done:", Sys.time(),sep=" "))
    save(km_mc_model,K_optimizer,stat_go_weight_vector, matrix_weight_random_split_inb, matrix_weight_random_split_oob, marker_index, file=paste(outdir,"/",km_name,"_delta_",delta,"_K_",kk_x,"_oob_",oob,"marker_cutoff",marker_cutoff_metrics,"_val_model.RData",sep="")) 
                                            
    # 2.3.8 identify the cluster in out-of bag samples               
    cluster_s3_oob_km_id<-cluster_label(matrix_weight_random_split_oob[,1:(marker_index-1)],cluster_s3_inb_km[["centers"]]) 
    cluster_s3_oob_mc_id<-predict(cluster_s3_inb_mc, matrix_weight_random_split_oob[,1:(marker_index-1)])$classification   
    names(cluster_s3_oob_km_id)<-rownames(matrix_weight_random_split_oob)
    names(cluster_s3_oob_mc_id)<-rownames(matrix_weight_random_split_oob)
    matrix_weight_random_split_oob$km_id<-cluster_s3_oob_km_id
    matrix_weight_random_split_oob$mc_id<-cluster_s3_oob_mc_id
 
                                            
    # 2.3.9 calculate and identify the best score_err for cross-validation oob data 
    score_err_oob_tmp <-score_err_all(matrix_weight_random_split_oob, marker_cutoff_metrics, kk_x,num_oob)
    flush.console() 
    score_err_oob_tmp_tmp<-score_err_oob_tmp[complete.cases(score_err_oob_tmp),]
    score_err_oob_tmp[is.na(score_err_oob_tmp)]<-999   
    
    # 2.3.10 concensus check between inb and oob  
    score_err_inb_tmp$inb<-oob                                            
    score_err_oob_tmp$oob<-oob
    score_err_oob_tmp$km_val<-0
    score_err_oob_tmp$mc_val<-0

    # checking surival/non-survival significance of oob and inb scores of Kmeans clustering
    if (is.element(1,score_err_oob_tmp$km_sig) & is.element(1,score_err_inb_tmp$km_sig) & is.element(2,score_err_oob_tmp$km_sig) & is.element(2,score_err_inb_tmp$km_sig))
        {
        score_err_oob_tmp$km_val<-12 

    # checking surival significance of oob and inb scores of Kmeans clustering 
    } else if (is.element(2,score_err_oob_tmp$km_sig) & is.element(2,score_err_inb_tmp$km_sig))
        {
        score_err_oob_tmp$km_val<-2   
        
    # checking non-survival significance of oob and inb scores of Kmeans clustering        
    } else if (is.element(1,score_err_oob_tmp$km_sig) & is.element(1,score_err_inb_tmp$km_sig))
        {
        score_err_oob_tmp$km_val<-1 
        }  
    
    # checking surival/non-survival significance of oob and inb scores of Mclust clustering
    if (is.element(1,score_err_oob_tmp$mc_sig) & is.element(1,score_err_inb_tmp$mc_sig) & is.element(2,score_err_oob_tmp$mc_sig) & is.element(2,score_err_inb_tmp$mc_sig))
        {
        score_err_oob_tmp$mc_val<-12 

    # checking surival significance of oob and inb scores of Mclust clustering    
    } else if (is.element(2,score_err_oob_tmp$mc_sig) & is.element(2,score_err_inb_tmp$mc_sig))
        {
        score_err_oob_tmp$mc_val<-2  
        
    # checking non-survival significance of oob and inb scores of Mclust clustering        
    } else if (is.element(1,score_err_oob_tmp$mc_sig) & is.element(1,score_err_inb_tmp$mc_sig))
        {
        score_err_oob_tmp$mc_val<-1 
        }     

    print(paste(oob,"folds cross-validation oob scoring done:", Sys.time(),sep=" "))   
    
    # 2.3.11 save and return cross-validation score matrix       
#save(stat_go_weight_vector,matrix_weight_stats_random_split_inb, matrix_weight_go_random_split_inb, ds_wt_matrix,go_wt_matrix, K_optimizer, matrix_weight_random_split_inb, cluster_s3_inb_km,cluster_s3_inb_mc, score_err_inb_tmp, file=paste(outdir,"/",km_name,"_delta_",delta,"_K_",kk_x,"_oob_",oob,"marker_cutoff",marker_cutoff_metrics,"_tmp.RData",sep="")) 
    return(list(score_err_inb_tmp,score_err_oob_tmp,cluster_s3_inb_km[["centers"]],cluster_s3_inb_km_id,cluster_s3_oob_km_id,cluster_s3_inb_mc_id,cluster_s3_oob_mc_id,silho_score,BIC_score_list))
    }
