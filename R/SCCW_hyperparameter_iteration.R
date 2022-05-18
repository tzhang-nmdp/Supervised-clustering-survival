################################################################################################################################################################################
#The main workflow function for supervised clustering with delta parameter and weight transformation 
                             
supervised_clustering<-function(input_matrix, marker_index, k_folds, marker_cutoff_metrics, km_name)
    {
    flush.console()
    print(paste("Data preparation starting:", Sys.time(),sep=" ")) 

    # set up s4 object for result summary    
    setClass("supervised_kmeans_model", slots = c(kmm_name="character", km_delta="list", km_mc_inb_score="list", km_mc_oob_score="list", km_inb_center="list", km_cluster_inb_oob_rel="list",  mc_cluster_inb_oob_rel="list", km_inb_cluster="list", mc_inb_cluster="list", km_oob_cluster="list", mc_oob_cluster="list"))   
    
    # hyperparameters: cluster number k / weight factor delta
    kk_x_list<-2:12
    delta_list<-c( 0.01,seq(0,1,by=0.1),0.99)
    num_inb<-10
    num_oob<-5
    cv_cutoff<-ceiling(k_folds/2)+1 
    
    # temporary metrics
    all_models<- list() 
    max_score_err<- 0.05
    score_err_list<-c() 
    sig_table<-list()

    
    # 1. random split 5 folds matrix  
    flush.console()    
    print(paste("Data random splitting:", Sys.time(),sep=" "))           
    matrix_random_split<-random_split(input_matrix,k_folds)
    #save(matrix_random_split, file = paste(outdir,"/",km_name,"_random_split_matrix.RData",sep="")) 
    print(paste("Data preparation done:", Sys.time(),sep=" "))  

    ## E step ##   
    # loop checking of cluster number k 
    for (kk_x in kk_x_list)     
        { 
        max_score_err_tmp<-c()
        sig_table_kk<-list()
        km_obj<-new("supervised_kmeans_model",kmm_name=km_name)
        i=1
                  
        # loop checking of weight factor delta
        for (delta in delta_list)          
            { 
            km_obj@km_cluster_inb_oob_rel[[i]]<-list()       
            km_obj@mc_cluster_inb_oob_rel[[i]]<-list()   
            km_obj@km_inb_center[[i]]<-list()
            km_obj@km_inb_cluster[[i]]<-list() 
            km_obj@km_oob_cluster[[i]]<-list()
            km_obj@mc_inb_cluster[[i]]<-list()
            km_obj@mc_oob_cluster[[i]]<-list()            
            flush.console()  
            
            print(paste("# 1. Five folds cross-validation loop at k cluster:",kk_x,"at delta:", delta, "at:", Sys.time(),sep=" "))       
            # parallel programming with multiple cores
            registerDoParallel(cores=k_folds)  
            cv.fit<-foreach(oob=1:k_folds) %dopar% clustering_cv_fit(matrix_random_split, marker_index, marker_cutoff_metrics, kk_x, delta, oob, num_inb, num_oob,km_name) 
            save(cv.fit, file = paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_cv_fit.RData",sep=""))           
         
            print(paste("# 2. save cross-validation score matrix into S4 object at k cluster:",kk_x,"at delta:", delta, "at:", Sys.time(),sep=" "))     
            km_obj@km_mc_inb_score[[i]]<-cv.fit[[1]][[1]]
            km_obj@km_mc_oob_score[[i]]<-cv.fit[[1]][[2]]            
            for ( m in 2:k_folds)
                {
                km_obj@km_mc_inb_score[[i]]<-rbind(km_obj@km_mc_inb_score[[i]],cv.fit[[m]][[1]])  
                km_obj@km_mc_oob_score[[i]]<-rbind(km_obj@km_mc_oob_score[[i]],cv.fit[[m]][[2]])                          
                }
            for ( m in 1:k_folds)
                {
                km_obj@km_inb_center[[i]][[m]]<-cv.fit[[m]][[3]]    
                km_obj@km_inb_cluster[[i]][[m]]<-cv.fit[[m]][[4]]   
                km_obj@km_oob_cluster[[i]][[m]]<-cv.fit[[m]][[5]]
                km_obj@mc_inb_cluster[[i]][[m]]<-cv.fit[[m]][[6]] 
                km_obj@mc_oob_cluster[[i]][[m]]<-cv.fit[[m]][[7]]   
                }         
            print(km_obj@km_mc_oob_score[[i]]) 
            
            print(paste("# 3. relabel of cluster id between k-fold cross-validation at k cluster:",kk_x,"at delta:", delta, "at:", Sys.time(),sep=" "))  
            # 4.1. get ref centroid for 4.2
            ref_center<-km_obj@km_inb_center[[i]][1]  
            km_obj@km_cluster_inb_oob_rel[[i]][[1]]<-c(km_obj@km_inb_cluster[[i]][[1]],km_obj@km_oob_cluster[[i]][[1]]) 
            
            for ( n in 2:k_folds)
                {
                cls_center<-as.data.frame(km_obj@km_inb_center[[i]][n])
                # 4.2. call cluster_relabel to compare cluster id   
                cluster_cor<-tryCatch({cluster_relabel(as.data.frame(ref_center),cls_center)
                         },error = function(e) { 
                    matrix(c(1:dim(cls_center)[1],1:dim(cls_center)[1]),ncol=2) })
                print(cluster_cor)
                # 4.3. relabel the cluster id in different oob data
                if (length(unique(cluster_cor[,2]))==1)
                    {
                    print('cluster id warning!!!')
                    #km_obj@km_mc_oob_score[[i]][(kk_x*(n-1)+1):(kk_x*n),4]<-0
                } else {               
                    km_obj@km_inb_cluster[[i]][[n]]<-sapply(km_obj@km_inb_cluster[[i]][[n]], function (x) cluster_cor[cluster_cor[,1]==x,2])   
                    km_obj@km_oob_cluster[[i]][[n]]<-sapply(km_obj@km_oob_cluster[[i]][[n]], function (x) cluster_cor[cluster_cor[,1]==x,2]) 
                    }   
                km_obj@km_cluster_inb_oob_rel[[i]][[n]]<-c(km_obj@km_inb_cluster[[i]][[n]],km_obj@km_oob_cluster[[i]][[n]])                                                            
                }
            for ( n in 1:k_folds)                                                        
                {km_obj@mc_cluster_inb_oob_rel[[i]][[n]]<-c(km_obj@mc_inb_cluster[[i]][[n]],km_obj@mc_oob_cluster[[i]][[n]])}
                                        
            print(paste("# 4. concensus of score between k-fold cross-validation at k cluster:",kk_x,"at delta:", delta, "at:", Sys.time(),sep=" ")) 
            # 5.1. get oob with sig cluster
            cluster_oob<- km_obj@km_mc_oob_score[[i]] 
            cluster_oob$kk_x<-kk_x
            cluster_oob$delta<-delta                                                        
            km_cluster_oob_sig<-cluster_oob[cluster_oob$km_val %in% c(1,2,12),'oob'] 
            mc_cluster_oob_sig<-cluster_oob[cluster_oob$mc_val %in% c(1,2,12),'oob']                                                         
            km_cluster_oob_sig2<-cluster_oob[cluster_oob$km_val %in% c(2,12),'oob'] 
            mc_cluster_oob_sig2<-cluster_oob[cluster_oob$mc_val %in% c(2,12),'oob']   
                                                            
            # 5.2. check if any matching cluster id or the mode of matching cluster id for each oob
            km_cluster_oob_sig_cv<-unique(km_cluster_oob_sig)                                                        
            km_cluster_oob_sig_cv_freq<-length(km_cluster_oob_sig_cv)
                                                            
            mc_cluster_oob_sig_cv<-unique(mc_cluster_oob_sig)                           
            mc_cluster_oob_sig_cv_freq<-length(mc_cluster_oob_sig_cv)
                                                        
            score_err_oob<-as.data.frame(km_obj@km_mc_oob_score[[i]])
            score_err_oob<-score_err_oob[complete.cases(score_err_oob),] 

            sig_table<-rbind(sig_table,cluster_oob)      
                                                        
            print(paste("# 6. checking and summarizing the clinically significant cluster at k cluster:",kk_x,"at delta:", delta, "at:", Sys.time(),sep=" "))   
            # checking non-survival significance of Kmeans clustering between oobs                                               
            if (km_cluster_oob_sig_cv_freq>= cv_cutoff)
                {                   
                pdf(paste(outdir,"/",km_name,"_deltaxxx_",delta,"_K_",kk_x,"marker_cutoff",marker_cutoff_metrics,"_tmp.oob_sig_val.pdf",sep=""), height=16, width=40) 
                grid.table(cluster_oob,rows = NULL)     
                dev.off()                
                for (oob_sig in unique(km_cluster_oob_sig))
                    {
                    # 6.0 retrieve sample id and clin info for each cluster ( all clusters)                     
                    cluster_inb_oob_rel<-km_obj@km_cluster_inb_oob_rel[[i]][[oob_sig]]
                    clin_data_tmp<<-clin_data[clin_data$formattedRID_LabCorpID %in% names(cluster_inb_oob_rel),]
                    clin_data_tmp$cl_id<-cluster_inb_oob_rel[match(clin_data_tmp$formattedRID_LabCorpID,names(cluster_inb_oob_rel))]
                    #print(clin_data_tmp[1:5,1:10])
                    print(clin_data_tmp$cl_id)   
                    
                    # 6.1 table or plot heatmap for different clinical phenotype                    
                    annotation_col<-as.data.frame(clin_data_tmp[,c('cl_id','mdstype','ipssr','HCT.CI')])
                    colnames(annotation_col)<-c('subgroup','mdstype','ipssr','HCT.CI')
                    annotation_col<-annotation_col[order(annotation_col$subgroup),]
                    options(repr.plot.width=12, repr.plot.height=6)
                    p<-pheatmap(t(annotation_col),cluster_row=F,cluster_col=F,fontsize=20)                  
                    save_pheatmap_pdf(p,paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob_sig,"_km_clin.pdf",sep=""))
               
                    # 6.2 plot coxph / survival / logrank
                    if (length(unique(clin_data_tmp$cl_id))>1 & oob_sig %in% km_cluster_oob_sig2)
                        {
                        print(paste("###### significant KM cluster for ",km_name, "k at:",kk_x,"and delta at:",delta,sep=" "))                         
                        clin_data_tmp$cl_id<-as.factor(clin_data_tmp$cl_id)
                        cox <- coxph(Surv(clin_data_tmp$intxsurv, clin_data_tmp$dead) ~ ipssr + mdstype + HMA + CHEMO +cl_id, data=clin_data_tmp)
                        p=ggforest(cox, data=clin_data_tmp,fontsize = 1)
                        options(repr.plot.width=8, repr.plot.height=8) 
                        tiff(paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob_sig,"_km_clin_cox.tiff",sep=""),width = 2000, height =1500,res=300) 
                        print(p)
                        dev.off()   
                        
                        # 6.3 keplan-Meier survival curve
                        surv<- survfit(Surv(clin_data_tmp$intxsurv, clin_data_tmp$dead) ~cl_id, data=clin_data_tmp)
                        options(repr.plot.width=8, repr.plot.height=12)  
                        p1<-ggsurvplot(surv, conf.int = FALSE, data = clin_data_tmp, pval = T,pval.method = TRUE, risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(), 
                        font.main = c(20, "bold", "black"),font.x = c(20, "bold.italic", "black"),font.y = c(20, "bold.italic", "black"),font.tickslab = c(20, "bold", "black"), font.table=c(20, "bold", "black"))
                        tiff(paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob_sig,"_km_clin_surv.tiff",sep=""),width = 2000, height =1500,res=300) 
                        print(p1)
                        dev.off()                          
                        }   
                    }
                }
# checking non-survival significance of Model based clustering between oobs                                            
            if (mc_cluster_oob_sig_cv_freq>= cv_cutoff) {                   
            pdf(paste(outdir,"/",km_name,"_deltaxxx_",delta,"_K_",kk_x,"marker_cutoff",marker_cutoff_metrics,"_tmp.oob_sig_val.pdf",sep=""), height=16, width=40) 
            grid.table(cluster_oob,rows = NULL)     
            dev.off()                   
            for (oob_sig in unique(mc_cluster_oob_sig))
                {
                # 6.0 retrieve sample id and clin info for each cluster ( all clusters)                           
                cluster_inb_oob_rel<-km_obj@mc_cluster_inb_oob_rel[[i]][[oob_sig]]
                clin_data_tmp<<-clin_data[clin_data$formattedRID_LabCorpID %in% names(cluster_inb_oob_rel),]
                clin_data_tmp$cl_id<-cluster_inb_oob_rel[match(clin_data_tmp$formattedRID_LabCorpID,names(cluster_inb_oob_rel))]
                #print(clin_data_tmp[1:5,1:10])
                print(clin_data_tmp$cl_id)   

                # 6.1 table or plot heatmap for different clinical phenotype                        
                annotation_col<-as.data.frame(clin_data_tmp[,c('cl_id','mdstype','ipssr','HCT.CI')])
                colnames(annotation_col)<-c('subgroup','mdstype','ipssr','HCT.CI')
                annotation_col<-annotation_col[order(annotation_col$subgroup),]
                options(repr.plot.width=12, repr.plot.height=6)
                p<-pheatmap(t(annotation_col),cluster_row=F,cluster_col=F,fontsize=20)                  
                save_pheatmap_pdf(p,paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob_sig,"_mc_clin.pdf",sep=""))

                # 6.2 plot coxph / survival / logrank
                if (length(unique(clin_data_tmp$cl_id))>1 & oob_sig %in% mc_cluster_oob_sig2)
                    {
                    print(paste("###### significant GMM cluster for ",km_name, "k at:",kk_x,"and delta at:",delta,sep=" "))                        
                    clin_data_tmp$cl_id<-as.factor(clin_data_tmp$cl_id)
                    cox <- coxph(Surv(clin_data_tmp$intxsurv, clin_data_tmp$dead) ~ ipssr + mdstype + HMA + CHEMO +cl_id, data=clin_data_tmp)
                    p=ggforest(cox, data=clin_data_tmp,fontsize = 1)
                    options(repr.plot.width=8, repr.plot.height=8) 
                    tiff(paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob_sig,"_mc_clin_cox.tiff",sep=""),width = 2000, height =1500,res=300) 
                    print(p)
                    dev.off()   

                    # 6.3 keplan-Meier survival curve
                    surv<- survfit(Surv(clin_data_tmp$intxsurv, clin_data_tmp$dead) ~cl_id, data=clin_data_tmp)
                    options(repr.plot.width=8, repr.plot.height=12)  
                    p1<-ggsurvplot(surv, conf.int = FALSE, data = clin_data_tmp, pval = T,pval.method = TRUE, risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),
                    font.main = c(20, "bold", "black"),font.x = c(20, "bold.italic", "black"),font.y = c(20, "bold.italic", "black"),font.tickslab = c(20, "bold", "black"), font.table=c(20, "bold", "black"))
                    tiff(paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_oob_",oob_sig,"_mc_clin_surv.tiff",sep=""),width = 2000, height =1500,res=300) 
                    print(p1)
                    dev.off()  
                    }    
               }
            }                                                            

            if (km_cluster_oob_sig_cv_freq < cv_cutoff & mc_cluster_oob_sig_cv_freq < cv_cutoff) {
                flush.console()                
                print(paste("No significant cluster for ",km_name, "k at:",kk_x,"and delta at:",delta,sep=" "))
                } 
                                                            
            # store significant score table between oobs                                                
            sig_table_kk<-rbind(sig_table_kk,cluster_oob)                                                  
            max_score_err_tmp[i]<-min(score_err_oob[score_err_oob$km_cluster_num>num_oob, 'km_surv_cox'],score_err_oob[score_err_oob$mc_cluster_num>num_oob, 'mc_surv_cox'])  
                                                            
            i<-i+1                                                               
            } # for delta loop
            
            # store significant score table between oobs                                                              
            write.table(sig_table_kk,paste(outdir,"/",km_name,"kk_x",kk_x,"marker_cutoff",marker_cutoff_metrics,"_tmp.oob_sig_val_all.csv",sep=""), sep="\t", quote=F,col.names=T,row.names=F)       
            all_models[[kk_x-1]]<-km_obj 
                 
            ## M step ##
            # mark the hyperparameters: cluster number k / weight factor delta with the best score                                                
            if (min(max_score_err_tmp)<max_score_err)
                {
                max_score_err<-min(max_score_err_tmp)
                kk_x_max<-kk_x
                delta_max<-delta_list[which(max_score_err_tmp==min(max_score_err_tmp))]
                flush.console()                  
                } 
                                                            
        }  # for kk_x loop 
        
        # store training model settings for validation                                                    
        write.table(sig_table,paste(outdir,"/",km_name,"marker_cutoff",marker_cutoff_metrics,"_tmp.oob_sig_val_all.csv",sep=""), sep="\t", quote=F,col.names=T,row.names=F)  
        print(paste("The hyperparameter with max kk_x at:", kk_x_max, "and max delta at:",delta_max, sep=" "))                                                   
        print(paste("The hyperparameter search ending:", Sys.time(),sep=" "))    
        print("#########################################################################################################################")        
        return(all_models)
    }     
           
