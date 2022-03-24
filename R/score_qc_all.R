# The main function for supervised score_err function proportion / count of multiple phenotypes ( determine the optimal result from the searching hyperparameter space)
                                              
score_qc_all <- function (outdir,km_name,marker_cutoff_metrics, kk_x_list, delta_list,k_folds)
    {                                              
    # Create matrix row k10Xdelta12=120  Col 3(c-index/logrank/siho-score)
    x<- length(kk_x_list)*length(delta_list)
    SVC_metics<-as.data.frame(matrix(NA, nrow=x,ncol=18))
    colnames(SVC_metics)<-c('c_index_max','c_index2_max','c_index_med','c_index2_med','log_rank_m','log_rank_mc','log_rank1_mc','log_rank_med','log_rank1_med','log_rank2_med','log_rank_min','log_rank2_mean','siho_score_med','siho_score2_med','siho_score_mc','siho_score2_mc','siho_score_min','siho_score2_min')

    SVC_metics2<-as.data.frame(matrix(NA, nrow=x,ncol=7))
    colnames(SVC_metics2)<-c('c_index','c_index2','log_rank','log_rank1','log_rank2','siho_score','siho_score2')
    # check metrics
    n<-1
    for (k in kk_x_list)
        {
        for (delta in delta_list)
            { 
            load(paste(outdir,"/",km_name,"_kk_x_",k,"_delta_",delta,"_cv_fit.RData",sep="")) 
            svc_cindex<-c()
            svc_cindex2<-c()        
            svc_siho<-c()
            svc_siho2<-c()  
            svc_logrank<-c()
            svc_logrank1<-c()            
            svc_logrank2<-c()  
            
            # loop for k_folds            
            for (i in 1:k_folds)
                { 
                flush.console() 
                print(Sys.time())            
                print(paste(k,delta,i,sep=":"))
                if (grepl('kn',opc))
                    {
                    clin_data_tmp<-clin_data_kn
                } else {
                    clin_data_tmp<-clin_data                    
                    } 
                svc_id<-rbind(melt(cv.fit[[i]][4]),melt(cv.fit[[i]][5]))
                svc_id1<-melt(cv.fit[[i]][4])                  
                svc_id2<-melt(cv.fit[[i]][5])                  
                clin_data_tmp$svc_id<-svc_id[match(clin_data_tmp$formattedRID_LabCorpID, rownames(svc_id)),1]
                
                # survival cindex and logranl
                svc_cox <- coxph(Surv(clin_data_tmp$intxsurv, clin_data_tmp$dead) ~ ipssr + mdstype + HMA + CHEMO + svc_id, data=clin_data_tmp)
                svc_cindex<-c(svc_cindex,as.numeric(svc_cox$concordance[6]))
                clin_data_tmp1<-clin_data_tmp[clin_data_tmp$formattedRID_LabCorpID %in% rownames(svc_id1),]  
                clin_data_tmp1$svc_id<-svc_id1[match(clin_data_tmp1$formattedRID_LabCorpID, rownames(svc_id1)),1]                  
                clin_data_tmp2<-clin_data_tmp[clin_data_tmp$formattedRID_LabCorpID %in% rownames(svc_id2),]
                clin_data_tmp2$svc_id<-svc_id2[match(clin_data_tmp2$formattedRID_LabCorpID, rownames(svc_id2)),1]   
         
                svc_survdiff1 <-survdiff(Surv(clin_data_tmp1$intxsurv, clin_data_tmp1$dead) ~ clin_data_tmp1$svc_id)
                print('1')
                svc_logrank1<-c(svc_logrank1, pchisq(svc_survdiff1$chisq, length(svc_survdiff1$n)-1, lower.tail = FALSE) )   
                if (length(unique(clin_data_tmp2$svc_id))>1)
                    {
                    svc_cox2<- coxph(Surv(clin_data_tmp2$intxsurv, clin_data_tmp2$dead) ~ ipssr + mdstype + HMA + CHEMO + svc_id, data=clin_data_tmp2)   
                    svc_cindex2<-c(svc_cindex2,as.numeric(svc_cox2$concordance[6]))                                      
                    svc_survdiff2 <-survdiff(Surv(clin_data_tmp2$intxsurv, clin_data_tmp2$dead) ~ clin_data_tmp2$svc_id)
                    print('2')
                    svc_logrank2<-c(svc_logrank2, pchisq(svc_survdiff2$chisq, length(svc_survdiff2$n)-1, lower.tail = FALSE) )        
                }  else {
                    svc_cindex2<-svc_cindex
                     svc_survdiff2 <-survdiff(Surv(clin_data_tmp$intxsurv, clin_data_tmp$dead) ~ clin_data_tmp$svc_id)              
                    svc_logrank2<-c(svc_logrank2, pchisq(svc_survdiff2$chisq, length(svc_survdiff2$n)-1, lower.tail = FALSE) )       
                     print('3')                   
                    }            
                
                # loop for oob clusters
                for (c in unique(melt(cv.fit[[i]][5])[,1]))
                    {
                     svc_clin<- clin_data_tmp[clin_data_tmp$svc_id==c,]
                     svc_clin$oi_id<-0
                    if (length(svc_clin[svc_clin$formattedRID_LabCorpID %in% rownames(melt(cv.fit[[i]][5])), 'oi_id'])>=10)  
                        {
                         svc_clin[svc_clin$formattedRID_LabCorpID %in% rownames(melt(cv.fit[[i]][5])), 'oi_id']<-1
                         svc_survdiff <-survdiff(Surv(svc_clin$intxsurv, svc_clin$dead) ~ svc_clin$oi_id)
                         print('4')
                         svc_logrank<-c(svc_logrank, pchisq(svc_survdiff$chisq, length(svc_survdiff$n)-1, lower.tail = FALSE) )                    
                    } else {
                         svc_clin_sub<-svc_clin[svc_clin$formattedRID_LabCorpID %in% rownames(melt(cv.fit[[i]][5])), ]
                         svc_clin[svc_clin$formattedRID_LabCorpID %in% rownames(melt(cv.fit[[i]][5])), 'oi_id']<-1  
                         svc_clin<-rbind(svc_clin,svc_clin_sub)
                         svc_survdiff <-survdiff(Surv(svc_clin$intxsurv, svc_clin$dead) ~ svc_clin$oi_id)
                         print('5 ')
                         svc_logrank<-c(svc_logrank, pchisq(svc_survdiff$chisq, length(svc_survdiff$n)-1, lower.tail = FALSE) )                                        
                         }
                    } # loop for c cluster of oob

                # check silhouette score
                load(paste(outdir,"/",km_name,"_delta_",delta,"_K_",k,"_oob_",i,"marker_cutoff",marker_cutoff_metrics,"_val_model.RData",sep=""))            
                options(warn = -1)
                set.seed(1)            
                silhouette_score_k<-sil.score(matrix_weight_random_split_inb[,1:(marker_index-1)],nb.clus = k, nb.run = 10, iter.max = 50,method = "euclidean")
                dis = dist(matrix_weight_random_split_inb[,1:(marker_index-1)])^2
                res = kmeans(matrix_weight_random_split_inb[,1:(marker_index-1)],k,iter.max = 50,nstart = 10)
                silhouette_score_k2<-silhouette(res$cluster, dis) 
                svc_siho<-c(svc_siho,silhouette_score_k[k])
                svc_siho2<-c(svc_siho2,mean(silhouette_score_k2[,3]))            
                }  # loop for i inb/oob 
            if (length(svc_logrank[svc_logrank<=0.05])>=round(length(svc_logrank)/2))
                {
                svc_logrank_m<- min(svc_logrank)
            }  else {
                svc_logrank_m<- median(svc_logrank) 
                }        
            sil_cls<-aggregate(silhouette_score_k2[,3],by=list(silhouette_score_k2[,1]),mean)[,2]
            if ((min(svc_siho)<0) && (median(svc_siho)<0) && (min(sil_cls)<0) && (median(sil_cls)<0))
                {
                print('clustering quality questionable!!!')
                }
            SVC_metics[n,]<-c(max(svc_cindex),max(svc_cindex2),median(svc_cindex),median(svc_cindex2),svc_logrank_m,length(svc_logrank[svc_logrank<=0.05])/length(svc_logrank),length(svc_logrank1[svc_logrank1<=0.05])/length(svc_logrank1),median(svc_logrank),median(svc_logrank1),median(svc_logrank2),min(svc_logrank),mean(svc_logrank2),median(svc_siho),median(svc_siho2),length(svc_logrank[svc_siho<0])/length(svc_siho),length(svc_logrank[svc_siho2<0])/length(svc_siho2),min(svc_siho),min(svc_siho2))
            rownames(SVC_metics)[n]<-paste(k,delta,n,sep='_')
            SVC_metics2[n,]<-c(paste(svc_cindex,collapse=":"),paste(svc_cindex2,collapse=":"),paste(svc_logrank,collapse=":"),paste(svc_logrank1,collapse=":"),paste(svc_logrank2,collapse=":"),paste(svc_siho,collapse=":"),paste(svc_siho2,collapse=":"))
            rownames(SVC_metics2)[n]<-paste(k,delta,n,sep='_')  
            print(SVC_metics[n,])
            n=n+1        
            } # loop for delta
        } # loop for kk_x
    write.table(SVC_metics,paste(outdir,"/",km_name,"_summary_metrics_new.csv",sep=""),quote=F,sep='\t',col.names=T,row.names=T)
    write.table(SVC_metics2,paste(outdir,"/",km_name,"_summary_metrics_all_new.csv",sep=""),quote=F,sep='\t',col.names=T,row.names=T)   
    return(SVC_metics)
    }             
           
