# The function for supervised score_err function proportion / count of multiple phenotypes ( mainly for inb/oob validation)
                                              
score_err_all <- function (cluster_matrix, marker_cutoff_metrics, kk_x, num_cutoff)
    {
    marker_cutoff_metrics_tmp<-unlist(str_split(marker_cutoff_metrics,"_"))
    cluster_marker_matrix_list<-list()
    i=1
    cluster_id<-1:kk_x
    
    # loop checking of the clustering from Kmeans and Mculst
    for (cls in c('km','mc'))
        {
        id_idx<-which(colnames(cluster_matrix)==paste(cls,'id',sep="_"))
        cluster_id_tmp<-unique(cluster_matrix[,id_idx])
        cluster_marker_matrix<-data.frame('cluster'=cluster_id[order(cluster_id,decreasing = F)],'cluster_num'=NA,'cluster_marker'=NA,'mdstype'=NA,'ipssr'=NA,'hct_ci'=NA,'surv'=NA,'surv_cox'=NA,'sig'=NA)
        j=1
        
        # checking in case of single cluster from oob / validation subset
        if (length(cluster_id_tmp)>1)
            {
            
            # loop checking for each cluster
            for ( cls_id in cluster_id_tmp[order(cluster_id_tmp,decreasing = F)])
                {
                
                # retrieve sample id and clin info for each cluster ( by 1 cluster vs other ) 
                cluster_sample<-rownames(cluster_matrix)[cluster_matrix[,id_idx]==cls_id]
                oth_sample<-rownames(cluster_matrix)[cluster_matrix[,id_idx]!=cls_id]
                clin_data_tmp<-clin_data
                clin_data_tmp[is.na(clin_data_tmp)]<-0
                clin_data_tmp[clin_data_tmp$ipssr==4,'ipssr']<-0
                rownames(clin_data_tmp)<-clin_data_tmp$formattedRID_LabCorpID
                cluster_clin<-clin_data[rownames(clin_data_tmp) %in% cluster_sample, c('mdstype','ipssr','HCT.CI','dead','intxsurv','HMA','CHEMO') ]
                cluster_clin$cl_id<-1
                oth_clin<-clin_data[rownames(clin_data_tmp) %in% oth_sample,  c('mdstype','ipssr','HCT.CI','dead','intxsurv','HMA','CHEMO')]
                oth_clin$cl_id<-0
                all_clin<-rbind(cluster_clin,oth_clin)
                clin_p<-c()
                
                # loop checking of multiple phenotypes for each cluster
                for (n in 1:5)
                    {    
                    # calculate the values of clin info by cluster       
                    marker_cutoff<-marker_cutoff_metrics_tmp[n] 
                    cluster_num<-dim(cluster_clin)[1]
                    oth_no<-dim(oth_clin)[1]
                    cluster_cutoff_num<-length(cluster_clin[cluster_clin[,n]>= marker_cutoff,n])
                    oth_cutoff_num<-length(oth_clin[oth_clin[,n]>= marker_cutoff,n])  
                    cluster_marker<-cluster_cutoff_num/cluster_num   
                    
                    # statistical scoring of multiple phenotypes for each cluster            
                    if (n<4)
                        {            
                        clin_p[n]<-chisq.test(matrix(c(cluster_cutoff_num,oth_cutoff_num,(cluster_num-cluster_cutoff_num),(oth_no-oth_cutoff_num)),ncol=2))$p.value                 
                    } else {
                        surv_diff<-survdiff(Surv(all_clin$intxsurv, all_clin$dead) ~ all_clin$cl_id)
                        clin_p[4]<-pchisq(surv_diff$chisq, length(surv_diff$n)-1, lower.tail = FALSE)  
                        all_clin$cl_id<-as.factor(all_clin$cl_id)
                        cox <- coxph(Surv(all_clin$intxsurv, all_clin$dead) ~ ipssr + mdstype + HMA + CHEMO +cl_id, data=all_clin)
                        clin_p[5]<-summary(cox)$coef[5,5]                        
                        }
                    } 
                
                # add the significant tag for each cluster
                if ((clin_p[4]<=0.05 | clin_p[5]<=0.05) & cluster_num>=num_cutoff)
                    {cluster_marker_matrix[j,]<-c(j,cluster_num,cluster_marker, clin_p,2)
                } else if (min(clin_p[1:3])<=0.05 & cluster_num>=num_cutoff) {
                    cluster_marker_matrix[j,]<-c(j,cluster_num,cluster_marker, clin_p,1)
                } else {
                    cluster_marker_matrix[j,]<-c(j,cluster_num,cluster_marker, clin_p,0)
                    }            
                j=j+1
                }
            
            } else {
            flush.console()            
            print('cluster id warning!!!')
            }
        
        colnames(cluster_marker_matrix)<-paste(cls,colnames(cluster_marker_matrix),sep="_")   
        cluster_marker_matrix[is.na(cluster_marker_matrix)]<-999
        cluster_marker_matrix_list[[i]]<-cluster_marker_matrix
        i=i+1
        }
    cluster_marker_matrix_list<-cbind(cluster_marker_matrix_list[[1]],cluster_marker_matrix_list[[2]])
    
    return (cluster_marker_matrix_list) 
    }    
         
