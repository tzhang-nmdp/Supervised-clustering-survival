################################################################################################################################################################################ # 1.1 # function for supervised score_err function ( proportion / count) of multiple phenotypes
# The function for survival metrics of supervised clustering                                         
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
                options(warn = -1)                  
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

################################################################################################################################################################################  
# The main function for supervised score_err function proportion / count of multiple phenotypes ( determine the optimal result from the searching hyperparameter space)
                                              
score_qc_all <- function (outdir,km_name,marker_cutoff_metrics, kk_x_list, delta_list,k_folds)
    {                                              
    # Create matrix row k10Xdelta12=120  Col 3(c-index/logrank/siho-score)
    x<- length(kk_x_list)*length(delta_list)
    SVC_metics<-as.data.frame(matrix(NA, nrow=x,ncol=5))
    colnames(SVC_metics)<-c('c_index_med','log_rank_mc','siho_score_min','log_rank1_med','siho_score_med')
    
    SVC_metics2<-as.data.frame(matrix(NA, nrow=x,ncol=7))
    colnames(SVC_metics2)<-c('c_index','log_rank','siho_score','c_index2','log_rank1','log_rank2','siho_score2')
    # check metrics
    n<-1
    options(warn = -1)   
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
                print(melt(cv.fit[[i]][5])[,1])
                # loop for oob clusters
                for (c in unique(melt(cv.fit[[i]][5])[,1]))
                    {
                     svc_clin<- clin_data_tmp[clin_data_tmp$svc_id==c,]
                     svc_clin$oi_id<-0
                     print(dim(svc_clin))
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
            SVC_metics[n,]<-c(median(svc_cindex),length(svc_logrank[svc_logrank<=0.05])/length(svc_logrank),min(svc_siho),median(svc_logrank1),median(svc_siho))
            rownames(SVC_metics)[n]<-paste(k,delta,n,sep='_')
            SVC_metics2[n,]<-c(paste(svc_cindex,collapse=":"),paste(svc_logrank,collapse=":"),paste(svc_siho,collapse=":"),paste(svc_cindex2,collapse=":"),paste(svc_logrank1,collapse=":"),paste(svc_logrank2,collapse=":"),paste(svc_siho2,collapse=":"))
            rownames(SVC_metics2)[n]<-paste(k,delta,n,sep='_')  
            print(SVC_metics[n,])
            n=n+1        
            } # loop for delta
        } # loop for kk_x
    write.table(SVC_metics,paste(outdir,"/",km_name,"_summary_metrics_new.csv",sep=""),quote=F,sep='\t',col.names=T,row.names=T)
    write.table(SVC_metics2,paste(outdir,"/",km_name,"_summary_metrics_all_new.csv",sep=""),quote=F,sep='\t',col.names=T,row.names=T)   
    return(SVC_metics)
    }                   
               
################################################################################################################################################################################
# The sihoutte score function of clustering optimal metrics 
                                              
# get.function for mode
getmode <- function(v) 
    {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
    }

# clustering optimal metrics function
K_delta_calculator<-function (matrix_weight_random_split_inb, marker_index, oob, kk_x, delta)
        {
        max_silhouette_score<-0
        max_k_list<-c()
        silhouette_score_matrix<-data.frame('delta'=rep(NA,12))  

        # 2. iterate delta parameter  
        flush.console()         
        print(paste("starting iteration for delta:",delta,"at", Sys.time(), sep=" "))   
    
        # calculate silhouette_score of kmeans clustering and BIC for model based clustering
        set.seed(1)
        silhouette_score_k<-sil.score(matrix_weight_random_split_inb[,1:(marker_index-1)],nb.clus = kk_x, nb.run = 100, iter.max = 100,method = "euclidean")
        silhouette_score_k<-silhouette_score_k[kk_x] 
        BIC_score <- mclustBIC(matrix_weight_random_split_inb[,1:(marker_index-1)],G=kk_x) 
    
        return(list(silhouette_score_k,BIC_score))
        }

################################################################################################################################################################################ # 2.1 function for Kmeans silhouette score 
# The clustering sihoutte score function                                 
# Rep a value N times
rep <- function(val, n) {
  if(n < 1)
    stop('n must be greater than 0')
  sapply(1:n, function(x) val)
}


## Row means
rowmean <- function(X) {
  X <- as.matrix(X)
  as.vector(apply(X, 1, mean))
}


## Row sums
rowsums <- function(X) {
  X <- as.matrix(X)
  as.vector(apply(X, 1, sum))
}


## Vector maximums
maximum <- function(x, y) {
  if(length(x) != length(y))
    stop('dim mismatch')
  
  sapply(1:length(x), function(i) max(x[i], y[i]))
}

## Vector minimums
minimum <- function(x, y) {
  if(length(x) != length(y))
    stop('dim mismatch')
  
  sapply(1:length(x), function(i) min(x[i], y[i]))
}

## Silhouette score function 
silhouette_score <- function(X, labels, samples = FALSE) {
  if(!(is.data.frame(X) | is.matrix(X)))
    stop('X must be a dataframe or matrix')
  if(!all(unlist(lapply(X, is.numeric))))
    stop('X must be numeric')

  X <- as.matrix(X)
  y <- as.numeric(as.factor(labels))
  
  ## Check dims
  if(nrow(X) != length(y))
    stop('dim mismatch between X and y')
  
  # Get dist mat and labels
  distances <- as.matrix(dist(X, diag = T, upper = T))
  unique_labels <- unique(y)

  # If unique len is one, we have to stop
  if(length(unique_labels) < 2)
    stop('must be at least two unique labels') 
 
  # For sample i, store the mean distance of the cluster to which
  # it belongs in intra_clust_dists[i]
  intra_clust_dists <- rep(1, nrow(distances))
  
  # For sample i, store the mean distance of the second closest
  # cluster in inter_clust_dists[i]
  inter_clust_dists <- Inf * intra_clust_dists
  
  for(curr_label in unique_labels) {
    # Find inter_clust_dist for all samples belonging to the same
    # label (extract the rows into current_distances).
    mask = y == curr_label
    current_distances = distances[mask,]
    count_y = sum(as.numeric(mask)) ## How many of this label?   
    
    ## If count_y is 1, we have an issue
    if(count_y == 1)
      stop('label only appears once')

    # Leave out current sample.
    n_samples_curr_lab = sum(as.numeric(mask)) - 1
    if(n_samples_curr_lab != 0) {
      intra_clust_dists[mask] = rowsums(current_distances[,mask]) / n_samples_curr_lab
    }
    
    # Now iterate over all other labels, finding the mean
    # cluster distance that is closest to every sample.
    for(other_label in unique_labels) {
      if(other_label != curr_label) {
        other_mask = y == other_label
        other_distances = rowmean(current_distances[,other_mask])
        inter_clust_dists[mask] = minimum(inter_clust_dists[mask], other_distances)
      }
    }
  }
  
  sil_samples = inter_clust_dists - intra_clust_dists
  sil_samples = sil_samples / maximum(intra_clust_dists, inter_clust_dists)
  
  ## Return
  ifelse(samples, sil_samples, mean(sil_samples))
}


## EX:
sil_ex <- function() 
    {
  X <- t(matrix(1:15, ncol=5))
  y <- c(0,0,0,1,1)
  silhouette_score(X, y)
    }     

 
################################################################################################################################################################################ 
 # survival c-index function        
concordance_index<-function(event_times, predicted_scores){
  event_times = abs(event_times)
  event_observed = (event_times > 0)
  concsumout<-concordance_summary_statistics(event_times, predicted_scores, event_observed)
  num_correct=concsumout[1]
  num_tied=concsumout[2]
  num_pairs =concsumout[3]
return((num_correct+num_tied/2)/num_pairs)
}

concordance_summary_statistics<-function(event_times, predicted_event_times, event_observed){
  valid_comparison<-function(time_a, time_b, event_a, event_b){
if (time_a == time_b){
  # Ties are only informative if exactly one event happened
  return(event_a != event_b)
} else if( event_a & event_b){
  return(TRUE)
} else if ((event_a & (time_a < time_b))){
  return (TRUE)
  } else if ((event_b & (time_b < time_a))){
  return (TRUE)
 } else {

return(FALSE)
}
}

concordance_value<-function(time_a, time_b, pred_a, pred_b, event_a, event_b){
 if (pred_a == pred_b) {
    return(c(FALSE, TRUE))
  } else if( pred_a < pred_b){
  return (c(((time_a < time_b) | ((time_a == time_b) & (event_a & !event_b))), FALSE))
} else {
 return( c((time_a > time_b) | ((time_a == time_b) &  (!event_a & event_b)), FALSE))
}
}

num_pairs = 0
num_correct = 0
num_tied = 0

for (a in 1:(length(event_times)-1)){
  time_a = event_times[a]
  pred_a = predicted_event_times[a]
  event_a = event_observed[a]
# Don't want to double count
for (b in ((a + 1):length(event_times))){
  time_b = event_times[b]
pred_b = predicted_event_times[b]
event_b = event_observed[b]

if (valid_comparison(time_a, time_b, event_a, event_b)){
num_pairs = num_pairs+1
crct_ties = concordance_value(time_a, time_b, pred_a, pred_b, event_a, event_b)
crct<-crct_ties[1]
ties<-crct_ties[2]

if (crct){
  num_correct =num_correct+1
}
if (ties){
  num_tied =num_tied+1
}
}
}
}
return (c(num_correct, num_tied, num_pairs))
}
         
################################################################################################################################################################################ 
# the function to save heatmap  
         
save_pheatmap_pdf <- function(x, filename, width=20, height=8) 
   {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
   }  
                 
