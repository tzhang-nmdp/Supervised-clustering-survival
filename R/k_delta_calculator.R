# The clustering model optimal metrics function (mainly for mclust)
                                              
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
           
