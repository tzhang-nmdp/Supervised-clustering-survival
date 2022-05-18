################################################################################################################################################################################  
# The function for Kmeans centroid computation
                                              
cluster_centroid <- function(input_matrix, marker_index, label) 
    {
    cluster_id<-unique(label)
    cluster_id<-cluster_id[order(cluster_id)]    
    centers<-as.data.frame(matrix(nrow=length(unique(label)),ncol=marker_index-1))
    colnames(centers)<-colnames(input_matrix)[1:(marker_index-1)]
    rownames(centers)<-cluster_id
    
    # compute variable mean for each cluster center
    for (i in cluster_id)
        {
        cluster_sub<-label[which(label==i)]
        cluster_sub_matrix<-input_matrix[rownames(input_matrix) %in% names(cluster_sub),1:(marker_index-1)]
        centers[i,]<-apply(cluster_sub_matrix,2,mean)
        }
    
    return (centers)     
    }                                                                                                                                                                                                                                         
################################################################################################################################################################################
# The function for Kmeans cluster id/centriod matching
                                              
cluster_label <- function(input_matrix, centers) 
    {
    
    # compute squared euclidean distance from each sample to each cluster center
    tmp <- sapply(seq_len(nrow(input_matrix)),
                function(i) apply(centers, 1,
                            function(v) sum((input_matrix[i, ]-v)^2)))
                                  
    # find index of min distance between sample and centroid
    max.col(-t(tmp)) 
    }
                                  
################################################################################################################################################################################
# The function for Kmeans cluster id/centriod matching
                                  
cluster_relabel <- function(centers1, centers2) 
    {
    # compute squared euclidean distance from each 2nd centroid to each 1st centroid
    tmp <- sapply(seq_len(nrow(centers2)),
                function(i) apply(centers1, 1,
                            function(v) sum((centers2[i, ]-v)^2)))
                                  
    # find index of min distance between 1st centroid and 2nd centroid
    cbind(1:dim(centers2)[1],max.col(-t(tmp)))  
    }                                  
 
