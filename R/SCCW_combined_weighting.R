################################################################################################################################################################################
# The function for calculate weighted matrix by weight factor delta 
                                              
matrix_weight_finals <- function (matrix_weight_stats, matrix_weight_go, delta, marker_index)
    {
    matrix_weight_final<-matrix_weight_stats
    # tranformed matrix by weight factor delta 
    matrix_weight_final[,1:(marker_index-1)]<-matrix_weight_stats[,1:(marker_index-1)]*(1-delta)+matrix_weight_go[,1:(marker_index-1)]*delta 
    return (matrix_weight_final)
    }      
          
################################################################################################################################################################################
# The function for calculate weighted matrix by weight factor delta 
                                              
matrix_weight_final <- function (matrix_weight_stats, matrix_weight_go, delta)
    {
    # tranformed matrix by weight factor delta 
    matrix_weight_final<-matrix_weight_stats*(1-delta)+matrix_weight_go*delta
    return (matrix_weight_final)
    }      
                       
