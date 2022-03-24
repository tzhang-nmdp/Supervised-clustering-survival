# The function for calculate weighted matrix by weight factor delta 
                                              
matrix_weight_final <- function (matrix_weight_stats, matrix_weight_go, delta)
    {
    # tranformed matrix by weight factor delta 
    matrix_weight_final<-matrix_weight_stats*(1-delta)+matrix_weight_go*delta
    return (matrix_weight_final)
    }      
      
