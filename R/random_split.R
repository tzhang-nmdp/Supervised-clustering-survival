# 1.6 function for randomally split 5 folds matrix (asign one as out-of-bag validation)
                                              
random_split <- function (input_matrix,k_folds)
    {
    ln<-dim(input_matrix)[1]
    sub<-round(ln/k_folds,0)
    input_matrix$io<-k_folds
    sample_list<-1:ln
    
    # sampling for k folds without replacement
    for (i in 1:(k_folds-1))
        {
        set.seed(1)        
        oob<-sample(sample_list, replace=F, size=sub)
        # track the k folds
        input_matrix[oob,'io']<-i
        
        # subset samples
        ln_tmp<-sample_list[sample_list %nin% oob]
        sample_list<-ln_tmp
        }
    
    return(input_matrix)
    }    
