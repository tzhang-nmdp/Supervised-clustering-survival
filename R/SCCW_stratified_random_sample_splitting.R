################################################################################################################################################################################ 
# The function for randomally split 5 folds matrix (asign one as out-of-bag validation)
                                                                                      
random_splits <- function (input_matrix,marker_index, k_folds)
    {
    ln0<-dim(input_matrix[input_matrix[,marker_index]==0, ])[1]
    ln1<-dim(input_matrix[input_matrix[,marker_index]==1, ])[1]    
    sub0<-round(ln0/(k_folds),0)
    sub1<-round(ln1/(k_folds),0)    
    input_matrix$io<-k_folds
    sample_list<-data.frame('index'=1:dim(input_matrix)[1],clin_var= input_matrix[,marker_index])
    
    # sampling for k folds without replacement
    for (i in 1:(k_folds-1))
        {
        set.seed(1)   
        #oob<-sample(sample_list, replace=F, size=sub)
        #oob_tmp<-sample_list %>% group_by(clin_var) %>% sample_n(sub)
        oob_tmp0<-sample(sample_list[sample_list$clin_var==0,'index'],sub0)
        oob_tmp1<-sample(sample_list[sample_list$clin_var==1,'index'],sub1)     
        # track the k folds
        input_matrix[c(oob_tmp0,oob_tmp1),'io']<-i
        # subset samples
        ln_tmp<-sample_list[sample_list$index %nin% c(oob_tmp0,oob_tmp1),]
        sample_list<-ln_tmp
        }
    
    return(input_matrix)
    }     
         
################################################################################################################################################################################ 
# The function for randomally split 5 folds matrix (asign one as out-of-bag validation)
                                              
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
