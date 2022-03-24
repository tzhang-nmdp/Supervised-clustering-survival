# The function for matrix regulated by alias/multicollinearity  

distance_regulation<- function(input_matrix,marker_index,n_lambda)
    {
    tmp <- cor(input_matrix)
    tmp[upper.tri(tmp)] <- 0
    diag(tmp) <- 0
    lc<-dim(input_matrix)[2]
    input_matrix <- input_matrix[, !apply(tmp, 2, function(x) any(abs(x) > 0.9, na.rm = TRUE))]
    marker_index<-dim(input_matrix)[2]-(lc-marker_index) 
    print(dim(input_matrix))
    # indentify alias/multicollinearity varaibles        
    ln_10k0<-ceiling((marker_index-1)/10000)  
    print(ln_10k0)
    if (ln_10k0>1)
        {
        rm<-c()
        for (i in 1:ln_10k0)
            {
            flush.console() 
            print(Sys.time())
            if (i!=ln_10k0)
                {
                x1<-input_matrix[,c((10000*(i-1)+1):(10000*i),marker_index)]
                } else {
                x1<-input_matrix[,c((10000*(i-1)+1):marker_index)]        
                }   
            cor_var<-alias(lm(intxsurv ~ ., data=x1))
            ld.vars <- attributes(cor_var$Complete)$dimnames[[1]]
            rm<-c(rm,ld.vars)
            }    
        
        #remove alias/multicollinearity varaibles         
        input_matrix_tmp<-input_matrix[,colnames(input_matrix) %nin% rm]
        set.seed(1)
        input_matrix_tmp<-input_matrix_tmp[,c(sample(1:(dim(input_matrix_tmp)[2]-2),dim(input_matrix_tmp)[2]-2),dim(input_matrix_tmp)[2]-1,dim(input_matrix_tmp)[2])]
        ln_10k<-ceiling((dim(input_matrix_tmp)[2]-2)/10000)    
        print(ln_10k) 
        
        # indentify and remove alias/multicollinearity varaibles with regrouping         
        if (ln_10k<ln_10k0) 
            {
            j=1       
            while (ln_10k>1 & j<=5)
                {
                rm<-c()      
                ln_10k<-ceiling((dim(input_matrix_tmp)[2]-2)/10000)  
                print(Sys.time())      
                print(ln_10k)
                for (i in 1:ln_10k)
                    {
                    if (i!=ln_10k)
                        {
                        x1<-input_matrix_tmp[,c((10000*(i-1)+1):(10000*i),dim(input_matrix_tmp)[2]-1)]
                        } else {
                        x1<-input_matrix_tmp[,c((10000*(i-1)+1):(dim(input_matrix_tmp)[2]-2),dim(input_matrix_tmp)[2]-1)]        
                        }   
                    cor_var<-alias(lm(intxsurv ~ ., data=x1))
                    ld.vars <- attributes(cor_var$Complete)$dimnames[[1]]
                    rm<-c(rm,ld.vars)
                    } 
                j=j+1
                input_matrix_tmp<-input_matrix_tmp[,colnames(input_matrix_tmp) %nin% rm]
                }
            }
        } else {
            x1<-input_matrix[,1:marker_index]
            cor_var<-alias(lm(intxsurv ~ ., data=x1))
            ld.vars <- attributes(cor_var$Complete)$dimnames[[1]]
            rm<-c(rm,ld.vars)  
            input_matrix_tmp<-input_matrix[,colnames(input_matrix) %nin% rm]
        }  
    return(input_matrix_tmp)
    }
