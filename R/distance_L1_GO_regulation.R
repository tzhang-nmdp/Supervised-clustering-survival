# The function for matrix regulated by Lasso and GO of gene based matrix

distance_L1_GO_regulation_g<- function(input_matrix,marker_index,n_lambda)
    {

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

    print(paste("# lasso regression of data matrix: ", Sys.time(),sep="")) 
    feature_glm_lasso_reg<-c()
    ln_10k<-ceiling((dim(input_matrix_tmp)[2]-2)/10000) 
    for (i in 1:ln_10k)    
        {
        if (i!=ln_10k)
            {
            x1<-input_matrix_tmp[,c((10000*(i-1)+1):(10000*i),dim(input_matrix_tmp)[2]-1,dim(input_matrix_tmp)[2])]
            } else {
            x1<-input_matrix_tmp[,c((10000*(i-1)+1):(dim(input_matrix_tmp)[2]-2),dim(input_matrix_tmp)[2]-1,dim(input_matrix_tmp)[2])]        
            }          
        x <- model.matrix( ~.,x1[,1:(dim(x1)[2]-2)])
        y <- Surv(x1$intxsurv, x1$dead)
        fit <- glmnet(x, y, family="cox", alpha=0.5)
        Coefficients <- coef(fit, s = exp(n_lambda))
        lasso_coeff<-names(which(Coefficients[,1]!=0))
        lasso_coeff<-sub("`", "",lasso_coeff, fixed = TRUE)
        lasso_coeff<-sub("`", "",lasso_coeff, fixed = TRUE)
        feature_glm_lasso_reg<-c(feature_glm_lasso_reg,lasso_coeff)        
        feature_glm_lasso_reg<-c(feature_glm_lasso_reg,names(which(Coefficients[,1]!=0)))
        }
                                      
    print(paste("# apply lasso and go regulation: ", Sys.time(),sep=""))    
    input_matrix_L1_reg<-input_matrix_tmp[,c(which(colnames(input_matrix_tmp) %in% feature_glm_lasso_reg),dim(input_matrix_tmp)[2]-1, dim(input_matrix_tmp)[2])]       
    print(dim(input_matrix_L1_reg))
                                          
    print(paste("# exclude high correlated features: ", Sys.time(),sep="")) 
    tmp <- cor(input_matrix_L1_reg)
    tmp[upper.tri(tmp)] <- 0
    diag(tmp) <- 0
    lc<-dim(input_matrix_L1_reg)[2]
    input_matrix_L1_reg <- input_matrix_L1_reg[, !apply(tmp, 2, function(x) any(abs(x) > 0.9, na.rm = TRUE))]
    marker_index<-dim(input_matrix_L1_reg)[2]-(lc-marker_index) 
    print(dim(input_matrix_L1_reg))
                                                        
    return(input_matrix_L1_reg)
    }
                                                        
################################################################################################################################################################################
# The function for matrix regulated by Lasso and GO of variant based matrix

colrs=brewer.pal(7,"Set1")
distance_L1_GO_regulation_v<- function(input_matrix,marker_index,n_lambda)
    {

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
   
    print(paste("# elastic net regression of data matrix: ", Sys.time(),sep="")) 
    feature_glm_lasso_reg<-c()
    ln_10k<-ceiling((dim(input_matrix_tmp)[2]-2)/10000) 
    for (i in 1:ln_10k)    
        {
        if (i!=ln_10k)
            {
            x1<-input_matrix_tmp[,c((10000*(i-1)+1):(10000*i),dim(input_matrix_tmp)[2]-1,dim(input_matrix_tmp)[2])]
            } else {
            x1<-input_matrix_tmp[,c((10000*(i-1)+1):(dim(input_matrix_tmp)[2]-2),dim(input_matrix_tmp)[2]-1,dim(input_matrix_tmp)[2])]        
            }          
        x <- model.matrix( ~.,x1[,1:(dim(x1)[2]-2)])
        y <- Surv(x1$intxsurv, x1$dead)
        fit <- glmnet(x, y, family="cox", alpha=0.5)
        Coefficients <- coef(fit, s = exp(n_lambda))
        lasso_coeff<-names(which(Coefficients[,1]!=0))
        lasso_coeff<-sub("`", "",lasso_coeff, fixed = TRUE)
        lasso_coeff<-sub("`", "",lasso_coeff, fixed = TRUE)
        feature_glm_lasso_reg<-c(feature_glm_lasso_reg,lasso_coeff)        
        feature_glm_lasso_reg<-c(feature_glm_lasso_reg,names(which(Coefficients[,1]!=0)))
        }
                                      
    print(paste("# apply lasso and go regulation: ", Sys.time(),sep=""))    
    input_matrix_L1_reg<-input_matrix_tmp[,c(which(colnames(input_matrix_tmp) %in% feature_glm_lasso_reg),dim(input_matrix_tmp)[2]-1, dim(input_matrix_tmp)[2])]       
    print(dim(input_matrix_L1_reg))
    
    print(paste("# exclude high correlated features: ", Sys.time(),sep="")) 
    tmp <- cor(input_matrix_L1_reg)
    tmp[upper.tri(tmp)] <- 0
    diag(tmp) <- 0
    lc<-dim(input_matrix_L1_reg)[2]
    input_matrix_L1_reg <- input_matrix_L1_reg[, !apply(tmp, 2, function(x) any(abs(x) > 0.9, na.rm = TRUE))]
    marker_index<-dim(input_matrix_L1_reg)[2]-(lc-marker_index) 
    print(dim(input_matrix_L1_reg))
                                                        
    return(input_matrix_L1_reg)
    }  
    
