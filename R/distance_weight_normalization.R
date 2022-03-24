# The function for calculate weight factors based on OR from statistical associations

squared_log_Odd_ratio<-function(input_matrix)
    {
    colnames(input_matrix)<-c('exp','res')
    
    # regression coefficients
    model1 <- glm(res ~ exp, data=input_matrix, family='gaussian')
    model1_tmp<-beta(model1)
    OR<-model1_tmp$coefficients[2,1]
    
    # log transform
    abs_OR<-abs(log(abs(OR),exp(1)))
    
    return(sqrt(abs_OR))  
    }

# The function for calculate weight factors based on HR from statistical associations

squared_log_Hazard_ratio<-function(input_matrix)
    {
    colnames(input_matrix)<-c('exp','intxsurv','dead')
    
    # regression coefficients
    model1 <- coxph(Surv(input_matrix$intxsurv, input_matrix$dead) ~ exp, data=input_matrix)
    abs_HR<-abs(summary(model1)$coef[1,1])
    
    return(sqrt(abs_HR))  
    }

################################################################################################################################################################################
# The function for calculate weighted matrix from statistical associations

distance_weight_normalization<- function(input_matrix,marker_index)
    {
    lr<-dim(input_matrix)[1]
    lc<-dim(input_matrix)[2]
    wc<-sum(input_matrix[,1:(marker_index-1)]==1)       
    feature_weight_idx<-c()
    dis_wt_matrix<-data.frame(matrix(nrow=1,ncol=dim(input_matrix)[2]))
    colnames(dis_wt_matrix)<-colnames(input_matrix)
    
    # calculate the feature weight and apply
    for (i in 1:(marker_index-1))
        {
        feature_weight_stat<-squared_log_Hazard_ratio(input_matrix[,c(i,marker_index,marker_index+1)])  
        if (is.na(feature_weight_stat))
            {feature_weight_stat<-0}            
        input_matrix[,i]<-sapply(input_matrix[,i],function(x) x*feature_weight_stat)
        dis_wt_matrix[,i]<-feature_weight_stat                         
        } 
                                 
    print("# calculate the sample stat norm and apply")    
    data_norm_stat<-sum(input_matrix[,1:(marker_index-1)])/wc #(lr*(lc-1))
    print(data_norm_stat)                                  
    input_matrix[,1:(marker_index-1)]<- apply(input_matrix[,1:(marker_index-1)],1:2, function(x) x/data_norm_stat)      
                                              
    return(list(input_matrix,dis_wt_matrix))
    }
        
