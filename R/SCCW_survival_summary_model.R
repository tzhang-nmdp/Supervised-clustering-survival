# survival summary model with permutations
survival_summary_model<-function(input_matrix,kk_x,delta, km_name)
{
    cindex_matrix<-data.frame(matrix(ncol = 3, nrow =3))
    colnames(cindex_matrix)<-c('RSF','coxph','xgboost')
    rownames(cindex_matrix)<-c('cindex_median','cindex_CI_lb','cindex_CI_ub')
    # run Random forest survival model with 100 permutation random status  
    print("# run the importance and prediction matrix with permutations and cross-validations")

    imp_dict_494_all<-data.frame(matrix(ncol = 1, nrow = dim(input_matrix)[2]-2))
    colnames_tmp<-colnames(input_matrix)[order(colnames(input_matrix))]
    imp_dict_494_all[,1]<-colnames_tmp[which(colnames_tmp %nin% c("intxsurv", "dead"))]
    cidx_494_rsf_all<-c()
    tempSubset_fs_tmp<-input_matrix
    for (i in 1:100)
        {
        cat(i)
        flush.console()
        set.seed(i)
        index <- sample(1:dim(tempSubset_fs_tmp)[1], replace=TRUE, size=dim(tempSubset_fs_tmp)[1])
        new.data <- tempSubset_fs_tmp[index,]   
        modRFSRC_494_all <- rfsrc(Surv(intxsurv, dead) ~ ., data=new.data,block.size=1,err.block=1,  na.action = "na.impute",importance = TRUE)
        imp_494_all<-as.data.frame(modRFSRC_494_all$importance)
        imp_494_all[,2]<-rownames(imp_494_all)
        imp_494_all<-imp_494_all[order(rownames(imp_494_all)),]
        imp_dict_494_all<-cbind(imp_dict_494_all,imp_494_all[,1]) 
        c_index<-c(1-median(modRFSRC_494_all$err.rate))
        cidx_494_rsf_all<-rbind(cidx_494_rsf_all,c_index) 
        }
    print(modRFSRC_494_all)

    # Random forest survival model performance check 
    cindex_matrix[,1]<-c(median(cidx_494_rsf_all),quantile(cidx_494_rsf_all, probs=c(0.05, 0.95)))

    # plot the feature importance
    imp_dict_494_all[,2]<-apply(imp_dict_494_all[,2:dim(imp_dict_494_all)[2]],1,mean)
    imp_dict_494_all<-imp_dict_494_all[,1:2]
    imp_dict_494_all$V3<-'mutational'
    colnames(imp_dict_494_all)<-c('feature_name','feature_importance','feature_category')
    imp_dict_494_all$feature_name<-as.character(imp_dict_494_all$feature_name)
    for (i in 1:dim(imp_dict_494_all)[1])
        {
    if  (grepl('_id',imp_dict_494_all[i,2]))
        {
        imp_dict_494_all[i,3]<-'clustering'
    }
        else if  ( imp_dict_494_all[i,2] %in% c('ipssr','HMA', 'CHEMO','mdstype'))
        {
        imp_dict_494_all[i,3]<-'clinical'
    }    
    }
    imp_dict_494_all$feature_importance_rev<-sapply(imp_dict_494_all$feature_importance, function(x) { x*(-1)})
    #imp_dict_494_all[imp_dict_494_all$feature_name=='km_id_com', 'feature_name']<-'cluster_common_variant'
    #imp_dict_494_all[imp_dict_494_all$feature_name=='km_id', 'feature_name']<-'cluster_rare_variant'
    #imp_dict_494_all[imp_dict_494_all$feature_name=='rc_id', 'feature_name']<-'cluster_recurrent_variant'
    options(repr.plot.width=20, repr.plot.height=10)
    p<-ggplot(imp_dict_494_all,aes(x=reorder(feature_name,feature_importance_rev),y=feature_importance,fill=feature_category)) + geom_bar(colour="black",stat="identity", width=.5, position = "stack",alpha=1)+theme(axis.text.x=element_text(angle=90,size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),legend.text=element_text(size=15))
    tiff(paste(outdir,"/",km_name,"_kk_x_",kk_x,"_delta_",delta,"_RSF_importance.tiff",sep=""),width = 2000, height =1500,res=300) 
    p+ylab('Feature importance')+ scale_fill_manual(values=c('indianred2','green3','skyblue'))
    dev.off()


    
    tempSubset_fs_tmp <- input_matrix[complete.cases(input_matrix), ] # doesn't handle missing values for other models

    # run coxph survival model with 100 permutation random status  
    print("# run coxph model survival with permutations")
    cidx_494_cox_all<-c()
    for (i in 1:100)
        {
        cat(i)
        flush.console()
        set.seed(i)
        index <- sample(1:dim(tempSubset_fs_tmp)[1], replace=T, size=dim(tempSubset_fs_tmp)[1])
        new.data <- tempSubset_fs_tmp[index,]   
        modRFSRC_494_all <- coxph(Surv(new.data$intxsurv, new.data$dead) ~ ., data=new.data)
        c_index<-modRFSRC_494_all$concordance[6]
        cidx_494_cox_all<-rbind(cidx_494_cox_all,c_index) 
        }
    cindex_matrix[,2]<-c(median(cidx_494_cox_all),quantile(cidx_494_cox_all, probs=c(0.05, 0.95)))

    # run xgboost survival model with 100 permutation random status  
    print("# run xgboost model survival with permutations")    
    tempSubset_fs_tmp[,colnames(tempSubset_fs_tmp)!='dead']<-apply(tempSubset_fs_tmp[,colnames(tempSubset_fs_tmp)!='dead'],2,as.numeric)    
    label <- ifelse(tempSubset_fs_tmp$dead == 1, tempSubset_fs_tmp$intxsurv, -tempSubset_fs_tmp$intxsurv)
    cidx_494_xgb_all<-c()
    for (i in 1:100)
        {    cat(i)
        flush.console()
    set.seed(i)
    val_ind <- sample.int(nrow(tempSubset_fs_tmp), 0.1 * nrow(tempSubset_fs_tmp))
    x_train <- as.matrix(tempSubset_fs_tmp[-val_ind, !names(tempSubset_fs_tmp) %in% c("intxsurv", "dead")])
    x_label <- label[-val_ind]
    x_val <- xgb.DMatrix(as.matrix(tempSubset_fs_tmp[val_ind, !names(tempSubset_fs_tmp) %in% c("intxsurv", "dead")]),label = label[val_ind])

    surv_xgboost_model <- xgb.train.surv(
    params = list(
        objective = "survival:cox",
        eval_metric = "cox-nloglik",
        eta = 0.05 # larger eta leads to algorithm not converging, resulting in NaN predictions
    ), data = x_train, label = x_label,
    watchlist = list(val2 = x_val),
    nrounds = 1000, early_stopping_rounds = 30
    )
    risk_scores <- predict(object = surv_xgboost_model, newdata = as.matrix(tempSubset_fs_tmp[, !names(tempSubset_fs_tmp) %in% c("intxsurv", "dead")]), type = "risk")
    cidx_494_xgb_all<-c(cidx_494_xgb_all,concordance_index(tempSubset_fs_tmp$intxsurv,1-risk_scores))
        }

    cindex_matrix[,3]<-c(median(cidx_494_xgb_all),quantile(cidx_494_xgb_all, probs=c(0.05, 0.95)))

    return(cindex_matrix)
        }


