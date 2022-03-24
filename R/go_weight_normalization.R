###############################################################################################################################################################################
# 1.3 function for calculate weight factors from gene ontology itemized overpresentation
                                              
hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="MF", computeIC=FALSE)  
gene_go_fraction<- function(cohort_go)
        {
        genes <- cohort_go
        gene_no<-length(cohort_go)
    
        # calculate gene semantic similarity score    
        cohort_go_similariy<-mgeneSim(genes, semData=hsGO2, measure="Wang", combine="BMA", verbose=FALSE)   
        gene_semantic_similariy<-apply(cohort_go_similariy,2,function(x) (sum(x)-1)/(gene_no-1))
                                       
        return(as.data.frame(t(gene_semantic_similariy)))
        } 
                                            
################################################################################################################################################################################            
# The function for calculate weighted matrix from statistical associations of gene based matrix
                                              
go_weight_normalization_g<- function(input_matrix,marker_index)
    {
    lr<-dim(input_matrix)[1]
    lc<-dim(input_matrix)[2]
    wc<-sum(input_matrix[,1:(marker_index-1)]==1)
    feature_weight_idx<-c()
    go_wt_matrix<-data.frame(matrix(nrow=1,ncol=marker_index-1)) 
    colnames(go_wt_matrix)<-colnames(input_matrix)[1:(marker_index-1)]   
    
    # calculate the feature weight and apply
    all_wgs_gene<-colnames(input_matrix[,1:(marker_index-1)])
    all_wgs_gene_go<-all_wgs_gene[which(all_wgs_gene %in% unique(hsGO2@geneAnno$SYMBOL))]
    nk_gene_list<-all_wgs_gene[which(all_wgs_gene %nin% unique(hsGO2@geneAnno$SYMBOL))]   
    all_wgs_gene_go_new <- hash(as.character(nk_gene_list),rep(0,length(nk_gene_list)))
    for (gene in nk_gene_list)
        {
        if (grepl('-',gene))
            {
            gene_list<-unlist(str_split(gene,'-'))
            gene_list_all<-c()
            for (g in 1:length(gene_list))
                {
                gene_list_all<-c(paste(gene_list[(g-1):g],collapse='-',gene_list[g]))
                }
            if (gene_list_all %in% unique(hsGO2@geneAnno$SYMBOL))
            {
            flush.console()   
            #print(gene_list)                
              all_wgs_gene_go_new[[gene]] <-gene_list_all[which(gene_list_all %in% unique(hsGO2@geneAnno$SYMBOL))][1]
            }
        }
    }
    for (g in keys(all_wgs_gene_go_new))
        {
            if (all_wgs_gene_go_new[[g]]!=0)
                {
                all_wgs_gene_go<-c(all_wgs_gene_go,all_wgs_gene_go_new[[g]])
            } else {
                del(g, all_wgs_gene_go_new)
            }
        }
    if (length(all_wgs_gene_go)<2)
        {
        flush.console()
        print("No enough GO genes!!!")
        stop()
        }
    all_wgs_gene_go<-unique(all_wgs_gene_go)
    
    # calculate the feature weight and apply      
    feature_weight_go<-gene_go_fraction(all_wgs_gene_go) 
    
    # check gene GO score availability
    nk_gene_go_list<-values(all_wgs_gene_go_new)[values(all_wgs_gene_go_new) %in% colnames(feature_weight_go)]    
    for (i in 1:(marker_index-1))
        {gg<-colnames(input_matrix)[i]
        if ( gg %in% colnames(feature_weight_go))
            {
            input_matrix[,i]<-sapply(input_matrix[,i],function(x) x*feature_weight_go[,which(colnames(feature_weight_go)==gg)])
            go_wt_matrix[1,i]<-feature_weight_go[,which(colnames(feature_weight_go)==gg)]                                  
        } else if (gg %in% keys(all_wgs_gene_go_new)) {
            if (all_wgs_gene_go_new[[gg]] %in% nk_gene_go_list) 
                {
                nk_gene<-all_wgs_gene_go_new[[gg]]
                input_matrix[,i]<-sapply(input_matrix[,i],function(x) x*feature_weight_go[,which(colnames(feature_weight_go)==nk_gene)])
                go_wt_matrix[1,i]<-feature_weight_go[,which(colnames(feature_weight_go)==nk_gene)]  
                }                         
        } else {
            input_matrix[,i]<-0
            go_wt_matrix[1,i]<-0        
            }
        }
    go_wt_matrix[is.na(go_wt_matrix)]<-0  
                                         
    print("# calculate the sample norm and apply")    
    data_norm_go<-sum(input_matrix[,1:(marker_index-2)])/wc #lc*lr
    print(data_norm_go)    
    # check non-GO gene norm                               
    if(data_norm_go==0)
        {
        data_norm_go<-0.0000001
        }                                   
    input_matrix[,1:(marker_index-1)]<- apply(input_matrix[,1:(marker_index-1)],2, function(x) x/data_norm_go)                    
    return(list(input_matrix,go_wt_matrix))
    }

################################################################################################################################################################################    
# The function for calculate weighted matrix from statistical associations of variant based matrix
                                       
go_weight_normalization_v<- function(input_matrix,marker_index)
    {
   lr<-dim(input_matrix)[1]
    lc<-dim(input_matrix)[2]
    wc<-sum(input_matrix[,1:(marker_index-1)]==1)
    feature_weight_idx<-c()
    go_wt_matrix<-data.frame(matrix(nrow=1,ncol=marker_index-1)) 
    colnames(go_wt_matrix)<-colnames(input_matrix)[1:(marker_index-1)]  
    
    # translate variant id into gene id
    all_wgs_variant_list<-colnames(input_matrix[,1:(marker_index-1)])
    if (grepl('\\:', all_wgs_variant_list[1]))
        {ix<-'\\:' 
    } else if (grepl('\\.', all_wgs_variant_list[1]))   
        {ix<-'\\.'}
    all_wgs_variant_list_sim<-sapply(all_wgs_variant_list,function(x) paste(unlist(str_split(x,ix))[1],unlist(str_split(x,ix))[2],sep=':'))  
    all_wgs_gene<-unique(variant_gene_id_dict[variant_gene_id_dict$V1 %in% all_wgs_variant_list_sim,3])      
    all_wgs_gene_go<-all_wgs_gene[which(all_wgs_gene %in% unique(hsGO2@geneAnno$SYMBOL))]
                                     
    # retreive gene id from intergenic variants                                 
    nk_gene_list<-all_wgs_gene[which(all_wgs_gene %nin% unique(hsGO2@geneAnno$SYMBOL))]   
    all_wgs_gene_go_new <- hash(as.character(nk_gene_list),rep(0,length(nk_gene_list)))                               
    for (gene in nk_gene_list)
        {
        if (grepl('-',gene))
            {
            gene_list<-unlist(str_split(gene,'-'))
            if (gene_list %in% unique(hsGO2@geneAnno$SYMBOL))
                {               
              all_wgs_gene_go_new[[gene]] <-gene_list[which(gene_list %in% unique(hsGO2@geneAnno$SYMBOL))][1]
                }
            }
        }
    for (g in keys(all_wgs_gene_go_new))
        {
            if (all_wgs_gene_go_new[[g]]!=0)
                {
                all_wgs_gene_go<-c(all_wgs_gene_go,all_wgs_gene_go_new[[g]])
            } else {
                del(g, all_wgs_gene_go_new)
            }
        }
    if (length(all_wgs_gene_go)<2)
        {
        flush.console()
        print("No enough GO genes!!!")
        stop()
        }
    all_wgs_gene_go<-unique(c(all_wgs_gene_go,values(all_wgs_gene_go_new)))
                                     
    # calculate the feature weight and apply                                     
    feature_weight_go<-gene_go_fraction(all_wgs_gene_go) 
                                     
    # check gene GO score availability
    nk_gene_go_list<-values(all_wgs_gene_go_new)[values(all_wgs_gene_go_new) %in% colnames(feature_weight_go)]                                     
    for (i in 1:(marker_index-1))
     {gg<-unique(variant_gene_id_dict[variant_gene_id_dict$V1==paste(unlist(str_split(colnames(input_matrix)[i],ix))[1],unlist(str_split(colnames(input_matrix)[i],ix))[2],sep=':'),3])
         if ( gg %in% colnames(feature_weight_go))
            {
            input_matrix[,i]<-sapply(input_matrix[,i],function(x) x*feature_weight_go[,which(colnames(feature_weight_go)==gg)])
            go_wt_matrix[1,i]<-feature_weight_go[1,which(colnames(feature_weight_go)==gg)]                                  
        } else if (gg %in% keys(all_wgs_gene_go_new)) {
            if (all_wgs_gene_go_new[[gg]] %in% nk_gene_go_list) {
                nk_gene<-all_wgs_gene_go_new[[gg]]
                input_matrix[,i]<-sapply(input_matrix[,i],function(x) x*feature_weight_go[,which(colnames(feature_weight_go)==nk_gene)])
                go_wt_matrix[1,i]<-feature_weight_go[1,which(colnames(feature_weight_go)==nk_gene)]  }                         
        } else {
            input_matrix[,i]<-0
            go_wt_matrix[1,i]<-0        
            }
        }
    go_wt_matrix[is.na(go_wt_matrix)]<-0  
                                         
    print("# calculate the sample go norm and apply")    
    data_norm_go<-sum(input_matrix[,1:(marker_index-2)])/wc #lc*lr
    print(data_norm_go)    
                                         
    # check and apply non-GO gene norm                               
    if(data_norm_go==0)
        {
        data_norm_go<-0.0000001
        }                                   
    input_matrix[,1:(marker_index-1)]<- apply(input_matrix[,1:(marker_index-1)],2, function(x) x/data_norm_go)        
                                              
    return(list(input_matrix,go_wt_matrix))
    }
