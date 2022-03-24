    kk_x_list<-2:12
    delta_list<-c( 0.01,seq(0,1,by=0.1),0.99)   
    km_name<-'genomic_vcf_'
    
# work flow   
if (opc=="g") {                                
    km_name<-paste(km_name,opc,sep="")  
    save.image(file=paste('germ_somatic_vcf_gene_',opc,'.RData' ,sep="") )        
    germ_somatic_vcf_reg<-distance_L1_GO_regulation_g(germ_somatic_vcf, marker_index, -6)    
    germ_somatic_vcf_reg_cluster<-supervised_clustering(germ_somatic_vcf_reg, marker_index, k_folds, marker_cutoff_metrics,km_name)      
    germ_somatic_vcf_reg_cluster_qc<-score_qc_all(outdir,km_name,marker_cutoff_metrics, kk_x_list, delta_list,k_folds)
    save.image(file=paste(outdir,'/',km_name,opc,'.RData' ,sep="") )        
} else if (opc=="v") {                                
    km_name<-paste(km_name,opc,sep="")  
    save.image(file=paste('germ_somatic_vcf_gene_',opc,'.RData' ,sep="") )        
    germ_somatic_vcf_reg<-distance_L1_GO_regulation_g(germ_somatic_vcf, marker_index, -6)    
    germ_somatic_vcf_reg_cluster<-supervised_clustering(germ_somatic_vcf_reg, marker_index, k_folds, marker_cutoff_metrics,km_name)      
    germ_somatic_vcf_reg_cluster_qc<-score_qc_all(outdir,km_name,marker_cutoff_metrics, kk_x_list, delta_list,k_folds)
    save.image(file=paste(outdir,'/',km_name,opc,'.RData' ,sep="") )         
} 
