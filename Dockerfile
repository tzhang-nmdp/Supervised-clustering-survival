FROM dockerhub.nmdp.org:8443/nmdp-shiny-base:3.6.1

MAINTAINER Tao Zhang “tzhang@nmdp.org"

RUN apt-get update && apt-get install -y openjdk-8-jdk libssl-dev && apt-get autoremove

RUN R CMD javareconf

RUN install2.r ggplot2 ggforce gridExtra BiocManager optparse reghelper glmnet pheatmap Hmisc survminer hash
RUN R -e 'BiocManager::install(c("bios2mds","GOSemSim","org.Hs.eg.db","biomaRt"))'

COPY /R /Supervised-clustering-survival/R/

EXPOSE 80

CMD ["/usr/local/bin/Rscript /Supervised-clustering-survival/R/SCCW_supervised_clustering.R"]
