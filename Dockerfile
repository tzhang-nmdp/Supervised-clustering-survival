FROM dockerhub.nmdp.org:8443/nmdp-shiny-base:3.6.1
#FROM r-base:3.6.1 

MAINTAINER Tao Zhang “tzhang@nmdp.org"

RUN apt-get update && apt-get install -y openjdk-8-jdk libssl-dev libxml2-dev && apt-get autoremove
RUN R CMD javareconf

RUN install2.r ggplot2 ggforce gridExtra BiocManager optparse reghelper glmnet pheatmap Hmisc survminer hash XML doParallel mclust fastmap htmltools randomForestSRC pacman RColorBrewer mlr
RUN R -e 'BiocManager::install(c("bios2mds","GOSemSim","org.Hs.eg.db"))'
RUN R -e 'BiocManager::install("biomaRt",ask = T,force=T)'
RUN R -e 'pacman::p_load_gh("IyarLin/survXgboost")'
RUN R -e 'pacman::p_load("survival")'
RUN R -e 'pacman::p_load("xgboost")'

COPY /R /Supervised-clustering-survival/R
COPY /Example /Supervised-clustering-survival/Example

EXPOSE 80
CMD ["/bin/sh"]
#CMD ["/Supervised-clustering-survival/Example/test.sh"]
