FROM asachet/shiny-base:R3.6.1

# MAINTAINER Tao Zhang “tzhang@nmdp.org"

#FROM rocker/r-base:3.6.3

RUN apt-get update && apt-get install -y openjdk-8-jdk libssl-dev libxml2-dev libz-dev libfreetype6-dev xorg libx11-dev libglu1-mesa-dev && apt-get autoremove
RUN R CMD javareconf
ENV PATH /usr/local/bin/R:$PATH 

RUN install2.r rgl data.table XML ggplot2 ggforce gridExtra BiocManager optparse reghelper glmnet pheatmap Hmisc survminer hash doParallel mclust fastmap htmltools randomForestSRC pacman RColorBrewer mlr
RUN R -e 'install.packages("XML")'
RUN R -e 'BiocManager::install("bios2mds")'
RUN R -e 'BiocManager::install(c("GOSemSim","org.Hs.eg.db"))'
RUN R -e 'BiocManager::install("biomaRt",ask = T,force=T)'
RUN R -e 'pacman::p_load_gh("IyarLin/survXgboost")'
RUN R -e 'pacman::p_load("survival")'
RUN R -e 'pacman::p_load("xgboost")'



COPY /R /Supervised-clustering-survival/R
COPY /Example /Supervised-clustering-survival/Example

EXPOSE 80
CMD ["/bin/sh"]
#CMD ["/Supervised-clustering-survival/Example/test.sh"]
