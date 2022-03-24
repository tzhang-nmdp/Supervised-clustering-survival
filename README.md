Supervised clustering survival
=======================

This is supervised clustering with integrated weighting of both GO semantic similarity and statistical association effect, and mine the mutational clusters with strong survival outcome stratifications. 

## Workflow:
![WGS_WORKFLOW](https://github.com/tzhang-nmdp/Supervised-clustering-survival/blob/main/Example/workflow.PNG)


## Algorithm:
![PSEUDOCODE](https://github.com/tzhang-nmdp/Supervised-clustering-survival/blob/main/Example/Pseudocode.PNG)

## Installation:

``` r
library(devtools)
install_github("tzhang-nmdp/Supervised-clustering-survival")
```

## Usage:

``` r
Rscript supervised_clustering_gene_variant-docu-L1-inb-Copy2.R \
-i  ${genomic_data}.RData \ # input matrix for genomic data
-o v/g \ # running model option ( 'v' for variant level of common variant analysis, 'g' for gene level of rare variant analysis)
-d ${outdir} \  # output directory
-k 5 # k_fold setting for cross-validation
```

## Example

### input data
![INPUT](https://github.com/tzhang-nmdp/Supervised-clustering-survival/blob/main/Example/Input.PNG)

### output data
1. Clustering model file: ***.model.RData
2. Model metrics file: ***.csv
3. Survival plot file: ***.tiff


