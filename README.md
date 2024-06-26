Supervised clustering survival
=======================

This is supervised clustering with integrated weighting of both GO semantic similarity and statistical association effect, and mine the mutational clusters with strong survival outcome stratifications. 

## Workflow:
![WGS_WORKFLOW](https://github.com/tzhang-nmdp/Supervised-clustering-survival/blob/main/Example/workflow.png)


## Algorithm:
![PSEUDOCODE](https://github.com/tzhang-nmdp/Supervised-clustering-survival/blob/main/Example/Pseudocode.PNG)

## Installation:

``` r
# Docker installation
docker build -t supervised-clustering-survival:main .

# R installation
library(devtools)
install_github("tzhang-nmdp/Supervised-clustering-survival@main")
```

## Usage:

``` r
Rscript Supervised-clustering-survival/R/SCCW_supervised_clustering.R \
-i  ${genomic_data}.RData \ # input matrix for genomic data
-o ***variant/gene \ # running model option ( 'variant' for variant level of common variant analysis, 'gene' for gene level of rare variant analysis)
-d ${outdir} \  # output directory
-k 5 # k_fold setting for cross-validation
```

## Example

``` r
Rscript Supervised-clustering-survival/R/SCCW_supervised_clustering.R \
-i Supervised-clustering-survival/Example/genomic_data_vcf.RData \ # input matrix for genomic data (simulated example for demo only)
-o test_gene \ # running model option ( 'variant' for variant level of common variant analysis, 'gene' for gene level of rare variant analysis)
-d Supervised-clustering-survival/Example \  # output directory
-k 5 # k_fold setting for cross-validation
```

### input data  (to be created from VCF and survival outcome information)
![INPUT](https://github.com/tzhang-nmdp/Supervised-clustering-survival/blob/main/Example/Input.PNG)
### Variant-gene-dictionary (to be created from gene variant annotation)
![Variant-gene-dictionary](https://github.com/tzhang-nmdp/Supervised-clustering-survival/blob/main/Example/variant_gene_dict.png)

### output data
1. Clustering model file: ***.model.RData
2. Model metrics file: ***.csv
3. Survival plot file: ***.tiff

### Citation
Please cite our paper in *Journal of Hematology & Oncology*:
[Whole-genome sequencing identifies novel predictors for hematopoietic cell transplant outcomes for patients with myelodysplastic syndrome: a CIBMTR study](https://jhoonline.biomedcentral.com/articles/10.1186/s13045-023-01431-7).
