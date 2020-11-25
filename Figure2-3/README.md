# Creating Figure 2 and 3

Download `20201020_All_motif_activations_HDF5_combined_matfet_pool120_5u5d_600mot.h5` from 
[this page](https://scover-figures-data.cog.sanger.ac.uk/browser.html) and place in `scover_output` folder.

Make sure the right R package versions are installed. See list of packages that I am using at the end of this page. 

Then, run `Rscript Figure2Figure3.R` to create the plots in the `output` directory. 


## sessionInfo()

I am using these R packages:

```
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 18.04.3 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas/liblapack.so.3
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] patchwork_1.0.0.9000        pheatmap_1.0.12            
##  [3] reshape2_1.4.4              tidytext_0.2.6             
##  [5] gghalves_0.1.0              ggthemes_4.2.0             
##  [7] ggrepel_0.8.2               stringr_1.4.0              
##  [9] rhdf5_2.30.1                igraph_1.2.5               
## [11] ggplot2_3.3.1               SingleCellExperiment_1.8.0 
## [13] SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
## [15] BiocParallel_1.20.1         matrixStats_0.56.0         
## [17] Biobase_2.46.0              GenomicRanges_1.38.0       
## [19] GenomeInfoDb_1.22.1         IRanges_2.20.2             
## [21] S4Vectors_0.24.4            BiocGenerics_0.32.0        
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.1.0       xfun_0.14              purrr_0.3.4           
##  [4] lattice_0.20-40        colorspace_1.4-1       vctrs_0.3.1           
##  [7] generics_0.0.2         SnowballC_0.7.0        htmltools_0.4.0       
## [10] yaml_2.2.1             rlang_0.4.6            pillar_1.4.4          
## [13] glue_1.4.1             withr_2.2.0            RColorBrewer_1.1-2    
## [16] plyr_1.8.6             GenomeInfoDbData_1.2.2 lifecycle_0.2.0       
## [19] zlibbioc_1.32.0        munsell_0.5.0          gtable_0.3.0          
## [22] evaluate_0.14          knitr_1.28             tokenizers_0.2.1      
## [25] Rcpp_1.0.4.6           scales_1.1.1           XVector_0.26.0        
## [28] digest_0.6.25          stringi_1.4.6          dplyr_1.0.0           
## [31] grid_3.6.3             tools_3.6.3            bitops_1.0-6          
## [34] magrittr_1.5           RCurl_1.98-1.2         tibble_3.0.1          
## [37] janeaustenr_0.1.5      crayon_1.3.4           pkgconfig_2.0.3       
## [40] ellipsis_0.3.1         Matrix_1.2-18          rmarkdown_2.2         
## [43] Rhdf5lib_1.8.0         R6_2.4.1               compiler_3.6.3
```
