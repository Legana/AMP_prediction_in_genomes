
# Benchmarking antimicrobial peptide (AMP) machine learning models in a genome-scanning context

This repository looks at how best to predict AMPs in genomes and how to
assess prediction performance.

-   [Comparison of training and test data used by AMP
    predictors](01_data_comparison.md) : Rmd file
    [01_data_comparison.Rmd](01_data_comparison.Rmd)
-   [Performance metrics in an omics-scanning
    context](02_performance_metrics.md) : Rmd file
    [02_performance_metrics.Rmd](02_performance_metrics.Rmd)
-   [Benchmarking AMP predictors](03_benchmarking.md) : Rmd file
    [03_benchmarking.Rmd](03_benchmarking.Rmd)

The files required to run the code in these Rmd files can be obtained
[here](https://cloudstor.aarnet.edu.au/plus/s/Hd51gUnXdCq0nEg) or by
using the command:

``` bash
wget 'https://cloudstor.aarnet.edu.au/plus/s/Hd51gUnXdCq0nEg/download' -O data.tgz
tar -zxvf data.tgz 
```

### `sessionInfo()`

    R version 4.1.2 (2021-11-01)
    Platform: x86_64-apple-darwin17.0 (64-bit)
    Running under: macOS Monterey 12.3.1

    Matrix products: default
    LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

    locale:
    [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] caret_6.0-90       lattice_0.20-45    stringdist_0.9.8   ComplexUpset_1.3.3 patchwork_1.1.1   
     [6] forcats_0.5.1      stringr_1.4.0      dplyr_1.0.7        purrr_0.3.4        readr_2.1.1       
    [11] tidyr_1.1.4        tibble_3.1.6       ggplot2_3.3.5      tidyverse_1.3.1    ampir_1.1.0       

    loaded via a namespace (and not attached):
     [1] httr_1.4.2           jsonlite_1.7.2       splines_4.1.2        foreach_1.5.1        prodlim_2019.11.13  
     [6] modelr_0.1.8         assertthat_0.2.1     stats4_4.1.2         cellranger_1.1.0     yaml_2.2.1          
    [11] globals_0.14.0       ipred_0.9-12         pillar_1.6.4         backports_1.4.1      glue_1.6.0          
    [16] pROC_1.18.0          digest_0.6.29        rvest_1.0.2          colorspace_2.0-2     recipes_0.1.17      
    [21] htmltools_0.5.2      Matrix_1.3-4         plyr_1.8.6           timeDate_3043.102    pkgconfig_2.0.3     
    [26] broom_0.7.11         listenv_0.8.0        haven_2.4.3          scales_1.1.1         gower_0.2.2         
    [31] lava_1.6.10          tzdb_0.2.0           generics_0.1.1       ellipsis_0.3.2       withr_2.4.3         
    [36] nnet_7.3-16          cli_3.1.0            survival_3.2-13      magrittr_2.0.1       crayon_1.4.2        
    [41] readxl_1.3.1         evaluate_0.14        fs_1.5.2             future_1.23.0        fansi_1.0.2         
    [46] parallelly_1.30.0    nlme_3.1-153         MASS_7.3-54          xml2_1.3.3           class_7.3-19        
    [51] tools_4.1.2          data.table_1.14.2    hms_1.1.1            lifecycle_1.0.1      reprex_2.0.1        
    [56] munsell_0.5.0        compiler_4.1.2       rlang_0.4.12         grid_4.1.2           rstudioapi_0.13     
    [61] iterators_1.0.13     rmarkdown_2.13       gtable_0.3.0         ModelMetrics_1.2.2.2 codetools_0.2-18    
    [66] DBI_1.1.2            reshape2_1.4.4       R6_2.5.1             lubridate_1.8.0      knitr_1.37          
    [71] fastmap_1.1.0        future.apply_1.8.1   utf8_1.2.2           stringi_1.7.6        parallel_4.1.2      
    [76] Peptides_2.4.4       Rcpp_1.0.8           vctrs_0.3.8          rpart_4.1-15         dbplyr_2.1.1        
    [81] tidyselect_1.1.1     xfun_0.30     
