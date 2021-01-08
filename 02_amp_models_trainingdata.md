Training datasets in other AMP predictors
================

``` r
library(ampir)
library(tidyverse)
```

**AmPEP Training Data**

The AmPEP 2018 AMP predictor provides its training data available
directly for download from
<https://cbbio.cis.um.edu.mo/software/AmPEP/>. The length distribution
of sequences in this database is interesting. It shows that sequences
classified as AMPs form a clear peak corresponding to mature peptides
whereas non-AMP (background) sequences are clearly larger and more
likely to represent full length proteins.

``` r
ampep_data <- read_faa("data/amp_predictors/amPEP/M_model_train_nonAMP_sequence.fasta") %>% add_column(class="non-AMP") %>% 
  rbind(read_faa("data/amp_predictors/amPEP/M_model_train_AMP_sequence.fasta") %>% add_column(class="AMP")) %>% 
  mutate(length = nchar(seq_aa)) %>% add_column(predictor="AmPEP") %>%
  mutate(seq_name = paste0("amPEP_trainset_neg", 1:n()))
```

AmPEP was redesigned in 2020 as Deep-AmPEP30 to focus on short AMPs ( \<
30 amino acids) [Yan et al](https://doi.org/10.1016/j.omtn.2020.05.006)
and its training and test data is available
[here](https://cbbio.online/AxPEP/?action=dataset).

``` r
deep_ampep_data <- read_faa("data/amp_predictors/deepamPEP30/train_ne.fasta") %>%
   rbind(read_faa("data/amp_predictors/deepamPEP30/test_ne.fasta")) %>%
  add_column(class = "non-AMP") %>%
  rbind(read_faa("data/amp_predictors/deepamPEP30/train_po.fasta") %>%
  rbind(read_faa("data/amp_predictors/deepamPEP30/test_po.fasta")) %>%
          add_column(class = "AMP")) %>%
  mutate(dataset = case_when(
    str_detect(seq_name, "^test") ~ "Test",
    str_detect(seq_name, "^uni") ~ "Test",
                             TRUE ~ "Train")) %>%
  mutate(length = nchar(seq_aa)) %>% 
  add_column(predictor="deep_AmPEP")
```

AmPEP was additionally created as a python application, amPEPpy, by
[Lawrence et al. 2020](https://doi.org/10.1093/bioinformatics/btaa917).
amPEPpy’s training data originated from amPEP and were obtained via
amPEPpy’s [GitHub page](https://github.com/tlawrence3/amPEPpy).

``` r
ampeppy_data <- read_faa("data/amp_predictors/amPEPpy/M_model_train_nonAMP_sequence.numbered.proplen.subsample.fasta") %>% add_column(class="non-AMP") %>% 
  rbind(read_faa("data/amp_predictors/amPEPpy/M_model_train_AMP_sequence.numbered.fasta") %>% add_column(class="AMP")) %>% 
  mutate(length = nchar(seq_aa)) %>% add_column(predictor="AmPEPpy")
```

**AMP Scanner v2 Data**

AMP Scanner ([Veltri et
al. 2018](https://doi.org/10.1093/bioinformatics/bty179%5D)) data used
for training, testing and evaluation are available directly for download
from <https://www.dveltri.com/ascan/v2/about.html>. AMP Scanner’s
training data consisted of 1,066 AMP and non-AMP sequences. Their
testing data consisted of 712 AMP and non-AMP sequences.

``` r
ampscan_train_data <- read_faa("data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/AMP.tr.fa") %>%
   rbind(read_faa("data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/AMP.eval.fa")) %>%
  rbind(read_faa("data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/DECOY.tr.fa")) %>%
  rbind(read_faa("data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/DECOY.eval.fa")) %>%
             add_column(dataset="Train")

ampscan_test_data <- rbind(read_faa("data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/AMP.te.fa")) %>%
  rbind(read_faa("data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/DECOY.te.fa")) %>%
             add_column(dataset="Test")

ampscan_data <- rbind(ampscan_train_data, ampscan_test_data) %>%
  mutate(class = case_when(str_detect(seq_name, "^Uni") ~ "non-AMP", TRUE ~ "AMP")) %>%
    mutate(length = nchar(seq_aa)) %>% 
    add_column(predictor="AMP Scanner v2")
```

**iAMP-2L Data**

The benchmark data provided by [Xiao et
al. 2013](https://doi.org/10.1016/j.ab.2013.01.019) used for
[iAMP-2L](http://www.jci-bioinfo.cn/iAMP/data.html) has been used in
several studies to provide a somewhat independent estimate of prediction
accuracy. Their training data, or benchmark dataset as they termed it,
comprises of 897 AMPs and 2,405 non-AMPs. Their test or independent
dataset comprises of 920 AMPs and 920 non-AMPs.

``` r
iAMP2L_train <- read_faa("data/amp_predictors/iAMP-2L/xiao_benchmark.fasta") %>% 
  mutate(class = ifelse(grepl(seq_name,pattern = "^AP"), "AMP", "non-AMP")) %>% 
  distinct(seq_name, .keep_all = TRUE) %>%
  add_column(dataset="Train") 
  
iAMP2L_test <- read_faa("data/amp_predictors/iAMP-2L/xiao_independent.fasta") %>% 
  mutate(class = ifelse(grepl(seq_name,pattern = "^AP"), "AMP", "non-AMP")) %>%
  add_column(dataset="Test") 
  
iAMP2L <- rbind(iAMP2L_train, iAMP2L_test) %>%
  mutate(length = nchar(seq_aa)) %>% add_column(database="iAMP-2L")
```
