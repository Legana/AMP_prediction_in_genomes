Training datasets in other AMP predictors
================

# Introduction

The effectiveness of supervised learning methods for predictive modeling
is highly affected by the data that were used to train the model with.
Supervised learning is a part of machine learning which finds patterns
between data characteristics and labels (positive and negative) that are
assigned to the training data used to facilitate learning of the model.
Once the model has been trained, it can then predict which label fits
best from new data it is used or tested on. For antimicrobial peptide
(AMP) prediction, positive data refers to AMP sequences and negative
data generally refers to sequences that presumably do not have
antimicrobial activity (non-AMPs).

# AMP predictor data

The training and test sets of six recent AMP predictors were examined to
assess their effectiveness in genome-wide AMP prediction.

**Table 2.1:** Summary table of the number of positive and negative
sequences present in the training and test set in six AMP predictors

| AMP predictor       | Train - AMPs | Train - non-AMPs | Test - AMPs | Test - non-AMPs |
|:--------------------|:------------:|:----------------:|:-----------:|:---------------:|
| iAMP-2L             |     897      |      2,405       |     920     |       920       |
| amPEP               |    3,268     |     166,791      |    *NS*     |      *NS*       |
| Deep-ampEP30        |    1,529     |      1,529       |     94      |       94        |
| amPEPpy             |    3,268     |      3,268       |    *NS*     |      *NS*       |
| AMP Scanner v2      |    1,066     |      1,066       |     712     |       712       |
| AMPlify             |    3,338     |      3,338       |     835     |       835       |
| AmpGram             |    2,216     |      2,216       |     247     |       247       |
| ampir\_precursor v1 |    1,187     |      11,864      |     296     |      2,966      |
| ampir\_mature v1    |    2,586     |      2,657       |     646     |       664       |

\*NS: not specified

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
  
iAMP2L_data <- rbind(iAMP2L_train, iAMP2L_test) %>%
  mutate(length = nchar(seq_aa)) %>% add_column(predictor="iAMP-2L")
```

**AmPEP Training Data**

The AmPEP 2018 AMP predictor provides its training data available
directly for download from
<https://cbbio.cis.um.edu.mo/software/AmPEP/>. The final training
dataset used by amPEP is a large dataset of 166,791 non-AMP sequences
and 3,268 AMPs. amPEP used the Xiao et al. 2013 dataset from the iAMP-2L
predictor as a test set (see above).

``` r
ampep_data <- read_faa("data/amp_predictors/amPEP/M_model_train_nonAMP_sequence.fasta") %>% add_column(class="non-AMP") %>% 
  rbind(read_faa("data/amp_predictors/amPEP/M_model_train_AMP_sequence.fasta") %>% add_column(class="AMP")) %>% 
  add_column(dataset = "Train") %>%
  mutate(length = nchar(seq_aa)) %>% add_column(predictor="AmPEP") %>%
  mutate(seq_name = paste0("amPEP_trainset_neg", 1:n()))
```

AmPEP was redesigned in 2020 as Deep-AmPEP30 to focus on short AMPs (
&lt; 30 amino acids) by [Yan et
al](https://doi.org/10.1016/j.omtn.2020.05.006) and its training and
test data is available
[here](https://cbbio.online/AxPEP/?action=dataset). Deep-AmPEP30’s
training set consists of 1,529 AMPs and non-AMPs and their test set
consists of 94 AMPs and non-AMPs.

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
  add_column(dataset = "Train") %>% mutate(length = nchar(seq_aa)) %>% add_column(predictor="AmPEPpy")
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
    add_column(predictor="AMP Scanner v2") %>%
  relocate(class, .before = dataset)
```

**AMPlify data**

AMPlify used multiple attention mechanisms and ensemble deep learning to
create an AMP prediction model ([Li et
al. 2020](https://doi.org/10.1101/2020.06.16.155705)). Similar to most
AMP predictors, a proportion of their AMP dataset originated from the
general [Antimicrobial Peptide Database](http://aps.unmc.edu/AP).
However, they also used the [Database of Anuran Defense
Peptides](http://split4.pmfst.hr/dadp/0) which focuses on AMPs from
frogs and toads. Their negative dataset originated from Swiss-Prot and,
like most AMP predictors, they excluded any proteins that had
annotations which referred to potential antimicrobial activity. However,
like ampir, they retain the secreted proteins. Their data is available
from the [AMPlify’s software GitHub
page](https://github.com/bcgsc/AMPlify). Their training set consists of
3,338 AMPs and 3,338 non-AMPs and their test set consists of 835 AMPs
and 835 non-AMPs.

``` r
amplify_data <- read_faa("data/amp_predictors/AMPlify/AMP_train_20190414.fa") %>%
   rbind(read_faa("data/amp_predictors/AMPlify/non_AMP_train_20190414.fa")) %>%
  add_column(dataset = "Train") %>%
  rbind(read_faa("data/amp_predictors/AMPlify/AMP_test_20190414.fa") %>%
  rbind(read_faa("data/amp_predictors/AMPlify/non_AMP_test_20190414.fa")) %>%
          add_column(dataset = "Test")) %>%
  mutate(class = case_when(str_detect(seq_name, "^trAMP|^teAMP") ~ "AMP", TRUE ~ "non-AMP")) %>%
  mutate(length = nchar(seq_aa)) %>% add_column(predictor="AMPlify") %>% relocate(class, .before = dataset)
```

AmpGram was created in 2020 and in addition to the standard AMPs, it
also focuses on predicting longer proteins that contain antimicrobial
activity, such as the milk protein, lactoferrin, and on non-AMP proteins
which produce antimicrobial proteolysis products, such as human
thrombin. AmpGram uses amino acid motifs and the random forest model to
classify AMPs [Burdukiewicz et
al. 2020](https://doi.org/10.3390/ijms21124310). AmpGram’s analysis
details and data can be obtained from [AmpGram’s analysis
repository](https://github.com/michbur/AmpGram-analysis). Their training
data is not specified as a file but their benchmark data can be found in
the `benchmark.fasta` file. This test set contains 247 AMP and non-AMP
sequences. AmpGram also used the dataset from [Gabere and Noble
2017](https://doi.org/10.1093/bioinformatics/btx081) which used AMPs
from the [DAMPD](https://dx.doi.org/10.1093%2Fnar%2Fgkr1063) and
[APD3](https://doi.org/10.1093/nar/gkv1278) AMP databases and non-AMPs
from UniProt.

``` r
ampgram_test_data <- read_faa("data/amp_predictors/AmpGram/benchmark.fasta") %>% mutate(class = case_when(str_detect(seq_name, "^dbAMP") ~ "AMP", TRUE ~ "non-AMP")) %>% add_column(dataset = "Test") %>% mutate(length = nchar(seq_aa)) %>% add_column(predictor = "AmpGram")
```

**ampir**

ampir was divided in two different models, precursor, which focuses on
longer sequences (between 60-300) and mature, which only contains short
sequences (between 10-50)

``` r
ampir_prec_feats_train <- readRDS("data/ampir_v1/featuresTrain_precursor_imbal.rds") %>% mutate(dataset = "Train")
ampir_prec_feats_test <- readRDS("data/ampir_v1/featuresTest_precursor_imbal.rds") %>% mutate(dataset = "Test")
ampir_prec_feats <- rbind(ampir_prec_feats_train, ampir_prec_feats_test) %>% mutate(predictor = "ampir_precursor")

ampir_mat_feats_train <- readRDS("data/ampir_v1/featuresTrain_mature.rds") %>% mutate(dataset = "Train")
ampir_mat_feats_test <- readRDS("data/ampir_v1/featuresTest_mature.rds") %>% mutate(dataset = "Test")
ampir_mat_feats <- rbind(ampir_mat_feats_train, ampir_mat_feats_test) %>% mutate(predictor = "ampir_mature")

ampir_feats <- rbind(ampir_prec_feats, ampir_mat_feats) %>% mutate(class = ifelse(Label == "Tg", "AMP","non-AMP"))

ampir_data <- ampir_feats %>% select(seq_name, seq_aa, class, dataset, length, predictor)
```

``` r
all_predictor_data <- rbind(iAMP2L_data, ampep_data, deep_ampep_data, ampeppy_data, ampscan_data, amplify_data, ampgram_test_data, ampir_data)

all_predictor_data <- all_predictor_data %>% filter(length >10)
```

``` r
all_predictor_data_wcounts <- all_predictor_data %>%
                                group_by(class, dataset, length, predictor) %>%
                                summarise(number = n()) %>%
    mutate(predictor = factor(predictor, levels = c("iAMP-2L", "AMP Scanner v2", "AmPEP", "AmPEPpy", "deep_AmPEP", "AMPlify", "AmpGram", "ampir_precursor", "ampir_mature")))
```

## PCA of AMP and non-AMP features on models and proteomes

Read in *Homo sapiens* and *Arabidopsis thaliana* proteomes and keep
only sequences longer than 10 amino acids and with standard amino acids
to match predictor data used for the previous PCA.

``` r
reference_proteomes <- read_tsv("data/proteomes/uniprot-proteome UP000005640.tab") %>%
  rbind(read_tsv("data/proteomes/uniprot-proteome UP000006548.tab")) %>%
  mutate(Label = case_when(str_detect(`Keyword ID`, "KW-0929") ~ "Pos", TRUE ~ "Neg")) %>%
  filter(Length >10) %>%
  filter(Length <3000) %>%
  filter(grepl(Sequence, pattern='^[ARNDCEQGHILKMFPSTWYV]+$'))
```

*Calculating features with ampir*

``` r
all_predictor_data_feats <- all_predictor_data %>% calculate_features(min_len = 10)

reference_proteomes_seqnames <- reference_proteomes %>%
  select(`Entry name`, Sequence) %>%
  as.data.frame()

reference_proteomes_feats <- reference_proteomes_seqnames %>% calculate_features()
```

*add class and names and combine model and proteome features*

``` r
all_predictor_data_feats <- all_predictor_data_feats %>% mutate(class = all_predictor_data$class) %>% mutate(name = all_predictor_data$predictor) %>% mutate(length = all_predictor_data$length)

reference_proteomes_feats <- reference_proteomes_feats %>% mutate(class = ifelse(reference_proteomes$Label == "Pos", "AMP","non-AMP")) %>% mutate(name = case_when(str_detect(reference_proteomes$Organism, "Homo") ~ "Homo sapiens", TRUE ~ "Arabidopsis thaliana")) %>% mutate(length = reference_proteomes$Length)

predictor_and_proteome_feats <- rbind(all_predictor_data_feats, reference_proteomes_feats) %>% mutate(name = factor(name, levels = c("iAMP-2L", "AMP Scanner v2", "AmPEP", "AmPEPpy", "deep_AmPEP", "AMPlify", "AmpGram", "ampir_precursor", "ampir_mature", "Homo sapiens", "Arabidopsis thaliana")))


predictor_and_proteome_counts <- predictor_and_proteome_feats %>%
                            group_by(class, length, name) %>%
                            summarise(number = n())
```

``` r
pca_features_pp <- predictor_and_proteome_feats %>% 
   select(c(Amphiphilicity:Xc2.lambda.2)) %>%
   prcomp(scale. = TRUE)

pca_values_pp <- pca_features_pp$x %>% 
   as.data.frame() %>%
   mutate(seq_name = predictor_and_proteome_feats$seq_name) %>%
   left_join(predictor_and_proteome_feats, by = "seq_name")
```

``` r
pca_prot_percentages <- round(pca_features_pp$sdev^2 / sum(pca_features_pp$sdev^2) * 100, 2)
pca_prot_percentages <- paste(colnames(pca_features_pp$x),"(",paste(as.character(pca_prot_percentages), "%",")", sep = ""))

pca_values_pp_amps <- filter(pca_values_pp, class == "AMP")
pca_values_pp_nonamps <- filter(pca_values_pp, class == "non-AMP")

model_seqlength <- ggplot(predictor_and_proteome_counts, aes(x = length, y = number)) +
  geom_col(aes(fill = factor(class, levels = c("non-AMP", "AMP")))) +
  facet_grid(name ~ . , scales = "free_y") +
  labs(x = "Sequence length", y = "Number of sequences", fill = "") +
  xlim(0,300) +
  theme(legend.position = "none",
        strip.text.y.right = element_text(angle = 0, hjust = 0),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(face = "italic")) +
  scale_fill_manual(values = c("AMP" = "darkblue", "non-AMP" = "brown4")) 

pca_models <-  ggplot(pca_values_pp) +
   geom_point(data = pca_values_pp_nonamps, aes(x = PC1, y = PC2, colour = class, shape = class), size = 0.7) +
   geom_point(data = pca_values_pp_amps, aes(x = PC1, y = PC2, colour = class, shape = class), size = 0.7) +
   facet_grid(name ~., scales = "free_y") +
   labs(x = pca_prot_percentages[1], y = pca_prot_percentages[2], shape = "", colour = "") +
   theme(legend.position = "bottom",
        #strip.text.y.right = element_text(angle = 0, hjust = 0),
        strip.text.y.right = element_blank(),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(face = "italic")) +
        scale_colour_manual(values = c("darkblue", "brown4")) +
        scale_shape_manual(values = c(1, 3)) +
   guides(colour = guide_legend(override.aes = list(size=1)))


pca1_models <- ggplot(pca_values_pp, aes(x = PC1)) +
   stat_density(aes(colour = class), geom = "line", position = "identity") +
   facet_grid(name ~., scales = "free_y") +
   labs(x = pca_prot_percentages[1], y = "Density", colour = "") +
   theme(legend.position = "bottom",
        strip.text.y.right = element_text(angle = 0, hjust = 0),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(face = "italic")) +
scale_colour_manual(values = c("brown4", "darkblue")) +
   guides(colour = guide_legend(override.aes = list(size=1)))

model_seqlength /  (pca_models | pca1_models) + plot_annotation(tag_levels = 'A')
```

![](02_amp_models_trainingdata_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

**Figure 2.1:** A) Sequence length, B) Scatterplot of the first two
principal components (PC), C) Density plot of the first PC of all AMP
and non-AMP sequences used in various AMP prediction models and in the
proteomes of *Homo sapiens* and *Arabidopsis thaliana*.
