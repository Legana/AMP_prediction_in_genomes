Animal vs plants
================

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.4     ✓ dplyr   1.0.2
    ## ✓ tidyr   1.1.2     ✓ stringr 1.4.0
    ## ✓ readr   1.4.0     ✓ forcats 0.5.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    ## Loading required package: lattice

    ## 
    ## Attaching package: 'caret'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     lift

    ## Type 'citation("pROC")' for a citation.

    ## 
    ## Attaching package: 'pROC'

    ## The following object is masked from 'package:precrec':
    ## 
    ##     auc

    ## The following objects are masked from 'package:stats':
    ## 
    ##     cov, smooth, var

# Metazoa vs Viridiplantae

``` r
swissprot_amps_2 <- read_tsv("data/uniprot-keyword Antimicrobial+[KW-0929] -filtered-reviewed yes.tab") %>%
    rename(Taxonomic_lineage = `Taxonomic lineage (all)`) %>%
    rename(Entry_name = `Entry name`) %>%
    mutate(Taxonomic_lineage = case_when(
       str_detect(Taxonomic_lineage, "Metazoa") ~ "Animals",
       str_detect(Taxonomic_lineage, "Viridiplantae") ~ "Plants",
                                        TRUE ~ Taxonomic_lineage))
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   Entry = col_character(),
    ##   `Entry name` = col_character(),
    ##   Status = col_character(),
    ##   `Protein names` = col_character(),
    ##   Organism = col_character(),
    ##   Length = col_double(),
    ##   `Taxonomic lineage (all)` = col_character(),
    ##   `Taxonomic lineage (CLASS)` = col_character(),
    ##   Sequence = col_character(),
    ##   Peptide = col_character(),
    ##   Propeptide = col_character(),
    ##   `Signal peptide` = col_character()
    ## )

``` r
metazoa_pos <- filter(swissprot_amps_2, Taxonomic_lineage == "Animals" & Length < 500) %>%
  select(Entry_name, Sequence) %>%
  as.data.frame() %>%
  remove_nonstandard_aa() %>%
  calculate_features(min_len = 5) %>%
  mutate(Label = "Pos")

plants_pos <- filter(swissprot_amps_2, Taxonomic_lineage == "Plants" & Length < 500) %>% 
   select(Entry_name, Sequence) %>% 
   as.data.frame() %>% 
   remove_nonstandard_aa() %>% 
   calculate_features(min_len = 5) %>% 
   mutate(Label = "Pos")
```

``` r
set.seed(9)
metazoa_neg <- read_faa("data/uniprot-NOT+keyword Antimicrobial+[KW-0929] +taxonomy Metazoa+%5--.fasta.gz") %>%
   mutate(Length = nchar(seq_aa)) %>%
   filter(Length < 500) %>% 
   as.data.frame() %>%
   remove_nonstandard_aa() %>%
   slice_sample(n = 10*nrow(metazoa_pos)) %>%
   calculate_features(min_len = 5) %>%
   mutate(Label = "Neg")

plants_neg <- read_faa("data/uniprot-NOT+keyword Antimicrobial+[KW-0929] +taxonomy Viridiplan--.fasta.gz") %>%
  mutate(Length = nchar(seq_aa)) %>%
  filter(Length < 500) %>% 
  as.data.frame() %>%
  remove_nonstandard_aa() %>%
  slice_sample(n = 10*nrow(plants_pos)) %>%
  calculate_features(min_len = 5) %>%
  mutate(Label = "Neg")
```

``` r
metazoa_feat <- rbind(metazoa_pos, metazoa_neg) %>% mutate(Label = as.factor(Label))
plants_feat <- rbind(plants_pos, plants_neg) %>% mutate(Label = as.factor(Label))
```

``` r
trainIndex <-createDataPartition(y = metazoa_feat$Label, p=.8, list = FALSE)
features_metazoaTrain <- metazoa_feat[trainIndex,]
features_metazoaTest <- metazoa_feat[-trainIndex,]
```

Metazoa total has 2381 AMPs and 23810 non-AMPs, training data consists
of 1905 AMPs and 19048 non-AMPs, test set consists of 475 AMPs and 4762
nonAMPs Plants has 515 AMPs and 5,150 non-AMPs

``` r
model_weights <- ifelse(features_metazoaTrain$Label == "Pos",
                        (1/table(features_metazoaTrain$Label)[1]) * 0.5,
                        (1/table(features_metazoaTrain$Label)[2]) * 0.5)

trctrl_prob <- trainControl(method = "repeatedcv", number = 10, repeats = 3, classProbs = TRUE)

metazoa_svm <- train(Label~.,
                       data = features_metazoaTrain[,c(2:28,46)],
                       method="svmRadial",
                       trControl = trctrl_prob,
                       preProcess = c("center", "scale"),
                       weights = model_weights,
                       metric = "Kappa",
                       tuneLength = 6)
```

``` r
metazoa_predict_and_actual <- predict(metazoa_svm, features_metazoaTest, type = "prob") %>%
  add_column(Label = features_metazoaTest$Label)


plants_predict_and_actual <- predict(metazoa_svm, plants_feat, type = "prob") %>% add_column(Label = plants_feat$Label)
```

## Calculating performance of both models on test sets

``` r
calculate_model_metrics <- function(df) {

  TP <- df %>% filter((Label=="Pos")) %>% filter(Pos >= 0.5) %>% n_distinct() %>% as.numeric()
  FP <- df %>% filter((Label=="Neg")) %>% filter(Pos >= 0.5) %>% n_distinct() %>% as.numeric()
  TN <- df %>% filter((Label=="Neg")) %>% filter(Pos < 0.5) %>% n_distinct() %>% as.numeric()
  FN <- df %>% filter((Label=="Pos")) %>% filter(Pos < 0.5) %>% n_distinct() %>% as.numeric()
  #as.numeric was necessary for the MCC calculation 
  #as otherwise it would result in a "NAs produced by integer overflow" error.
  
  Specificity <- round(TN / (TN + FP), digits = 3) 
  Accuracy <- round((TP + TN) / (TP + TN + FP + FN), digits = 3) 
  Recall <- round(TP / (TP + FN), digits = 3)
  Precision <- round(TP/ (TP + FP), digits = 3) 
  FPR <- round(FP / (TN + FP), digits = 3) 
  F1 <- round((2 * Precision * Recall) / (Precision + Recall), digits = 3) 
  MCC <- round(((TP*TN) - (FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)), digits = 3) 
  
  df1 <- data.frame(FPR, Accuracy, Specificity, Recall, Precision, F1, MCC)
  
  df2 <- evalmod(scores = df$Pos, labels = df$Label, mode = "rocprc") %>% precrec::auc() %>% select(curvetypes, aucs) %>% pivot_wider(names_from = curvetypes, values_from = aucs) %>% rename(AUROC = "ROC", AUPRC = "PRC") %>% round(digits = 3)
  
  cbind(df1, df2)
  
}
```

``` r
metazoa_pred1 <- calculate_model_metrics(metazoa_predict_and_actual ) %>% mutate(Model = "Metazoa1") %>% mutate(Data = "Metazoa1")

plants_pred1 <- calculate_model_metrics(plants_predict_and_actual) %>% mutate(Model = "Metazoa1") %>% mutate(Data = "Plants1")

metazoa_models_performance <- rbind(metazoa_pred1, plants_pred1 ) 
```

### second metazoa model

Similar to the deuterostome models.. a second metazoa model was created
where the metazoa training data equals the plant dataset (515 AMPs and
5,150 non-AMPs)

``` r
trainIndex2 <-createDataPartition(y = metazoa_feat$Label, p=.21626, list = FALSE)
features_metazoaTrain2 <- metazoa_feat[trainIndex2,]
features_metazoaTest2 <- metazoa_feat[-trainIndex2,] 

features_metazoaTest2 <-  features_metazoaTest2 %>% slice(c(1:515, 15376:20525)) #subset test set to same size as training data 
```

``` r
model_weights2 <- ifelse(features_metazoaTrain2$Label == "Pos",
                        (1/table(features_metazoaTrain2$Label)[1]) * 0.5,
                        (1/table(features_metazoaTrain2$Label)[2]) * 0.5)

trctrl_prob2 <- trainControl(method = "repeatedcv", number = 10, repeats = 3, classProbs = TRUE)

metazoa_svm2 <- train(Label~.,
                       data = features_metazoaTrain2[,c(2:28,46)],
                       method="svmRadial",
                       trControl = trctrl_prob2,
                       preProcess = c("center", "scale"),
                       weights = model_weights2,
                       metric = "Kappa",
                       tuneLength = 6)
```

``` r
metazoa_svm2 <- readRDS("cache/metazoa_svm2.rds")
```

``` r
metazoa_predict_and_actual2 <- predict(metazoa_svm2, features_metazoaTest2, type = "prob") %>%
  add_column(Label = features_metazoaTest2$Label)


plants_predict_and_actual2 <- predict(metazoa_svm2, plants_feat, type = "prob") %>% add_column(Label = plants_feat$Label)
```

## Calculating performance of both models on test sets

``` r
metazoa_pred1 <- calculate_model_metrics(metazoa_predict_and_actual ) %>% mutate(Model = "Metazoa1") %>% mutate(Data = "Metazoa1")

plants_pred1 <- calculate_model_metrics(plants_predict_and_actual) %>% mutate(Model = "Metazoa1") %>% mutate(Data = "Plants1")

metazoa_models_performance <- rbind(metazoa_pred1, plants_pred1 ) 


metazoa_pred2 <- calculate_model_metrics(metazoa_predict_and_actual2) %>% mutate(Model = "Metazoa2") %>% mutate(Data = "Metazoa2")

plants_pred2 <- calculate_model_metrics(plants_predict_and_actual2) %>% mutate(Model = "Metazoa2") %>% mutate(Data = "Plants1")

metazoa_models_performance <- rbind(metazoa_models_performance, metazoa_pred2, plants_pred2)
```

**Table 5.5:** Performance of both deuterostome models on all test sets

| Accuracy | Specificity | Recall | Precision |    F1 |   MCC | AUROC | AUPRC | Model    | Data     |
|---------:|------------:|-------:|----------:|------:|------:|------:|------:|:---------|:---------|
|    0.963 |       0.989 |  0.708 |     0.866 | 0.779 | 0.764 | 0.961 | 0.846 | Metazoa1 | Metazoa1 |
|    0.896 |       0.980 |  0.095 |     0.338 | 0.148 | 0.137 | 0.717 | 0.221 | Metazoa1 | Plants1  |
|    0.949 |       0.988 |  0.561 |     0.821 | 0.667 | 0.654 | 0.924 | 0.735 | Metazoa2 | Metazoa2 |
|    0.886 |       0.975 |  0.047 |     0.161 | 0.073 | 0.038 | 0.677 | 0.166 | Metazoa2 | Plants1  |

``` r
m1_roc <- roc(metazoa_predict_and_actual$Label, metazoa_predict_and_actual$Pos)
```

    ## Setting levels: control = Neg, case = Pos

    ## Setting direction: controls < cases

``` r
m2_roc <- roc(metazoa_predict_and_actual2$Label, metazoa_predict_and_actual2$Pos)
```

    ## Setting levels: control = Neg, case = Pos
    ## Setting direction: controls < cases

``` r
p1_roc <- roc(plants_predict_and_actual$Label, plants_predict_and_actual$Pos)
```

    ## Setting levels: control = Neg, case = Pos
    ## Setting direction: controls < cases

``` r
p2_roc <- roc(plants_predict_and_actual2$Label, plants_predict_and_actual2$Pos)
```

    ## Setting levels: control = Neg, case = Pos
    ## Setting direction: controls < cases

``` r
roc.test(m1_roc, m2_roc, method = "delong")
```

    ## 
    ##  DeLong's test for two ROC curves
    ## 
    ## data:  m1_roc and m2_roc
    ## D = 4.0015, df = 10064, p-value = 6.341e-05
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.9610771   0.9236901

``` r
roc.test(m1_roc, p1_roc, method = "delong")
```

    ## 
    ##  DeLong's test for two ROC curves
    ## 
    ## data:  m1_roc and p1_roc
    ## D = 17.228, df = 7528, p-value < 2.2e-16
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.9610771   0.7167469

``` r
roc.test(m2_roc, p2_roc, method = "delong")
```

    ## 
    ##  DeLong's test for two correlated ROC curves
    ## 
    ## data:  m2_roc and p2_roc
    ## Z = 16.377, p-value < 2.2e-16
    ## alternative hypothesis: true difference in AUC is not equal to 0
    ## sample estimates:
    ## AUC of roc1 AUC of roc2 
    ##   0.9236901   0.6772224

Model and data 1 vs. model and data 2 \| P-value \|

Metazoa\_svm1 on metazoa1 data vs. Metazoa\_svm2 on metazoa2 data \| 6.3
x 10<sup>-05 \|  
Metazoa\_svm1 on metazoa1 data vs. Metazoa\_svm1 on plant data \| &lt;
2.2 x 10<sup>-16</sup> \|  
Metazoa\_svm2 on metatazoa2 data vs. Metazoa\_svm2 on plant data \| &lt;
2.2 x 10<sup>-16</sup> \|

## PCA and tSNE of metazoa and plants AMPs vs nonAMPs

``` r
set.seed(3)

metaplant_feat <- metazoa_feat %>% mutate(group = "Animals") %>% rbind(plants_feat %>% mutate(group = "Plants"))
```

**PCA**

``` r
pca_metaplant_feat <- metaplant_feat %>% 
  column_to_rownames("seq_name") %>% 
  select(Amphiphilicity:Xc2.lambda.4) %>% 
  prcomp(scale. = TRUE)


percentage_meta <- round(pca_metaplant_feat$sdev^2 / sum(pca_metaplant_feat$sdev^2) * 100, 1)
percentage_meta <- paste(colnames(pca_metaplant_feat$x),"(",paste(as.character(percentage_meta), "%",")", sep = ""))

pca_metaplant_feat.x <- pca_metaplant_feat$x %>% 
   as.data.frame() %>% 
   rownames_to_column("seq_name") %>% 
   left_join(metaplant_feat %>% select(seq_name, Label, group), by = "seq_name") %>%
   mutate(taxgroup = as.factor(group))
```

**tSNE**

``` r
pca_metaplant_feat_plot <- ggplot(pca_metaplant_feat.x, aes(x = PC1, y = PC2)) +
   geom_point(aes(colour = taxgroup, shape = taxgroup), size = 2, alpha = 0.9) +
   facet_wrap(~factor(Label, levels = c("Pos", "Neg"), labels = c("AMPs", "non-AMPs"))) +
   labs(x = percentage_meta[1], y = percentage_meta[2], colour = "", shape = "") +
   scale_shape_manual(labels = c("Animals" , "Plants"), values=c(16, 5)) +
   scale_colour_manual(labels = c("Animals" , "Plants"), values = c("blueviolet", "forestgreen")) +
   theme_classic() +
   theme(legend.position = "none",
         strip.background = element_blank(),
         strip.text = element_text(size = 10, face = "bold")) 


pca1_metaplant_feat_plot <- ggplot(pca_metaplant_feat.x, aes(x = PC1)) +
   stat_density(aes(colour = group), geom = "line", position = "identity") +
   facet_wrap(~factor(Label, levels = c("Pos", "Neg"), labels = c("AMPs", "non-AMPs"))) +
   labs(x = percentage_meta[1], y = "Density", colour = "") +
   scale_colour_manual(values = c("blueviolet", "forestgreen")) +
   theme_classic() +
   theme(legend.position = "none",
         strip.background = element_blank(),
         strip.text = element_text(size = 10, face = "bold"))

tsne_metaplant_feat_plot <- ggplot(metaplant_tsne_annotated, aes(x = tSNE_1, y = tSNE_2)) +
   geom_point(aes(colour = group, shape = group), size = 2, alpha = 0.9) +
   facet_wrap(~factor(Label, levels = c("Pos", "Neg"), labels = c("AMPs", "non-AMPs"))) +
   labs(x = "tSNE 1", y = "tSNE 2", colour = "", shape = "") +
   scale_shape_manual(labels = c("Animals" , "Plants"), values=c(16, 5)) +
   scale_colour_manual(labels = c("Animals" , "Plants"), values = c("blueviolet", "forestgreen")) +
   theme_classic() +
   theme(legend.position = "bottom",
         strip.background = element_blank(),
         strip.text = element_text(size = 10, face = "bold")) 

(pca_metaplant_feat_plot | pca1_metaplant_feat_plot) / tsne_metaplant_feat_plot + plot_annotation(tag_levels = 'A')
```

![](06_animals_vs_plants_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->
