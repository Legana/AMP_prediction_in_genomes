# This function calculates multiple performance metrics for a classification model with probabilities output

## Input to the function is a dataframe with "Neg" and Pos" named columns that contain the predicted probabilities from the
##   caret R package `predict(MODELNAME, TESTSET, type = "prob")` function. 
## The third required column is a column named "Label" which contains the actual class types ("Pos" or "Neg") for each prediction result
##   ( added from the test set, e.g. `add_column(Label = Testset$Label)` )

## Dependencies: this function depends on `tidyverse` for data wrangling and `precrec` for AUC calculations. 

## Output: Output is a dataframe with 9 columns

library(tidyverse)
library(precrec)

calculate_model_metrics <- function(df) {
  
  TP <- df %>% filter((Label=="Pos")) %>% filter(Pos >= 0.5) %>% n_distinct() %>% as.numeric()
  FP <- df %>% filter((Label=="Neg")) %>% filter(Pos >= 0.5) %>% n_distinct() %>% as.numeric()
  TN <- df %>% filter((Label=="Neg")) %>% filter(Pos < 0.5) %>% n_distinct() %>% as.numeric()
  FN <- df %>% filter((Label=="Pos")) %>% filter(Pos < 0.5) %>% n_distinct() %>% as.numeric()
  #as.numeric was necessary for the MCC calculation 
  #as otherwise it would result in a "NAs produced by integer overflow" error.
  
  Specificity <- round(TN / (TN + FP), digits = 3) #aka TNR
  Accuracy <- round((TP + TN) / (TP + TN + FP + FN), digits = 3)  # all correct / all | (misclassification = 1-accuracy)
  Recall <- round(TP / (TP + FN), digits = 3) # aka sensitivity, TPR 
  Precision <- round(TP/ (TP + FP), digits = 3) # positive predictive value
  FPR <- round(FP / (TN + FP), digits = 3) # false positive rate
  F1 <- round((2 * Precision * Recall) / (Precision + Recall), digits = 3) # F1 score
  MCC <- round(((TP*TN) - (FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)), digits = 3) # Matthews correlation coefficient
  
  df1 <- data.frame(FPR, Accuracy, Specificity, Recall, Precision, F1, MCC)
  
  df2 <- evalmod(scores = df$Pos, labels = df$Label, mode = "rocprc") %>% 
    precrec::auc() %>% 
    select(curvetypes, aucs) %>% 
    pivot_wider(names_from = curvetypes, values_from = aucs) %>% 
    rename(AUROC = "ROC", AUPRC = "PRC") %>% 
    round(digits = 3)
  
  cbind(df1, df2)
  
}