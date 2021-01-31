#function_to_calculate_cm_stats

calc_cm_metrics <- function(p_threshold, df) {
  
  
  # TP <- df %>% filter((actual=="Tg")) %>% filter(Tg > p_threshold) %>% n_distinct()
  # FP <- df %>% filter((actual=="Bg")) %>% filter(Tg > p_threshold) %>% n_distinct()
  # TN <- df %>% filter((actual=="Bg")) %>% filter(Tg < p_threshold) %>% n_distinct()
  # FN <- df %>% filter((actual=="Tg")) %>% filter(Tg < p_threshold) %>% n_distinct()
  
  # for when the column names are different ... 
  
  TP <- df %>% filter((Label=="Pos")) %>% filter(prob_AMP > p_threshold) %>% n_distinct()
  FP <- df %>% filter((Label=="Neg")) %>% filter(prob_AMP > p_threshold) %>% n_distinct()
  TN <- df %>% filter((Label=="Neg")) %>% filter(prob_AMP < p_threshold) %>% n_distinct()
  FN <- df %>% filter((Label=="Pos")) %>% filter(prob_AMP < p_threshold) %>% n_distinct()
  
  Specificity <- TN / (TN + FP) #aka TNR
  Recall <- TP / (TP + FN) # aka sensitivity, TPR
  Precision <- TP / (TP + FP)  # positive predictive value
  FPR <- FP / (TN + FP)
  #MCC <- (TP*TN) - (FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  
  cm <- c(TP, FP, TN, FN, Specificity, Recall, Precision, FPR, p_threshold)
  names(cm) <-c("TP", "FP", "TN", "FN", "Specificity", "Recall", "Precision", "FPR", "p_threshold") 
  cm
}