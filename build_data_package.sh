# ampir v_0.1.0 model related data

echo data/ampir_0.1.0_data/svm_Radial98_final.rds >> data_list
echo data/ampir_0.1.0_data/bg_tg98.rds >> data_list
echo data/ampir_0.1.0_data/features98.rds >> data_list
echo data/ampir_0.1.0_data/features98TrainNov19.rds >> data_list
echo data/ampir_0.1.0_data/features98TestNov19.rds >> data_list
echo data/ampir_0.1.0_data/ampir_prob_data.rds >> data_list

# ampir v_1 related data

echo data/ampir_v1/tuned_precursor_imbal_nobench.rds >> data_list


# Data from AMP predictors 

echo data/amp_predictors/amPEP/M_model_train_nonAMP_sequence.fasta >> data_list
echo data/amp_predictors/amPEP/M_model_train_AMP_sequence.fasta >> data_list
echo data/amp_predictors/deepamPEP30/train_ne.fasta >> data_list
echo data/amp_predictors/deepamPEP30/test_ne.fasta >> data_list
echo data/amp_predictors/deepamPEP30/train_po.fasta >> data_list
echo data/amp_predictors/deepamPEP30/test_po.fasta >> data_list
echo data/amp_predictors/amPEPpy/M_model_train_nonAMP_sequence.numbered.proplen.subsample.fasta >> data_list
echo data/amp_predictors/amPEPpy/M_model_train_AMP_sequence.numbered.fasta >> data_list
echo data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/AMP.tr.fa >> data_list
echo data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/AMP.eval.fa >> data_list
echo data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/DECOY.tr.fa >> data_list
echo data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/DECOY.eval.fa >> data_list
echo data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/AMP.te.fa >> data_list
echo data/amp_predictors/AMP_Scan2_OrigPaper_Dataset/DECOY.te.fa >> data_list
echo data/amp_predictors/iAMP-2L/xiao_benchmark.fasta >> data_list
echo data/amp_predictors/iAMP-2L/xiao_independent.fasta >> data_list
echo data/amp_predictors/AMPlify/AMP_test_20190414.fa >> data_list
echo data/amp_predictors/AMPlify/AMP_train_20190414.fa >> data_list
echo data/amp_predictors/AMPlify/non_AMP_test_20190414.fa >> data_list
echo data/amp_predictors/AMPlify/non_AMP_train_20190414.fa >> data_list
echo data/amp_predictors/AmpGram/benchmark.fasta >> data_list

# Positive and negative datasets from the Swiss-Prot database (accessed January 2021)
echo data/uniprot-keyword Antimicrobial+[KW-0929] -filtered-reviewed yes.tab >> data_list
echo data/uniprot-taxonomy Deuterostomia+[33511] +reviewed yes+NOT+keyword--.fasta >> data_list
echo data/uniprot-taxonomy Protostomia+[33317] +reviewed yes+NOT+keyword--.fasta >> data_list

# Deuterostome and protostome model data

echo cache/deuterostome_neg.rds >> data_list
echo cache/protostome_neg.rds >> data_list
echo cache/features_deutTrain.rds >> data_list
echo cache/features_deutTrain2.rds >> data_list
echo cache/features_deutTest.rds >> data_list
echo cache/features_deutTest2.rds >> data_list
echo cache/features_protostome.rds >> data_list
echo cache/deut_svm.rds >> data_list
echo cache/deut_svm2.rds >> data_list
echo cache/deutprottsne.rds >> data_list

# Arabidopsis thaliana and Homo sapiens proteomes without and with standard AA

echo data/proteomes/arabidopsis-proteomeUP000006548.fasta >> data_list
echo data/proteomes/human-proteomeUP000005640.fasta >> data_list
echo data/proteomes/uniprot-proteome UP000005640.tab >> data_list
echo data/proteomes/uniprot-proteome UP000006548.tab >> data_list
echo cache/arab_proteome_standardaa.fasta >> data_list
echo cache/homo_proteome_standardaa.fasta >> data_list

# Benchmark prediction results

echo data/prediction_results/ampir/ampir_refproteomes_predictions.rds >> data_list
echo data/prediction_results/ampscanner_v2/arab1_1611719550167_Prediction_Summary.csv >> data_list
echo data/prediction_results/ampscanner_v2/arab2_1611719690212_Prediction_Summary.csv >> data_list
echo data/prediction_results/ampscanner_v2/homo1_1611718866504_Prediction_Summary.csv >> data_list
echo data/prediction_results/ampscanner_v2/homo2_1611719132616_Prediction_Summary.csv >> data_list
echo data/prediction_results/ampscanner_v2/homo3_1611719332326_Prediction_Summary.csv >> data_list
echo data/prediction_results/ampep/arab_proteome_standardaa_ampep.txt >> data_list
echo data/prediction_results/ampep/homo_proteome_standardaa_ampep.txt >> data_list



tar -zcvf data_amp_pred.tgz -T data_list

rm data_list