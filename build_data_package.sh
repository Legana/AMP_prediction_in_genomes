# ampir v_0.1.0 model related data

echo data/ampir_0.1.0_data/svm_Radial98_final.rds >> data_list
echo data/ampir_0.1.0_data/bg_tg98.rds >> data_list
echo data/ampir_0.1.0_data/features98.rds >> data_list
echo data/ampir_0.1.0_data/features98TrainNov19.rds >> data_list
echo data/ampir_0.1.0_data/features98TestNov19.rds >> data_list
echo data/ampir_0.1.0_data/ampir_prob_data.rds >> data_list

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


tar -zcvf data_amp_pred.tgz -T data_list

rm data_list