library(AmpGram)
library(AmpGramModel)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE)

print("args is:")
args


input_fasta <- args[1]

print("inputfasta is:")
input_fasta

outputname <- paste(input_fasta, ".rds" , sep="")


arab_ampgram <- read_txt(input_fasta)

arab_ampgram <- arab_ampgram

arab_ampgram_pred <- predict(AmpGram_model, arab_ampgram)

print("prediction done!")

saveRDS(arab_ampgram_pred, outputname)

print("prediction saved")