
#' This function first tranform FPKM data into log2-scale transformation, followed by quantile normalization, where training data values used as target matrix both for training and test data for quantile  normalization 
#' @param train_data :args1 - training data (TSV) (Patients data with clinical and gene expression, where samples are in rows and features/genes are in columns)
#' @param test_data_data :args2 - test data (TSV) (Patients data with clinical and gene expression, where samples are in rows and features/genes are in columns)
#' @param col_num :args3 - column number in data at where clinical info ends 
#' @param train_clin_data: args4- name of training data output which stores only clinical information 
#' @param test_clin_data: args5 - name of test data output which stores only clinical information 
#' @param train_Normalized_data_clin_data: args6 -  output file name - outfile file containing training clinical data and normalized data  (which is log-scaled followed by quantile normalized data).
#' @param test_Normalized_data_clin_data: args7- output filename - outfile file containing test clinical data and normalized data  (which is log-scaled followed by quantile normalized data).
#' @examples
#' SurvPredPipe::train_test_normalization_f(train_data="train_FPKM.txt",test_data="test_FPKM.txt", col_num=21, train_clin_data="Train_Clin.txt", test_clin_data="TestClin.txt", train_Normalized_data_clin_data="Train_Norm_data.txt", test_Normalized_data_clin_data="Test_Norm_data.txt")
#' Usage: train_test_normalization_f(train_data,  test_data, col_num, train_clin_data, test_clin_data, train_Normalized_data_clin_data, test_Normalized_data_clin_data)

#' @export
#setwd: /Users/kaurh8/Documents/GDC_TCGA_Biolinks/GDC_All_samples_Data/TCGA-LGG/LGG_Survival_package/Final_package_functions

train_test_normalization_f <- function(train_data,  test_data, col_num, train_clin_data, test_clin_data, train_Normalized_data_clin_data, test_Normalized_data_clin_data) {
  
  
  # Check if any input variable is empty
  if (length(train_data) == 0 ||  length(test_data) == 0|| length(col_num) == 0 ||  length(train_clin_data) == 0 ||  length(test_clin_data) == 0 ||  length(train_Normalized_data_clin_data) == 0 ||  length(test_Normalized_data_clin_data) == 0) {
    stop("Error: Empty input variable detected.")
  }
  
  # Check if any input variable is missing
  if (any(is.na(train_data)) || any(is.na(test_data))  || any(is.na(col_num))) {
    stop("Error: Missing values in input variables.")
  }
  
  #training and test data contains clinical data in first 21 columns and rest columns repression gene expression values 
  
 # training <- read.table("Training_set_FPKM_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  #testing <- read.table("Test_set_FPKM_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)  
  
  
  training <- read.table(train_data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  testing <- read.table(test_data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)  
  
  
  #head(training [1:22],2)
  
  n<- col_num - 1
  #extract clinical data  
  #tr_clin <- training[1:20]
  #te_clin <- testing[1:20]
  
  tr_clin <- training[1:n]
  te_clin <- testing[1:n]
  
  # write Clinical data into files
  #write.table(cbind("ID"=rownames(Train_norm_data), Train_norm_data),file="Log_scaled_and_QN_train_exp_data.txt",sep="\t",quote=F, row.names=FALSE)
  write.table(cbind("ID"=rownames(tr_clin), tr_clin),file=train_clin_data,sep="\t",quote=F, row.names=FALSE)
  
  #write.table(cbind("ID"=rownames(Test_norm_data), Test_norm_data),file="Log_scaled_and_QN_test_exp_data.txt",sep="\t",quote=F, row.names=FALSE)
  write.table(cbind("ID"=rownames(te_clin), te_clin),file=test_clin_data,sep="\t",quote=F, row.names=FALSE)
  
  
  ##Extract Expression dat
  tr_exp <- training[col_num:ncol(training)]
  te_exp<- testing[col_num:ncol(testing)]
  #tr_exp <- training[21:ncol(training)]
  #te_exp<- testing[21:ncol(testing)]
  
  head(tr_exp[1:5],2)
  
  ##log scale transformation
  
  #library(MASS);
  tr_exp1 <-  (tr_exp+1);
  tr_log_mat = round(log2(tr_exp1 ),3);
  
  
  te_exp1 <-  (te_exp+1);
  te_log_mat = round(log2(te_exp1 ),3);
  
  #write.table(cbind("ID"=rownames(log_mat), log_mat),file=args[2],sep="\t",quote=F, row.names=F)
  
  
  
  #quantile normalization
  ####Usage Rscript quantile_2_datamatrices_according_1.R train_mat.csv test_mat.csv
  #args <- commandArgs(TRUE);
  
  
  #samples in columns and genes in the rows
  #transpose train data
  ref<- as.data.frame(t(tr_exp))
  head(ref[1:10],2)
  
  #transpose test data
  test<-as.data.frame(t(te_exp))
  head(test[1:10],2)
  
  ref1<-as.matrix(ref)
  test1<-as.matrix(test)
  
  head(ref1[1:10],2)
  norm_ref<- round(normalize.quantiles(ref1),3)
  #colnames(norm_ref) <- colnames(ref1)
  #rownames(norm_ref) <- rownames(ref1)
  
  #colnames(norm_ref) <- rownames(tr_exp)
  #rownames(norm_ref) <- colnames(tr_exp)
  target <- normalize.quantiles.determine.target(norm_ref)
  #target
  tt <- round(normalize.quantiles.use.target(test1,target),3)
  
  head(norm_ref)
  dim(norm_ref)
  #tt
  #rownames(tt) <- colnames(te_exp)
  #colnames(tt) <- rownames(te_exp)
  
  
  ##transpose quantile train data
  norm_ref_t <- as.data.frame(t(norm_ref))
  colnames(norm_ref_t) <- colnames(tr_exp)
  rownames(norm_ref_t) <- rownames(tr_exp)
  head(norm_ref_t[1:10],2)
  #transpose quantile test data
  test_t <- as.data.frame(t(tt))
  colnames(test_t) <- colnames(te_exp)
  rownames(test_t) <- rownames(te_exp)
  head(test_t[1:10],2)
  #tt
  
  #combine clin and normalzed data
  Train_norm_data <- cbind(tr_clin , norm_ref_t)
  Test_norm_data <- cbind(te_clin , test_t)
  
  #write normalized data into files
  #write.table(cbind("ID"=rownames(Train_norm_data), Train_norm_data),file="Log_scaled_and_QN_train_exp_data.txt",sep="\t",quote=F, row.names=FALSE)
  
  write.table(cbind("ID"=rownames(Train_norm_data), Train_norm_data),file=train_Normalized_data_clin_data,sep="\t",quote=F, row.names=FALSE)
  
  #write.table(cbind("ID"=rownames(Test_norm_data), Test_norm_data),file="Log_scaled_and_QN_test_exp_data.txt",sep="\t",quote=F, row.names=FALSE)
  write.table(cbind("ID"=rownames(Test_norm_data), Test_norm_data),file=test_Normalized_data_clin_data,sep="\t",quote=F, row.names=FALSE)
  
  
  
  
}
