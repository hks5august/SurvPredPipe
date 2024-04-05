
#' This function selects significant features (p-value <0.05)s based on median expression values of features using Univariate Survival analysis.
#' @param train_data :args1 - training data (Patients data with clinical and gene expression, where samples are in rows and features/genes are in columns)
#' @param test_data_data :args2 - training data (Patients data with clinical and gene expression, where samples are in rows and features/genes are in columns)
#' @param col_num :args3 - column number in data at where clinical info ends
#' @param surv_time :arg4 - name of column which contain survival time (in days) information
#' @param surv_event :arg5 - name of column which contain survival eventinformation
#' @param output_univariate_train :args6- name of output to store univariate selected features for training data
#' @param output_univariate_test :args7- name of output to store univariate selected features for test data
#' @import dplyr
#' @import survival
#' @import survminer
#' @import ggplot2
#' @examples
#' Univariate_sig_features_f(train_data="../inst/extdata/Train_Norm_data.txt", test_data="../inst/extdata/Test_Norm_data.txt", col_num=21, surv_time="OS_month" , surv_event="OS" ,output_univariate_train="Train_Uni_sig_data.txt", output_univariate_test="Test_Uni_sig_data.txt") 
#' Usage: Univariate_sig_features_f(train_data, test_data, col_num, surv_time, surv_event, output_univariate_train, output_univariate_test)
#' @export



#setwd("/Users/kaurh8/Documents/GDC_TCGA_Biolinks/GDC_All_samples_Data/updated_folder2/LGG")
Univariate_sig_features_f <- function(train_data, test_data, col_num,surv_time, surv_event,  output_univariate_train, output_univariate_test)  {

  
  # Check if any input variable is empty
  if (length(train_data) == 0 ||  length(test_data) == 0|| length(col_num) == 0 ||length(surv_time) == 0 ||length(surv_event) == 0 ||length(output_univariate_train) == 0 ||length(output_univariate_test) == 0  ) {
    stop("Error: Empty input variable detected.")
  }
  
  # Check if any input variable is missing
  if (any(is.na(train_data)) || any(is.na(test_data))  || any(is.na(col_num))  || any(is.na(surv_time)) || any(is.na(surv_event)) || any(is.na(output_univariate_train)) || any(is.na(output_univariate_test)) ) {
    stop("Error: Missing values in input variables.")
  }



#load data
tr_data1 <- read.table(train_data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
te_data1 <- read.table(test_data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)


# rename survival time and event column name
colnames(tr_data1)[colnames(tr_data1) == surv_time] <- "OS_month"
colnames(tr_data1)[colnames(tr_data1) == surv_event] <- "OS"



# rename survival time and event column name
colnames(te_data1)[colnames(te_data1) == surv_time] <- "OS_month"
colnames(te_data1)[colnames(te_data1) == surv_event] <- "OS"


###################### Create Survival Object #####################
surv_object <- Surv(time = tr_data1$OS_month, event = tr_data1$OS)
# create a file to store results


n= col_num -1
### Significant results

write.table(cbind("ID","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="./Univariate_Survival_Significant_genes_List.txt",row.names=F,col.names=F,sep = '\t');


#perform univariate survival analysis for each feature based on median expression
for(i in seq(from=col_num, to=length(tr_data1), by=1))

{
#create survival object
  surv_object <- Surv(time = tr_data1$OS_month, event = tr_data1$OS)
  
  #survival analysis: fits cox ph model to find HR for median cut
  fit1 <- survfit(surv_object ~ (tr_data1[,i])>(median(tr_data1[1,i])), data=tr_data1);
  summary(fit1);
  fit1.coxph <- coxph(surv_object ~ (tr_data1[,i])>(median(tr_data1[1,i])), data = tr_data1)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  
  #check whether the pvalue is significant (< 0.05) or not
  if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))

#append results into file
  {write.table(cbind(colnames(tr_data1[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="./Univariate_Survival_Significant_genes_List.txt",row.names=F,col.names=F,sep = '\t',append = T);#output file
  }
}


#load selected Univariate features
sel_univariate_features<-read.table("Univariate_Survival_Significant_genes_List.txt",header =TRUE, sep = "\t", row.names=1, check.names = FALSE)


#training data preparation with selected features 
sel_univ_train <- as.data.frame(tr_data1[,colnames(tr_data1) %in% c(row.names(sel_univariate_features)), ])

#test data preparation with selected features
sel_univ_test<- as.data.frame(te_data1[,colnames(te_data1) %in% c(row.names(sel_univariate_features)), ])


#write output into files
#write.table(sel_univ_train , file = "sel_univariate_train_data.txt", sep="\t", quote=F, row.names = T)

write.table(data.frame('ID'=rownames(sel_univ_train), sel_univ_train) , file = output_univariate_train, sep="\t", quote=F, row.names = FALSE)
#write.table(sel_univ_test , file = "sel_univariate_test_data.txt", sep="\t", quote=F, row.names = T)

write.table(data.frame('ID'=rownames(sel_univ_test), sel_univ_test) , file = output_univariate_test, sep="\t", quote=F, row.names = FALSE)



}

