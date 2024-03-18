
#' This function converting OS time (in days) into months and removing samples where OS time information is missing. Note: In the Example data OS time labelled is "OS.time" and OS event labelled as "OS"
#' @param data:args1 - data (Patients data with clinical and gene expression, where samples are in rows and features/genes are in columns)
#' @param col_num:args2 - column number in data at where clinical info ends 
#' @param surv_time:arg3 - name of column which contain survival time (in days) information
#' @param output:arg4 - name of output file, in chich user want to store data
 #' @examples
#' SurvPredPipe::data_process_f(data="TCGA-LGG_protein_coding_FPKM_data_with_clin_data.txt",col_num=20, surv_time= "OS.time",output="New_data.txt") 
#' Usage: data_process_f(data,  col_num,surv_time, output)
#' @export
#' 

data_process_f <- function(data,  col_num, surv_time, output) {
  
    
    # Check if any input variable is empty
    if (length(data) == 0 || length(col_num) == 0 || length(surv_time) == 0 || length(output) == 0) {
      stop("Error: Empty input variable detected.")
    }
    
    # Check if any input variable is missing
    if (any(is.na(data)) || any(is.na(col_num)) || any(is.na(surv_time)) || any(is.na(output)) ) {
      stop("Error: Missing values in input variables.")
    }
    
  
  
  data <- read.table(data, header = T, sep="\t", row.names=1, check.names = F)
 
  
  n<- col_num - 1
  
  #Extract Clinical data
  #tr_clin <- training[1:19]
  
  
  data_clin <- data[1:n]
  
 
  # Column name
  #OS.time1 <- "OS.time"
  OS.time1 <- surv_time
  
  # Access column using column name
  column_data <- data_clin[[OS.time1]]
  
  #Convert days into months (OS time)
  data_clin$OS_month <- round(column_data /30.417, 0)

  
  
  ##Extract Expression data
  data_exp <- data[col_num:ncol(data)]
  
    
  # Combine clinical and Expression data
  data1 <- cbind(data_clin , data_exp)
  
  dim(data1)
  
  #remove samples with missing OS data or negative OS time
 #data2  <- subset(data1, OS!="NA")
  data2 <- subset(data1, OS_month > 0)
  dim(data2)
  
  
  
  ###################################### Write into files ########################################
  
  #write.table(data.frame('ID'=rownames(data2), data2) , file = "Data_after_removing_samples_without_OS_time_information.txt", sep="\t", quote=F, row.names = T)
  
  write.table(data.frame('ID'=rownames(data2), data2) , file = output, sep="\t", quote=F, row.names = FALSE)
  
  
}


