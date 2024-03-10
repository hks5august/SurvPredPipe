
#' This function converting OS time (in days) into months and removing samples where OS time information is missing. Note: In the Example data OS time labelled is "OS.time" and OS event labelled as "OS"
#' @param data:args1 - data (Patients data with clinical and gene expression, where samples are in rows and features/genes are in columns)
#' @param col_num:args2 - column number in data at where clinical info ends 
#' @param output : arg3 - name of output file, in chich user want to store data
 #' @examples
#' SurvPredPipe::data_process_f(data="TCGA-LGG_protein_coding_FPKM_data_with_clin_data.txt",col_num=20, output="New_data.txt") 
#' Usage: data_process_f(data,  col_num, output)
#' @export
#' 

data_process_f <- function(data,  col_num, output) {
  
    
    # Check if any input variable is empty
    if (length(data) == 0 || length(col_num) == 0 || length(output) == 0) {
      stop("Error: Empty input variable detected.")
    }
    
    # Check if any input variable is missing
    if (any(is.na(data)) || any(is.na(col_num))) {
      stop("Error: Missing values in input variables.")
    }
    
  
  
  data <- read.table(data, header = T, sep="\t", row.names=1, check.names = F)
  
  n<- col_num - 1
  
  #Extract Clinical data
  #tr_clin <- training[1:19]
  
  
  data_clin <- data[1:n]
  

  #Convert days into months (OS time)
  data_clin1 <- data_clin %>% mutate(OS_month = round(OS.time/30.417, digit=0))
  
  
  
  ##Extract Expression data
  data_exp <- data[col_num:ncol(data)]
  
    
  # Combine clinical and Expression data
  data1 <- cbind(data_clin1 , data_exp)
  
  dim( data1)
  
  #remove samples with missing OS data or negative OS time
  data2  <- subset(data1, OS!="NA")
  data2 <- subset(data2, OS_month > 0)
  dim(data2)
  
  
  
  ###################################### Write into files ########################################
  
  #write.table(data.frame('ID'=rownames(data2), data2) , file = "Data_after_removing_samples_without_OS_time_information.txt", sep="\t", quote=F, row.names = T)
  
  write.table(data.frame('ID'=rownames(data2), data2) , file = output, sep="\t", quote=F, row.names = FALSE)
  
  
}


