
#' This function split data into training and test data as per user defined ratio
#'
#' @param data: a dataframe where features in the columns and samples must be in rows
#' @param  fraction, by which user want to split data for training data; e.g. 90% training, so fraction would be 0.9
#' @param train_data, name of training set that user can provide
#'  @param test_data, name of test set that user can provide
#' @return  training and test data 
#'
#' @examples
#'SurvPredPipe::tr_test_f(data="New_data.txt",fraction=0.9, train_data="train_FPKM.txt", test_data="test_FPKM.txt") 
#'Usage: tr_test_f(data, fraction, train_data, test_data)
#' @export


tr_test_f <- function(data, fraction, train_data, test_data)
  # tr_test_f <- function(data, 0.8)
  #  tr_test_f <- function(args[1], args[2])
{
  
  
  # Check if any input variable is empty
  if (length(data) == 0 || length(fraction) == 0 || length(train_data) == 0 || length(test_data) == 0) {
    stop("Error: Empty input variable detected.")
  }
  
  # Check if any input variable is missing
  if (any(is.na(data)) || any(is.na(fraction))) {
    stop("Error: Missing values in input variables.")
  }
  
  #set.seed(7)
  
  #data <- read.table("GBM.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  data <- read.table(data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  numberTrain <- floor(nrow(data)* fraction)
  
  trInd <- sample(1:nrow(data), numberTrain)
  training <- data[trInd,]
  testing <- data[-trInd,]
  
  
  ######################### Write into files  #######################################
  
  #write.table(training , file = "train_data.txt", sep="\t", quote=F, row.names = T)
  write.table(data.frame('ID'=rownames(training), training) , file = train_data, sep="\t", quote=F, row.names = FALSE)
  #write.table(testing , file = "test_data.txt", sep="\t", quote=F, row.names = T)
  write.table(data.frame('ID'=rownames(testing), testing) , file = test_data, sep="\t", quote=F, row.names = FALSE)
  
  
}


