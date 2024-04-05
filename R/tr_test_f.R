
#' This function split data into training and test data as per user defined ratio
#'
#' @param data: a dataframe where features in the columns and samples must be in rows
#' @param  fraction, by which user want to split data for training data; e.g. 90% training, so fraction would be 0.9
#' @param train_data, name of training set that user can provide
#' @param test_data, name of test set that user can provide
#' @return  training and test data 
#' @import MASS
#' @import dplyr
#' @examples
#' tr_test_f(data="../inst/extdata/New_data.txt",fraction=0.9, train_data="train_FPKM.txt",test_data="test_FPKM.txt")
#'Usage: tr_test_f(data, fraction, train_data, test_data)
#' @export


tr_test_f <- function(data, fraction, train_data, test_data)
{
  
  
  # Check if any input variable is empty
  if (length(data) == 0 || length(fraction) == 0 || length(train_data) == 0 || length(test_data) == 0) {
    stop("Error: Empty input variable detected.")
  }
  
  # Check if any input variable is missing
  if (any(is.na(data)) || any(is.na(fraction))) {
    stop("Error: Missing values in input variables.")
  }
  
#load data  
data <- read.table(data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

 #define train-test split fraction
 numberTrain <- floor(nrow(data)* fraction)
  
 #extract training index
  trInd <- sample(1:nrow(data), numberTrain)
  
  #create training data
  training <- data[trInd,]
  #create test data
  testing <- data[-trInd,]
  
  
  ######################### Write into files  #######################################
  
  #write.table(training , file = "train_data.txt", sep="\t", quote=F, row.names = T)
  write.table(data.frame('ID'=rownames(training), training) , file = train_data, sep="\t", quote=F, row.names = FALSE)
  #write.table(testing , file = "test_data.txt", sep="\t", quote=F, row.names = T)
  write.table(data.frame('ID'=rownames(testing), testing) , file = test_data, sep="\t", quote=F, row.names = FALSE)
  
  
}


