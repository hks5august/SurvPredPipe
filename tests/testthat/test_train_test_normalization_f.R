# Load necessary libraries
library(testthat)
library(SurvPredPipe)  # Load your package

# Load your function
source("/Users/kaurh8/Documents/Survival_Pred_Package/Surv_Pred_Functions/SurvPredPipe/R/train_test_normalization_f.R")

# Define test cases
test_that("train_test_normalization_f function works as expected", {
  
  # Test case 1: Ensure function works with valid input
  train_data <- "train_FPKM.txt"
  test_data <- "test_FPKM.txt"
  col_num <- 21
  train_clin_data <- "Train_Clin.txt"
  test_clin_data <- "TestClin.txt"
  train_Normalized_data_clin_data <- "Train_Norm_data.txt"
  test_Normalized_data_clin_data <- "Test_Norm_data.txt"
  
  # Run the function
  train_test_normalization_f(train_data, test_data, col_num, train_clin_data, test_clin_data, train_Normalized_data_clin_data, test_Normalized_data_clin_data)
  
  # Assert that the output files were created
  expect_true(file.exists(train_clin_data))
  expect_true(file.exists(test_clin_data))
  expect_true(file.exists(train_Normalized_data_clin_data))
  expect_true(file.exists(test_Normalized_data_clin_data))
  
  # Cleanup: remove the output files
  file.remove(train_clin_data, test_clin_data, train_Normalized_data_clin_data, test_Normalized_data_clin_data)
  
  # Test case 2: Ensure function raises an error with empty input
  expect_error(train_test_normalization_f(train_data = NULL, test_data = NULL, col_num = NULL, train_clin_data = NULL, test_clin_data = NULL, train_Normalized_data_clin_data = NULL, test_Normalized_data_clin_data = NULL),
               "Error: Empty input variable detected.")
  
  # Test case 3: Ensure function raises an error with missing values in input
  expect_error(train_test_normalization_f(train_data = train_data, test_data = test_data, col_num = col_num, train_clin_data = NULL, test_clin_data = NULL, train_Normalized_data_clin_data = NULL, test_Normalized_data_clin_data = NULL),
               "Error: Missing values in input variables.")
})

# Run the tests
test_file("test_train_test_normalization_f.R")

