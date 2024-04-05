# Load necessary libraries
library(testthat)
library(SurvPredPipe)  # Load your package
# Load your function
source("/Users/kaurh8/Documents/Survival_Pred_Package/Surv_Pred_Functions/SurvPredPipe/R/tr_test_f.R")

# Define test cases
test_that("tr_test_f function works as expected", {
  
  # Test case 1: Ensure function works with valid input
  data <- "./extdata/New_data.txt"
  fraction <- 0.9
  train_data <- "train_data.txt"
  test_data <- "test_data.txt"
  
  # Run the function
  tr_test_f(data, fraction, train_data, test_data)
  
  # Assert that the output files were created
  expect_true(file.exists(train_data))
  expect_true(file.exists(test_data))
  
  # Cleanup: remove the output files
  file.remove(train_data)
  file.remove(test_data)
  
  # Test case 2: Ensure function raises an error with empty input
  expect_error(tr_test_f(data = NULL, fraction = NULL, train_data = NULL, test_data = NULL),
               "Error: Empty input variable detected.")
  
  # Test case 3: Ensure function raises an error with missing values in input
  expect_error(tr_test_f(data = data, fraction = fraction, train_data = NULL, test_data = NULL),
               "Error: Missing values in input variables.")
})

# Run the tests
test_file("test_tr_test_f.R")

