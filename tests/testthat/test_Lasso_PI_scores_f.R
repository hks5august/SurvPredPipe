# Load necessary libraries
library(testthat)
library(SurvPredPipe)  # Load your package
# Load your function
source("/Users/kaurh8/Documents/Survival_Pred_Package/Surv_Pred_Functions/SurvPredPipe/R/Lasso_PI_scores_f.R")

# Define test cases
test_that("Lasso_PI_scores_f function works as expected", {
  
  # Test case 1: Ensure function works with valid input
  train_data <- "Train_Norm_data.txt"
  test_data <- "Test_Norm_data.txt"
  nfolds <- 5
  col_num <- 21
  surv_time <- "OS_month"
  surv_event <- "OS"
  train_PI_data <- "Train_PI_data.txt"
  test_PI_data <- "Test_PI_data.txt"
  
  # Run the function
  Lasso_PI_scores_f(train_data, test_data, nfolds, col_num, surv_time, surv_event, train_PI_data, test_PI_data)
  
  # Assert that the output files were created
  expect_true(file.exists(train_PI_data))
  expect_true(file.exists(test_PI_data))
  
  # Cleanup: remove the output files
  file.remove(train_PI_data)
  file.remove(test_PI_data)
  
  # Test case 2: Ensure function raises an error with empty input
  expect_error(Lasso_PI_scores_f(train_data = NULL, test_data = NULL, nfolds = NULL, col_num = NULL, surv_time = NULL, surv_event = NULL, train_PI_data = NULL, test_PI_data = NULL),
               "Error: Empty input variable detected.")
  
  # Test case 3: Ensure function raises an error with missing values in input
  expect_error(Lasso_PI_scores_f(train_data = train_data, test_data = test_data, nfolds = nfolds, col_num = col_num, surv_time = surv_time, surv_event = surv_event, train_PI_data = NULL, test_PI_data = NULL),
               "Error: Missing values in input variables.")
})

# Run the tests
test_file("test_Lasso_PI_scores_f.R")

