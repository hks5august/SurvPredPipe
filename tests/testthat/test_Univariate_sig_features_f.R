library(testthat)

library(SurvPredPipe)  # Load your package

# Load the function
source("/Users/kaurh8/Documents/Survival_Pred_Package/Surv_Pred_Functions/SurvPredPipe/R/Univariate_sig_features_f.R")
# Define the test cases
test_that("Univariate_sig_features_f function test", {
  # Test case 1: Check for empty input variables
  expect_error(Univariate_sig_features_f(
    train_data = NULL,
    test_data = "Test_Norm_data.txt",
    col_num = 21,
    surv_time = "OS_month",
    surv_event = "OS",
    output_univariate_train = "Train_Uni_sig_data.txt",
    output_univariate_test = "Test_Uni_sig_data.txt"
  ), "Error: Empty input variable detected.")

  # Test case 2: Check for missing values in input variables
  expect_error(Univariate_sig_features_f(
    train_data = "Train_Norm_data.txt",
    test_data = NA,
    col_num = 21,
    surv_time = "OS_month",
    surv_event = "OS",
    output_univariate_train = "Train_Uni_sig_data.txt",
    output_univariate_test = "Test_Uni_sig_data.txt"
  ), "Error: Missing values in input variables.")

  # Add more test cases as needed
})

test_file("test_Univariate_sig_features_f.R")
