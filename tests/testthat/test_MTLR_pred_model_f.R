library(testthat)
library(SurvPredPipe)  # Load your package

# Load the function
source("/Users/kaurh8/Documents/Survival_Pred_Package/Surv_Pred_Functions/SurvPredPipe/R/MTLR_pred_model_f.R")
# Define the test cases
test_that("MTLR_pred_model_f function test", {
  # Test case 1: Check for empty input variables
  expect_error(MTLR_pred_model_f(
    train_clin_data = NULL,
    test_clin_data = "TestClin.txt",
    Model_type = 2,
    train_features_data = "Train_PI_data.txt",
    test_features_data = "Test_PI_data.txt",
    Clin_Feature_List = "Key_PI_list.txt",
    surv_time = "OS_month",
    surv_event = "OS"
  ), "Error: Empty input variable detected.")

  # Test case 2: Check for missing values in input variables
  expect_error(MTLR_pred_model_f(
    train_clin_data = "Train_Clin.txt",
    test_clin_data = NA,
    Model_type = 2,
    train_features_data = "Train_PI_data.txt",
    test_features_data = "Test_PI_data.txt",
    Clin_Feature_List = "Key_PI_list.txt",
    surv_time = "OS_month",
    surv_event = "OS"
  ), "Error: Missing values in input variables.")

  # Add more test cases as needed
})
test_file("test_MTLR_pred_model_f.R")
