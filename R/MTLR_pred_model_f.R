
#' This function can generate 5 Prediction models types based on :(1)Clinical features,  (2) PI score , (3) PI score + Clin, (4) Significant Univariate features, (5) Significant Univariate features + clin features, using MTLR method.
#' Further, this model will predict the survival of test data in terms of survival probability over different time points, mean and median survival time pf patients in test data.

#' #'@param train_clin_data :args1 - training data with Clin features (Patients data with clinical and gene expression, where samples are in rows and features/genes are in columns)
#' @param test_clin_data :args2 -  test data with Clin features  (Patients data with clinical and gene expression, where samples are in rows and features/genes are in columns)
#' @param  Model_type: args3 - Prediction Model type : 1 - Clinical features, 2 - PI score features, 3-  PI score + Clinical features ; 4 - Univariate features, 5 -  Univariate features + Clinical features
#' @param train_features_data :args4 - Training data with selected features only (Significant Univariate/LASSO PI score)
#' @param test_features_data :args5 - Training data with selected features only (Significant Univariate/LASSO PI score)
#' @param Clin_Feature_List :args6 - List of key clinical features which contain non-unique values (e.g. Key_Clin_feature_list.txt, Key_PI_list.txt Key_Clin_features_with_PI_list.txt, Key_univariate_features_list.txt, Key_univariate_features_with_Clin_list.txt)
#' 
#' @examples
#' SurvPredPipe::MTLR_pred_model_f(train_clin_data = "Train_Clin.txt", test_clin_data = "TestClin.txt", Model_type = 2, train_features_data = "Train_PI_data.txt", test_features_data = "Test_PI_data.txt" ) 
#' Usage: MTLR_pred_model_f(train_clin_data, test_clin_data, Model_type, train_features_data, test_features_data )
#' @export


MTLR_pred_model_f <- function(train_clin_data, test_clin_data, Model_type, train_features_data, test_features_data, Clin_Feature_List )  {
#set.seed(7)

  # Check if any input variable is empty
  if (length(train_clin_data) == 0 ||  length(test_clin_data) == 0 || length(Model_type) == 0   || length(train_features_data) == 0 ||length(test_features_data ) == 0 ||length(Clin_Feature_List ) == 0) {
    stop("Error: Empty input variable detected.")
  }
  
  # Check if any input variable is missing
  if (any(is.na(train_clin_data)) || any(is.na(test_clin_data)) || any(is.na(Model_type))  || any(is.na(train_features_data)) || any(is.na(test_features_data )) || any(is.na(Clin_Feature_List))) {
    stop("Error: Missing values in input variables.")
  }
  
  
  #load data
 # tr_clin1 <- read.table("Train_Clin.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  #te_clin1 <- read.table("Test_Clin.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  tr_clin1 <- read.table(train_clin_data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  te_clin1 <- read.table(test_clin_data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  
  #train_features_data1 <- read.table("Train_PI_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  #test_features_data1 <- read.table("Test_PI_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  
  #train_features_data1 <- read.table("Train_Uni_sig_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  #test_features_data1 <- read.table("Test_Uni_sig_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  
  train_features_data1 <- read.table(train_features_data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  test_features_data1 <- read.table(test_features_data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  
  #combine clinical and feature data
  
  tr_data2 <- cbind(tr_clin1, train_features_data1)
  te_data2 <- cbind(te_clin1, test_features_data1 )  
  
  # Load user defined a list of features forclinical data
  #ftr_list <- read.table("Key_Clin_feature_list.txt", header = TRUE, sep = "\t",check.names = FALSE)
  #ftr_list <- read.table("Key_PI_list.txt", header = TRUE, sep = "\t",check.names = FALSE)
  #ftr_list <- read.table("Key_Clin_features_with_PI_list.txt", header = TRUE, sep = "\t",check.names = FALSE)

  #ftr_list <- read.table("Key_univariate_features_list.txt", header = TRUE, sep = "\t",check.names = FALSE)
  #ftr_list <- read.table("Key_univariate_features_with_Clin_list.txt", header = TRUE, sep = "\t",check.names = FALSE)
  ftr_list <- read.table(Clin_Feature_List, header = TRUE, sep = "\t",check.names = FALSE)
  
    ##############
  
 # Model1 == 1 #Clin_model
 # Model2 == 2 #PI_model
 
 # Model3 == 3 #PI_with_Clin_model
 # Model4 == 4 #Univ_model
 # Model5 == 5 #PI_with_Clin_model
  
if(Model_type == 1){
  

############################# Selected Clin features models #######################
############################# MTLR Clinical features ##############################

# Load user defined a list of features forclinical data
#ftr_list <- read.table("Key_Clin_feature_list.txt", header = TRUE, sep = "\t",check.names = FALSE)
ftr_list <- read.table(Clin_Feature_List, header = TRUE, sep = "\t",check.names = FALSE)

#create data frame with selected features (user provided list)
sel_clin_tr <- as.data.frame(tr_clin1[,colnames(tr_clin1) %in% c(ftr_list$ID), ])

sel_clin_te<- as.data.frame(te_clin1[,colnames(te_clin1) %in% c(ftr_list$ID), ])

#add survival information
sel_clin_tr1  <- cbind(tr_clin1["OS"], tr_clin1["OS_month"],  sel_clin_tr)
sel_clin_te1  <- cbind(te_clin1["OS"], te_clin1["OS_month"],  sel_clin_te)

#sel_clin_tr1$surv_obj <- Surv(sel_clin_tr1$OS_month,sel_clin_tr1$OS)
formula1 <- Surv(OS_month,OS) ~ .
#Next, we just need the data argument which in our case is training. We can finally make our first model!
formula1
#Next, we just need the data argument which in our case is training. We can finally make our first model!   
Mod1<- mtlr(formula = formula1, data =  sel_clin_tr1)
#Mod1 <- mtlr(surv_obj ~ ., data = sel_clin_tr1)
Mod1 

#save model
save(Mod1, file = "Model_with_Clin_ftrs.RData")
#load saved model
#load("Model_with_Clin_ftrs.RData")

#mtlr outputs the weight matrix for the model – 
#these are the weights corresponding to each feature at each time point (additionally notice that we include the bias weights). 
#The row names correspond to the time point for which the feature weight belongs. 
#If you would like to access these weights, they are saved in the model object as weight_matrix so you can access them using smallMod$weight_matrix.

#We can also plot the weights for a mtlr model. Before we printed the small model but here we will look at the weights for the complete model.



sel_clin_te2 <- na.omit(sel_clin_te1)
############  Prediction on Test Data #############

dim(sel_clin_te2)

survCurves1 <- predict(Mod1, sel_clin_te1, type = "survivalcurve")
#survCurves is pretty large so we will look at the first 5 rows/columns.

dim(survCurves1)
head(survCurves1)
dim(sel_clin_te2)

colnames(survCurves1)<- c("time_point", rownames(sel_clin_te2))


survCurves1_df <- as.data.frame(survCurves1)
head(survCurves1_df )

#write.table(survCurves1_df , file = "survCurves1_data.txt", sep="\t", quote=F, row.names = F)
#When we use the predict function for survival curves we will be returned a matrix where the first column (time) is the list of time points that the model evaluated the survival probability for each observation (these will be the time points used by mtlr and an additional 0 point). Every following column will correspond to the row number of the data passed in, e.g. column 2 (named 1) corresponds to row 1 of testing. Each row of this matrix gives the probabilities of survival at the corresponding time point (given by the time column). For example, testing observation 1 has a survival probability of 0.919 at time 60.625.

#Since these curves may be hard to digest by observing a matrix of survival probabilities we can also choose to plot them.


#jpeg(file="survival_curves_test_data_based_on_selected_clin.jpeg", units="in", width=10, height=10, res=300)
#plotcurves(survCurves, 1:63)
#dev.off()

#Mean/Median Survival Time: In addition to the entire survival curve one may also be interested in the average survival time. 
#This is again available from the predict function.

#survivalcurve
survivalcurve1 <- predict(Mod1, sel_clin_te2, type = "survivalcurve")
head(survivalcurve1)



#Prob_Event
Survival_prob_event1 <- predict(Mod1, sel_clin_te2, type = "prob_event")
head(Survival_prob_event1)


#Mean
meanSurv1 <- predict(Mod1, sel_clin_te2, type = "mean_time")
head(meanSurv1)

#> [1] 319.0312 407.2113 262.0245 317.3220 327.1613 336.3179
#Median
medianSurv1 <- predict(Mod1, sel_clin_te2, type = "median_time")
head(medianSurv1)


meanSurv_d1  <- as.data.frame(meanSurv1)
medianSurv_d1  <- as.data.frame(medianSurv1)

meanSurv_d1
medianSurv_d1

names1 <- as.data.frame(rownames(sel_clin_te2))
names1
#colnames(names1) <- c("ID")
mean_median_surv1_d <- cbind(names1, meanSurv_d1 , medianSurv_d1 )
rownames(mean_median_surv1_d)<- c(rownames(sel_clin_te2))

colnames(mean_median_surv1_d) <- c("IDs", "Mean", "Median")
head(mean_median_surv1_d )


#write.table(mean_median_surv1_d , file = "mean_median_time_data.txt", sep="\t", quote=F, row.names = F)



#Survival Probability at Event Time: The last prediction type supported is acquiring the observations survival probability at the respective event time. However, in order to use this prediction, the event time (whether censored or uncensored) must be included in the features passed into the predict function.

#Survival_Probs_event1<- predict(Mod1, sel_clin_te, type = "prob_event")
#head(Survival_Probs_event)

survivalProbs_p1 <- predict(Mod1, sel_clin_te1, type = "prob_times")
head(survivalProbs_p1,2)



#extract prob times at diff times points
survivalProbs_t_mat <-as.matrix(survivalProbs_p1 )
head(survivalProbs_t_mat,2)
survivalProbs_t_mat1 <- survivalProbs_t_mat[,-1]
head(survivalProbs_t_mat1) 
survivalProbs_t_mat1_t<- t(survivalProbs_t_mat1 )
head(survivalProbs_t_mat1_t,2)
survivalProbs_t_mat1_t2 <- survivalProbs_t_mat1_t[,-1]
head(survivalProbs_t_mat1_t2)

meanSurv1
medianSurv1

surv_res1 <- cbind(meanSurv1, medianSurv1, Survival_prob_event1 ,  sel_clin_te2$OS_month, sel_clin_te2$OS )
dim(surv_res1)
colnames(surv_res1) <- c("Mean_Surv", "Median_surv",  "Prob_event", "Actual_OS_time", "OS_event")
rownames(surv_res1) <- rownames(sel_clin_te2)
head(surv_res1,10)

#write.table( surv_res1 , file = "surv_resMTLR_selected_clin.txt", sep="\t", quote=F, row.names = T)


surv_obj1 <- Surv(sel_clin_te2$OS_month, sel_clin_te2$OS)
surv_obj1


IBS_1 <- round(IBS(surv_obj1, sp_matrix = survivalProbs_t_mat1_t2, survivalProbs_p1$time[-1]),3)
IBS_1

#Concordance Index

#SurvMetrics::Cindex(surv_obj1, predicted = medianSurv)
c_index1<- round(SurvMetrics::Cindex(surv_obj1, predicted = medianSurv1),2)
c_index1


Error_mat_1 <- cbind(IBS_1,  c_index1)
colnames(Error_mat_1) <- c("IBS_score", "c_index")
Error_mat_1

#write.table(Error_mat , file = "Error_mat_sel_clin.txt", sep="\t", quote=F, row.names = T)
}
  
  
########################## Model with only PI score ###############################

  else if (Model_type == 2) {

    
    #create data frame with selected features (user provided list)
    #sel_clin_tr <- as.data.frame(tr_data2[,colnames(tr_data2) %in% c(ftr_list$ID), ])
    
    #sel_clin_te<- as.data.frame(te_data2[,colnames(te_data2) %in% c(ftr_list$ID), ])
    
    #add survival information
    sel_clin_tr1  <- cbind(tr_clin1["OS"], tr_clin1["OS_month"],  tr_data2["PI"])
    sel_clin_te1  <- cbind(te_clin1["OS"], te_clin1["OS_month"],  te_data2["PI"])
    

  
#formula2<- Surv(OS_month,OS) ~ PI
    formula2 <- Surv(OS_month,OS) ~ .
    #Next, we just need the data argument which in our case is training. We can finally make our first model!
    formula2
    #Next, we just need the data argument which in our case is training. We can finally make our first model!   
    Mod2<- mtlr(formula = formula2, data =  sel_clin_tr1)
    
    

  #  Mod2 <- mtlr(surv_obj ~ ., data = sel_clin_tr1)

#smallMod <- mtlr(formula = formulaSmall, data = tr_clin3)
#We will print the small model so the output is more compact.
Mod2
save(Mod2, file = "Model_with_PI.RData")
#load saved model
#load("Model_with_PI.RData")

#mtlr outputs the weight matrix for the model – 
#these are the weights corresponding to each feature at each time point (additionally notice that we include the bias weights). 
#The row names correspond to the time point for which the feature weight belongs. 
#If you would like to access these weights, they are saved in the model object as weight_matrix so you can access them using smallMod$weight_matrix.

#We can also plot the weights for a mtlr model. Before we printed the small model but here we will look at the weights for the complete model.

plot(Mod2)

#Model Predictions: Now that we have trained a MTLR model we should make some predictions! 
#This is where our testing set and the predict function will come into play.
#Note that there are a number of predictions we may be interested in acquiring. 
#First, we may want to view the survival curves of our test observations.

sel_clin_te2 <- na.omit(sel_clin_te1)

survCurves2 <- predict(Mod2, sel_clin_te2  , type = "survivalcurve")
#survCurves is pretty large so we will look at the first 5 rows/columns.

colnames(survCurves2)<- c("time_point", rownames(sel_clin_te2))


survCurves2_df <- as.data.frame(survCurves2)

#write.table(survCurves2_df , file = "survCurves2_data.txt", sep="\t", quote=F, row.names = F)
#When we use the predict function for survival curves we will be returned a matrix where the first column (time) is the list of time points that the model evaluated the survival probability for each observation (these will be the time points used by mtlr and an additional 0 point). Every following column will correspond to the row number of the data passed in, e.g. column 2 (named 1) corresponds to row 1 of testing. Each row of this matrix gives the probabilities of survival at the corresponding time point (given by the time column). For example, testing observation 1 has a survival probability of 0.919 at time 60.625.

#Since these curves may be hard to digest by observing a matrix of survival probabilities we can also choose to plot them.

#plotcurves(survCurves2,1:10)


#jpeg(file="survival_curves_test_data_based_on_PI_score.jpeg", units="in", width=10, height=10, res=300)
#plotcurves(survCurves2, 1:10)
#dev.off()

#Mean/Median Survival Time: In addition to the entire survival curve one may also be interested in the average survival time. 
#This is again available from the predict function.


#survivalcurve
survivalcurve2 <- predict(Mod2, sel_clin_te2, type = "survivalcurve")
head(survivalcurve2)
dim(survivalcurve2)


#Prob_Event
Survival_prob_event2 <- predict(Mod2, sel_clin_te2 , type = "prob_event")
head(Survival_prob_event2)
dim(Survival_prob_event2)

#Mean
meanSurv2 <- predict(Mod2, sel_clin_te2, type = "mean_time")
head(meanSurv2)
dim(meanSurv2 )


#Median
medianSurv2 <- predict(Mod2, sel_clin_te2, type = "median_time")
head(medianSurv2)
dim(medianSurv2)


meanSurv2_d  <- as.data.frame(meanSurv2)
medianSurv2_d  <- as.data.frame(medianSurv2)

names_2 <- as.data.frame(rownames(sel_clin_te2))
mean_median_surv2_d <- cbind(names_2, meanSurv2_d , medianSurv2_d )
rownames(mean_median_surv2_d)<- c(rownames(sel_clin_te2))

colnames(mean_median_surv2_d) <- c("IDs", "Mean", "Median")
head(mean_median_surv2_d )


#write.table(mean_median_surv2_d , file = "mean_median_time_data2.txt", sep="\t", quote=F, row.names = F)


#Survival Probability at Event Time: The last prediction type supported is acquiring the observations survival probability at the respective event time. However, in order to use this prediction, the event time (whether censored or uncensored) must be included in the features passed into the predict function.

survivalProbs_p2 <- predict(Mod2, sel_clin_te2, type = "prob_times")
head(survivalProbs_p2,2)



#extract prob times at diff times points
survivalProbs_t_mat_2 <-as.matrix(survivalProbs_p2 )
head(survivalProbs_t_mat_2,2)
survivalProbs_t_mat1_2 <- survivalProbs_t_mat_2[,-1]
head(survivalProbs_t_mat1_2) 
survivalProbs_t_mat1_t_2<- t(survivalProbs_t_mat1_2 )
head(survivalProbs_t_mat1_t_2,2)
survivalProbs_t_mat1_t2_2 <- survivalProbs_t_mat1_t_2[,-1]
head(survivalProbs_t_mat1_t2_2)




surv_res2 <- cbind(meanSurv2, medianSurv2, Survival_prob_event2 , sel_clin_te2$OS_month, sel_clin_te2$OS )
dim(surv_res2)
colnames(surv_res2) <- c("Mean_Surv","Prob_event",  "Median_surv", "Actual_OS_time", "OS_event")
rownames(surv_res2) <- rownames(sel_clin_te2)
head(surv_res2,10)

#write.table( surv_res2 , file = "surv_resMTLR_with_PI.txt", sep="\t", quote=F, row.names = T)


surv_obj_2 <- Surv(sel_clin_te2$OS_month,sel_clin_te2$OS)
surv_obj_2


IBS1_2 <- round(IBS(surv_obj_2, sp_matrix = survivalProbs_t_mat1_t2_2, survivalProbs_p2$time[-1]),3)
IBS1_2

#Concordance Index

c_index_2<- round(SurvMetrics::Cindex(surv_obj_2, predicted = medianSurv2),2)
c_index_2

Error_mat_2 <- cbind( IBS1_2, c_index_2)
colnames(Error_mat_2) <- c("IBS_Score", "c_index")
Error_mat_2

#write.table(Error_mat_2 , file = "Error_mat2_with_PI.txt", sep="\t", quote=F, row.names = T)

####################

}
  
  else if (Model_type == 3) {  
  
#################################### PI with Clin features ################
 #create data frame with selected features (user provided list)
  sel_clin_tr <- as.data.frame(tr_data2[,colnames(tr_data2) %in% c(ftr_list$ID), ])
  sel_clin_te<- as.data.frame(te_data2[,colnames(te_data2) %in% c(ftr_list$ID), ])
    
  #add survival information
  sel_clin_tr1  <- cbind(tr_clin1["OS"], tr_clin1["OS_month"],  sel_clin_tr)
  sel_clin_te1  <- cbind(te_clin1["OS"], te_clin1["OS_month"],  sel_clin_te)
  
  head(sel_clin_tr1,2)
  
  #remove samples where informtion missing for any of selected feature
  sel_clin_te2 <- na.omit(sel_clin_te1)


formula3 <- Surv(OS_month,OS) ~ .
#Next, we just need the data argument which in our case is training. We can finally make our first model!
formula3

Mod3<- mtlr(formula = formula3, data = sel_clin_tr1 )
#smallMod <- mtlr(formula = formulaSmall, data = tr_clin3)
#We will print the small model so the output is more compact.
Mod3

save(Mod3, file = "Model_with_PI_and_Clin.RData")
#load saved model
#load("Model_with_PI_and_Clin.RData")
#mtlr outputs the weight matrix for the model – 
#these are the weights corresponding to each feature at each time point (additionally notice that we include the bias weights). 
#The row names correspond to the time point for which the feature weight belongs. 
#If you would like to access these weights, they are saved in the model object as weight_matrix so you can access them using smallMod$weight_matrix.

#We can also plot the weights for a mtlr model. Before we printed the small model but here we will look at the weights for the complete model.

plot(Mod3, numfeatures = 10)


jpeg(file="survival_curves_features_weight_PI_clin_ftes.jpeg", units="in", width=10, height=10, res=350)
plot(Mod3, numfeatures = 10)
dev.off()


#Model Predictions: Now that we have trained a MTLR model we should make some predictions! 
#This is where our testing set and the predict function will come into play.
#Note that there are a number of predictions we may be interested in acquiring. 
#First, we may want to view the survival curves of our test observations.



survCurves3 <- predict(Mod3, sel_clin_te2, type = "survivalcurve")
#survCurves is pretty large so we will look at the first 5 rows/columns.
dim(survCurves3 )
head(survCurves3,2)
survCurves3[1:5,1:5]
colnames(survCurves3)<- c("time_point", rownames(sel_clin_te2))

survCurves3_df <- as.data.frame(survCurves3)
head(survCurves3,2)
#write.table(survCurves3_df , file = "survCurves3_data.txt", sep="\t", quote=F, row.names = F)


#Mean/Median Survival Time: In addition to the entire survival curve one may also be interested in the average survival time. 
#This is again available from the predict function.



#Mean
meanSurv3 <- predict(Mod3, sel_clin_te2, type = "mean_time")
head(meanSurv3)
dim(meanSurv3)

#Median
medianSurv3 <- predict(Mod3, sel_clin_te2, type = "median_time")
head(medianSurv3)
dim(medianSurv3)

meanSurv3_d  <- as.data.frame(meanSurv3)
medianSurv3_d  <- as.data.frame(medianSurv3)

names_3 <- as.data.frame(rownames(sel_clin_te2))
mean_median_surv3_d <- cbind(names_3, meanSurv3_d , medianSurv3_d )
rownames(mean_median_surv3_d)<- c(rownames(sel_clin_te2))

colnames(mean_median_surv3_d) <- c("IDs", "Mean", "Median")
head(mean_median_surv3_d )


#write.table(mean_median_surv3_d , file = "mean_median_time_data3.txt", sep="\t", quote=F, row.names = F)

#Survival Probability at Event Time: The last prediction type supported is acquiring the observations survival probability at the respective event time. However, in order to use this prediction, the event time (whether censored or uncensored) must be included in the features passed into the predict function.

Survival_Probs_event3 <- predict(Mod3, sel_clin_te2, type = "prob_event")
head(Survival_Probs_event3)
dim(Survival_Probs_event3)


survivalProbs_p3 <- predict(Mod3, sel_clin_te2, type = "prob_times")
head(survivalProbs_p3)
dim(survivalProbs_p3)

surv_res3 <- cbind(meanSurv3, medianSurv3, Survival_Probs_event3, sel_clin_te2$OS_month, sel_clin_te2$OS)
surv_res3
colnames(surv_res3) <- c("Mean", "Median_surv", "Survival_Prob_event", "Actual_OS_time", "Event")

head(surv_res3)

#write.table( surv_res3 , file = "surv_res1_MTLR_based_on_PI_and_Clin.txt", sep="\t", quote=F, row.names = T)


#extract prob times at diff times points
survivalProbs_t_mat_3 <-as.matrix(survivalProbs_p3 )
head(survivalProbs_t_mat_3,2)
survivalProbs_t_mat1_3 <- survivalProbs_t_mat_3[,-1]
head(survivalProbs_t_mat1_3) 
survivalProbs_t_mat1_t_3<- t(survivalProbs_t_mat1_3 )
head(survivalProbs_t_mat1_t_3,2)
survivalProbs_t_mat1_t2_3 <- survivalProbs_t_mat1_t_3[,-1]
head(survivalProbs_t_mat1_t2_3)




surv_obj_3 <- Surv(sel_clin_te2$OS_month,sel_clin_te2$OS)
surv_obj_3

IBS1_3<- round(IBS(surv_obj_3, sp_matrix = survivalProbs_t_mat1_t2_3, survivalProbs_p3$time[-1]),3)
IBS1_3

#Concordance Index

c_index_3<- round(SurvMetrics::Cindex(surv_obj_3, predicted = medianSurv3),2)
c_index_3


Error_mat_3 <- cbind(IBS1_3, c_index_3)
colnames(Error_mat_3) <- c("IBS_Score", "c_index")
Error_mat_3

#write.table(Error_mat_3 , file = "Error_mat3_with_PI_and Clin.txt", sep="\t", quote=F, row.names = T)
}

  else if (Model_type == 4) {
############ Models based on Univariate Significant Genes #######

 #create data frame with selected features (user provided list)
 sel_clin_tr <- as.data.frame(tr_data2[,colnames(tr_data2) %in% c(ftr_list$ID), ])
    
 sel_clin_te<- as.data.frame(te_data2[,colnames(te_data2) %in% c(ftr_list$ID), ])
    
 #add survival information
 sel_clin_tr1  <- cbind(tr_clin1["OS"], tr_clin1["OS_month"],  sel_clin_tr)
 sel_clin_te1  <- cbind(te_clin1["OS"], te_clin1["OS_month"],  sel_clin_te)
 
 head(sel_clin_tr1,2)
 

#formula
formula4 <- Surv(OS_month,OS) ~ .
#Next, we just need the data argument which in our case is training. We can finally make our first model!
formula4

Mod4 <- mtlr(formula = formula4, data = sel_clin_tr1)
#smallMod <- mtlr(formula = formulaSmall, data = tr_clin3)
#We will print the small model so the output is more compact.
Mod4


#mtlr outputs the weight matrix for the model – 
#these are the weights corresponding to each feature at each time point (additionally notice that we include the bias weights). 
#The row names correspond to the time point for which the feature weight belongs. 
#If you would like to access these weights, they are saved in the model object as weight_matrix so you can access them using smallMod$weight_matrix.

#We can also plot the weights for a mtlr model. Before we printed the small model but here we will look at the weights for the complete model.

# remove samples where features aree missing
sel_clin_te2 <- na.omit(sel_clin_te1)

survCurves4 <- predict(Mod4, sel_clin_te2, type = "survivalcurve")
#survCurves is pretty large so we will look at the first 5 rows/columns.

colnames(survCurves4)<- c("time_point", rownames(sel_clin_te2))


survCurves4_df <- as.data.frame(survCurves4)

head(survCurves4_df,3)

#write.table(survCurves4_df , file = "survCurves4_data.txt", sep="\t", quote=F, row.names = F)
#When we use the predict function for survival curves we will be returned a matrix where the first column (time) is the list of time points that the model evaluated the survival probability for each observation (these will be the time points used by mtlr and an additional 0 point). Every following column will correspond to the row number of the data passed in, e.g. column 2 (named 1) corresponds to row 1 of testing. Each row of this matrix gives the probabilities of survival at the corresponding time point (given by the time column). For example, testing observation 1 has a survival probability of 0.919 at time 60.625.

#Since these curves may be hard to digest by observing a matrix of survival probabilities we can also choose to plot them.

#plotcurves(survCurves4,1:10)


#jpeg(file="survival_curves_test_data_based_on_Univariate_features.jpeg", units="in", width=10, height=10, res=300)
#plotcurves(survCurves4, 1:10)
#dev.off()

#Mean/Median Survival Time: In addition to the entire survival curve one may also be interested in the average survival time. 
#This is again available from the predict function.


#survivalcurve
survivalcurve4 <- predict(Mod4, sel_clin_te2, type = "survivalcurve")
head(survivalcurve4)
dim(survivalcurve4)


#Prob_Event
Survival_prob_event4 <- predict(Mod4, sel_clin_te2 , type = "prob_event")
head(Survival_prob_event4)
dim(Survival_prob_event4)

#Mean
meanSurv4 <- predict(Mod4, sel_clin_te2, type = "mean_time")
head(meanSurv4)
dim(meanSurv4 )


#Median
medianSurv4 <- predict(Mod4, sel_clin_te2, type = "median_time")
head(medianSurv4)
dim(medianSurv4)


meanSurv4_d  <- as.data.frame(meanSurv4)
medianSurv4_d  <- as.data.frame(medianSurv4)

names_4 <- as.data.frame(rownames(sel_clin_te2))
mean_median_surv4_d <- cbind(names_4, meanSurv4_d , medianSurv4_d )
rownames(mean_median_surv4_d)<- c(rownames(sel_clin_te2))

colnames(mean_median_surv4_d) <- c("IDs", "Mean", "Median")
head(mean_median_surv4_d )


#write.table(mean_median_surv4_d , file = "mean_median_time_data4.txt", sep="\t", quote=F, row.names = F)


#Survival Probability at Event Time: The last prediction type supported is acquiring the observations survival probability at the respective event time. However, in order to use this prediction, the event time (whether censored or uncensored) must be included in the features passed into the predict function.

survivalProbs_p4 <- predict(Mod4, sel_clin_te2, type = "prob_times")
head(survivalProbs_p4,2)



#extract prob times at diff times points
survivalProbs_t_mat_4 <-as.matrix(survivalProbs_p4 )
head(survivalProbs_t_mat_4,2)
survivalProbs_t_mat1_4 <- survivalProbs_t_mat_4[,-1]
head(survivalProbs_t_mat1_4) 
survivalProbs_t_mat1_t_4<- t(survivalProbs_t_mat1_4 )
head(survivalProbs_t_mat1_t_4,2)
survivalProbs_t_mat1_t2_4 <- survivalProbs_t_mat1_t_4[,-1]
head(survivalProbs_t_mat1_t2_4)




surv_res4 <- cbind(meanSurv4, medianSurv4, Survival_prob_event4, sel_clin_te2$OS_month, sel_clin_te2$OS )
dim(surv_res4)
colnames(surv_res4) <- c("Mean_Surv", "Median_surv", "Survival_Prob", "Actual_OS_time", "OS_event")
rownames(surv_res4) <- rownames(sel_clin_te2)
head(surv_res4,10)

#write.table( surv_res4 , file = "surv_resMTLR_with_Univariate_ftrs.txt", sep="\t", quote=F, row.names = T)


surv_obj_4 <- Surv(sel_clin_te2$OS_month,sel_clin_te2$OS)
surv_obj_4


IBS1_4 <- round(IBS(surv_obj_4, sp_matrix = survivalProbs_t_mat1_t2_4, survivalProbs_p4$time[-1]),3)
IBS1_4

#Concordance Index

c_index_4<- round(SurvMetrics::Cindex(surv_obj_4, predicted = medianSurv4),2)
c_index_4

Error_mat_4 <- cbind( IBS1_4, c_index_4)
colnames(Error_mat_4) <- c("IBS_Score", "c_index")
Error_mat_4

#write.table(Error_mat_4 , file = "Error_mat4_with_Univariate_ftrs.txt", sep="\t", quote=F, row.names = T)
############################################################################################################################# 

  }
  
  else if (Model_type == 5) {

########################################## Univariate with Clin features ###############################

    
 #create data frame with selected features (user provided list)
 sel_clin_tr <- as.data.frame(tr_data2[,colnames(tr_data2) %in% c(ftr_list$ID), ])
    
 sel_clin_te<- as.data.frame(te_data2[,colnames(te_data2) %in% c(ftr_list$ID), ])
    
 #add survival information
 sel_clin_tr1  <- cbind(tr_clin1["OS"], tr_clin1["OS_month"],  sel_clin_tr)
 sel_clin_te1  <- cbind(te_clin1["OS"], te_clin1["OS_month"],  sel_clin_te)
    
 head(sel_clin_tr1,2)  
formula5 <- Surv(OS_month,OS) ~ .
#Next, we just need the data argument which in our case is training. We can finally make our first model!
formula5

Mod5<- mtlr(formula = formula5, data =  sel_clin_tr1)
#smallMod <- mtlr(formula = formulaSmall, data = tr_clin3)
#We will print the small model so the output is more compact.
Mod5

save(Mod5, file = "Model_with_univ_and_Clin.RData")
#load saved model
#load("Model_with_PI_and_Clin.RData")
#mtlr outputs the weight matrix for the model – 
#these are the weights corresponding to each feature at each time point (additionally notice that we include the bias weights). 
#The row names correspond to the time point for which the feature weight belongs. 
#If you would like to access these weights, they are saved in the model object as weight_matrix so you can access them using smallMod$weight_matrix.

#We can also plot the weights for a mtlr model. Before we printed the small model but here we will look at the weights for the complete model.

#plot(Mod5, numfeatures = 10)


#jpeg(file="survival_curves_features_weight_PI_clin_ftes.jpeg", units="in", width=10, height=10, res=350)
#plot(Mod3, numfeatures = 10)
#dev.off()


#Model Predictions: Now that we have trained a MTLR model we should make some predictions! 
#This is where our testing set and the predict function will come into play.
#Note that there are a number of predictions we may be interested in acquiring. 
#First, we may want to view the survival curves of our test observations.

sel_clin_te2 <- na.omit(sel_clin_te1)

survCurves5 <- predict(Mod5, sel_clin_te2, type = "survivalcurve")
#survCurves is pretty large so we will look at the first 5 rows/columns.
dim(survCurves5 )
head(survCurves5,2)
survCurves5[1:5,1:5]
colnames(survCurves5)<- c("time_point", rownames(sel_clin_te2))

survCurves5_df <- as.data.frame(survCurves5)
head(survCurves5,2)
#write.table(survCurves5_df , file = "survCurves5_data.txt", sep="\t", quote=F, row.names = F)


#Mean/Median Survival Time: In addition to the entire survival curve one may also be interested in the average survival time. 
#This is again available from the predict function.



#Mean
meanSurv5 <- predict(Mod5, sel_clin_te2, type = "mean_time")
head(meanSurv5)
dim(meanSurv5)

#Median
medianSurv5 <- predict(Mod5, sel_clin_te2, type = "median_time")
head(medianSurv5)
dim(medianSurv5)

meanSurv5_d  <- as.data.frame(meanSurv5)
medianSurv5_d  <- as.data.frame(medianSurv5)

names_5 <- as.data.frame(rownames(sel_clin_te2))
mean_median_surv5_d <- cbind(names_5, meanSurv5_d , medianSurv5_d )
rownames(mean_median_surv5_d)<- c(rownames(sel_clin_te2))

colnames(mean_median_surv5_d) <- c("IDs", "Mean", "Median")
head(mean_median_surv5_d )


#write.table(mean_median_surv5_d , file = "mean_median_time_data5.txt", sep="\t", quote=F, row.names = F)

#Survival Probability at Event Time: The last prediction type supported is acquiring the observations survival probability at the respective event time. However, in order to use this prediction, the event time (whether censored or uncensored) must be included in the features passed into the predict function.

Survival_Probs_event5 <- predict(Mod5, sel_clin_te2, type = "prob_event")
head(Survival_Probs_event5)
dim(Survival_Probs_event5)


survivalProbs_p5 <- predict(Mod5, sel_clin_te2, type = "prob_times")
head(survivalProbs_p5)
dim(survivalProbs_p5)


surv_res5 <- cbind(meanSurv5, medianSurv5, Survival_Probs_event5, sel_clin_te2$OS_month, sel_clin_te2$OS)
surv_res5
colnames(surv_res5) <- c("Mean", "Median_surv", "Survival_Prob_event", "Actual_OS_time", "Event")
rownames(surv_res5) <- rownames(sel_clin_te2)
head(surv_res5)

#write.table( surv_res5 , file = "surv_res5_MTLR_based_on_Univ_and_Clin.txt", sep="\t", quote=F, row.names = T)


#extract prob times at diff times points
survivalProbs_t_mat_5 <-as.matrix(survivalProbs_p5 )
head(survivalProbs_t_mat_5,2)
survivalProbs_t_mat1_5 <- survivalProbs_t_mat_5[,-1]
head(survivalProbs_t_mat1_5) 
survivalProbs_t_mat1_t_5<- t(survivalProbs_t_mat1_5 )
head(survivalProbs_t_mat1_t_5,2)
survivalProbs_t_mat1_t2_5 <- survivalProbs_t_mat1_t_5[,-1]
head(survivalProbs_t_mat1_t2_5)




surv_obj_5 <- Surv(sel_clin_te2$OS_month,sel_clin_te2$OS)
surv_obj_5

IBS1_5<- round(IBS(surv_obj_5, sp_matrix = survivalProbs_t_mat1_t2_5, survivalProbs_p5$time[-1]),3)
IBS1_5

#Concordance Index

c_index_5<- round(SurvMetrics::Cindex(surv_obj_5, predicted = medianSurv5),2)
c_index_5


Error_mat_5 <- cbind(IBS1_5, c_index_5)
colnames(Error_mat_5) <- c("IBS_Score", "c_index")
Error_mat_5

#write.table(Error_mat_5 , file = "Error_mat5_with_Univ_and Clin.txt", sep="\t", quote=F, row.names = T)

  }
  
  if(Model_type == 1){
    
    
    survCurves_df <- survCurves1_df
    mean_median_surv_d  <- mean_median_surv1_d 
    Error_mat <- Error_mat_1
    surv_res <- surv_res1
  
    
  }
  
  
  if(Model_type == 2){
    
    
    survCurves_df <- survCurves2_df
    mean_median_surv_d  <- mean_median_surv2_d 
    Error_mat <- Error_mat_2
    surv_res <- surv_res2
    
    
  }
  
  
  if(Model_type == 3){
    
    
    survCurves_df <- survCurves3_df
    mean_median_surv_d  <- mean_median_surv3_d 
    Error_mat <- Error_mat_3
    surv_res <- surv_res3
    
    
  }
  
  
  if(Model_type == 4){
    
    
    survCurves_df <- survCurves4_df
    mean_median_surv_d  <- mean_median_surv4_d 
    Error_mat <- Error_mat_4
    surv_res <- surv_res4
    
    
  }
  
  
  
  if(Model_type == 5){
    
    
    survCurves_df <- survCurves5_df
    mean_median_surv_d  <- mean_median_surv5_d 
    Error_mat <- Error_mat_5
    surv_res <- surv_res5
    
    
  }
  

  write.table(survCurves_df , file = "survCurves_data.txt", sep="\t", quote=F, row.names = F)
  write.table(mean_median_surv_d , file = "mean_median_survival_time_data.txt", sep="\t", quote=F, row.names = F)
  write.table( surv_res , file = "survival_result_based_on_MTLR.txt", sep="\t", quote=F, row.names = T)
  write.table(Error_mat , file = "Error_mat_for_Model.txt", sep="\t", quote=F, row.names = T)
  
  
}

  
