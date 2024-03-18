#' Subsequently, PI score will be calculated by multiplying expression of genes with their beta coeff values. 
#' 
#' #' @param train_data :args1 - training data (Patients data with clinical and gene expression, where samples are in rows and features/genes are in columns)
#' @param test_data_data :args2 - training data (Patients data with clinical and gene expression, where samples are in rows and features/genes are in columns)
#' @param  nfolds : args3 - Number of folds for cross in cvglmnet model to select top features
#' @param col_num :args4 - column number in data at where clinical info ends
#' @param surv_time :arg5 - name of column which contain survival time (in days) information
#' @param surv_event :arg6 - name of column which contain survival eventinformation
#' @param train_PI_data :args7 - name of output file containing LASSO selected genes and PI score data for training data
#' @param test_PI_data :args8 - name of output file containing LASSO selected genes and PI score data for test data

#' @examples
#' SurvPredPipe::Lasso_PI_scores_f(train_data="Train_Norm_data.txt",test_data="Test_Norm_data.txt", nfolds=5, col_num=21, surv_time="OS_month", surv_event="OS", train_PI_data="Train_PI_data.txt", test_PI_data="Test_PI_data.txt" )
#' Usage: Lasso_PI_scores_f(train_data, test_data, nfolds, col_num, surv_time, surv_event, train_PI_data, test_PI_data)
#' @export



Lasso_PI_scores_f <- function(train_data, test_data, nfolds, col_num, surv_time, surv_event, train_PI_data, test_PI_data) {
  
  # Check if any input variable is empty
  if (length(train_data) == 0 ||  length(test_data) == 0|| length(nfolds) == 0 ||length(col_num) == 0 ||length(surv_time) == 0 ||length(surv_event) == 0 ||length(train_PI_data) == 0 ||length(test_PI_data) == 0 ) {
    stop("Error: Empty input variable detected.")
  }
  
  # Check if any input variable is missing
  if (any(is.na(train_data)) || any(is.na(test_data)) || any(is.na(nfolds)) || any(is.na(col_num)) || any(is.na(surv_time)) || any(is.na(surv_event)) || any(is.na(train_PI_data)) || any(is.na(test_PI_data))    ) {
    stop("Error: Missing values in input variables.")
  }
  

#load data
tr_data1 <- read.table(train_data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
te_data1 <- read.table(test_data, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)


# rename survival time and event column name
colnames(tr_data1)[colnames(tr_data1) == surv_time] <- "OS_month"
colnames(tr_data1)[colnames(tr_data1) == surv_event] <- "OS"


# rename survival time and event column name
colnames(te_data1)[colnames(te_data1) == surv_time] <- "OS_month"
colnames(te_data1)[colnames(te_data1) == surv_event] <- "OS"

########################################### LASSO COX ##########################################
surv_object <- Surv(time = tr_data1$OS_month, event = tr_data1$OS)

set.seed(7)
cvfit1 <- cv.glmnet(
  #as.matrix(tr_data1[21:ncol(tr_data1)]),
  as.matrix(tr_data1[col_num:ncol(tr_data1)]),
  surv_object, # create survival object from the data
  family = "cox", # specify Cox PH model
  type.measure = "C",
  #nfolds = 5,
  nfolds =  nfolds,
  alpha = 1, # lasso: alpha = 1; ridge: alpha=0
  maxit = 1000)

lambda_min <- cvfit1$lambda.min
lambda_min


est.coef = coef(cvfit1, s = cvfit1$lambda.min) # returns the p length coefficient vector
head(est.coef)
# of the solution corresponding to lambda
est.coef1 <- as.numeric(est.coef)
active.k = which(est.coef1 != 0)
#active.k = which(!is.na(est.coef) & est.coef != 0)
active.k

#Extract coefficinet values
active.k.vals = est.coef[active.k]
active.k.vals

key_variables <- as.data.frame(est.coef[est.coef[,1]!=0,])
colnames(key_variables) <- c("coeff")
key_variables <- round(key_variables,3)
dim(key_variables)

#write.table(key_variables,file="Lasso_key_variables.txt",row.names=T,col.names=T,sep = '\t', quote = F)
write.table(cbind("ID"=rownames(key_variables), key_variables),file="Train_Lasso_key_variables.txt",sep="\t",quote=F, row.names=F)


jpeg(file="Train_Cox_Lasso_Regression_lamda_plot.jpeg", units="in", width=10, height=10, res=350)
plot(cvfit1)
dev.off()


  
svg(file="Train_Cox_Lasso_Regression_lamda_plot.svg")
plot(cvfit1) 
dev.off()
################## Create PI Index for training data #######################

sel_features2 <- read.table("Train_Lasso_key_variables.txt",header =TRUE, sep = "\t", row.names=1, check.names = FALSE)

head(sel_features2)
dim(sel_features2)


#training data preparation with selected features 
#sel_train <- as.data.frame(data1[,colnames(data1) %in% c(sel_features_results$ID), ])
sel_train2 <- as.data.frame(tr_data1[,colnames(tr_data1) %in% c(row.names(sel_features2)), ])
head(sel_train2,2)
dim(sel_train2)

######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
train_feature_mat2<- cbind(tr_data1["OS_month"],tr_data1["OS"],sel_train2)
head(train_feature_mat2,2)
dim(train_feature_mat2)


######################## PI Index ########################
dim(train_feature_mat2)
dim(sel_features2)
#Create prognostic index 

tr_PI <- train_feature_mat2[3:ncol(train_feature_mat2)]

E= length(tr_PI)
E

head(sel_features2)
sel_features2[2,1]


PI_tr = 0
for(i in seq(from=1, to=E ,by=1))
{
  PI_tr= PI_tr+((tr_PI[,i])*(sel_features2[i,1]))
  #PI_tr= PI_tr+((tr[,i])* 1)
}



# add PI as new column to the data
tr_PI$PI<-PI_tr
head(tr_PI)


#tr_PI <- as.data.frame(train_feature_mat2$PI)
#rownames(tr_PI) <- rownames(train_feature_mat2)
#colnames(tr_PI) <- c("PI")
#head(tr_PI)

train_PI <- cbind(train_feature_mat2["OS"], train_feature_mat2["OS_month"], tr_PI)
head(train_PI)
dim(train_PI)
#save selected data with PI value
write.table(cbind("ID"=rownames(train_PI ), train_PI),file=train_PI_data,sep="\t",quote=F, row.names=F)


#######################################################

###################### PI for Test data ###############

sel_test <- as.data.frame(te_data1[,colnames(te_data1) %in% c(row.names(sel_features2)), ])
head(sel_test,2)
dim(sel_test)

######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
test_feature_mat2<- cbind(te_data1["OS_month"],te_data1["OS"],sel_test)
head(test_feature_mat2,2)
dim(test_feature_mat2)


dim(test_feature_mat2)
dim(sel_features2)
#Create prognostic index 

te_PI <- test_feature_mat2[3:ncol(test_feature_mat2)]

E= length(te_PI)
E

head(sel_features2)
sel_features2[2,1]

PI_te = 0
for(i in seq(from=1, to=E ,by=1))
{
  PI_te= PI_te+((te_PI[,i])*(sel_features2[i,1]))
  #PI_tr= PI_tr+((tr[,i])* 1)
}


# add PI as new column to the data
te_PI$PI<-PI_te



test_PI <- as.data.frame(te_PI$PI)
rownames(test_PI) <- rownames(te_PI)
colnames(test_PI) <- c("PI")

test_PI <- cbind(test_feature_mat2["OS"], test_feature_mat2["OS_month"], te_PI)
#save selected data with PI value
head(test_PI)

#write into a file
#write.table(cbind("ID"=rownames(test_PI ), test_PI),file="test_Lasso_genes_with_PI.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(test_PI ), test_PI),file=test_PI_data,sep="\t",quote=F, row.names=F)



}


