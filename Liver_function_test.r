## Hui Xin Sim
## MovieLens Project 
## HarvardX: PH125.9x - Capstone Project
## https://github.com/zenasimlab

#############################################
# Role of Liver Function Test in Diagnosis of Liver Diseases Code
##############################################

###Introduction###

## Dataset ##

knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(caret)


## Retrieve the dataset
columnNames <- c('age', # Age of the patient 
                'sex', # Sex of the patient 
                'tb', # Total Bilirubin
                'db', # Direct Bilirubin 
                'alkphos', # Alkaline Phosphotase
                'sgpt', # Alamine Aminotransferase
                'sgot', # Aspartate Aminotransferase
                'tp', # Total Protiens
                'alb', # Albumin
                'ag', # Ratio	Albumin and Globulin Ratio 
                'outcome') # Selector field used to split the data into two sets

fullData <- read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/00225/Indian%20Liver%20Patient%20Dataset%20(ILPD).csv",
                       sep=',',
                       header=FALSE,
                       col.names=columnNames)

## Data pre-processing and exploratory analysis

## Remove any incomplete rows
fullData <- subset(fullData, complete.cases(fullData))

## Format it in a more human-understandable way

fullData <- fullData %>% 
  mutate(outcome = as.character(outcome)) %>% 
  mutate(outcome = replace(outcome, outcome == '1', 'Care')) %>%
  mutate(outcome = replace(outcome, outcome == '2', 'Control')) %>%
  mutate(outcome = as.factor(outcome))

rm(columnNames)
head(fullData)

## Create train and test sets
set.seed(1)
trainIndex <- createDataPartition(fullData$outcome, p=.7, list=FALSE)
train <- fullData[trainIndex,]
test <- fullData[-trainIndex,]
rm(trainIndex)

## Age distribution
train %>% 
  ggplot(aes(x = age)) + 
  geom_histogram(binwidth = 20) + 
  theme_minimal() +
  facet_grid(~ outcome)

## Sex distribution
train %>% 
ggplot(aes(x = sex)) + 
  geom_bar() + 
  theme_minimal() +
  facet_grid(~ outcome)

## Bilirubin level
## Total bilirubin
train %>% 
  ggplot(aes(x = tb)) + 
  geom_histogram(binwidth = 5) + 
  theme_minimal() +
  facet_grid(~ outcome)

## Direct bilirubin
train %>% 
  ggplot(aes(x = db)) + 
  geom_histogram(binwidth = 1) + 
  theme_minimal() +
  facet_grid(~ outcome)

## Enzymes

## Alkaline phosphatase
train %>% 
  ggplot(aes(x = alkphos)) + 
  geom_histogram(binwidth = 200) + 
  theme_minimal() +
  facet_grid(~ outcome)

## Alanine aminotransferase
train %>% 
  ggplot(aes(x = sgpt)) + 
  geom_histogram(binwidth = 200) + 
  theme_minimal() +
  facet_grid(~ outcome)

## Aspartate aminotransferase
train %>% 
  ggplot(aes(x = sgpt)) + 
  geom_histogram(binwidth = 200) + 
  theme_minimal() +
  facet_grid(~ outcome)

## Alkaline phosphatase vs Total bilirubin
train %>% 
  ggplot(aes(x = alkphos, y = tb, shape = outcome, color = outcome)) + 
  geom_point() +
  scale_y_log10() + 
  scale_x_log10() +
  geom_vline(xintercept = 700) +
  geom_hline(yintercept = 7.5) +
  theme_minimal()

## Protein
## Albumin
train %>% 
  ggplot(aes(x = alb)) + 
  geom_histogram(binwidth = 0.5) + 
  theme_minimal() +
  facet_grid(~ outcome)

## Total Protein
train %>% 
  ggplot(aes(x = tp)) + 
  geom_histogram(binwidth = 0.5) + 
  theme_minimal() +
  facet_grid(~ outcome)

## Ratio between albumin & globulin
subset(train, !is.na(ag)) %>% 
  ggplot(aes(x = ag)) + 
  geom_histogram(binwidth = 0.1) + 
  theme_minimal() +
  facet_grid(~ outcome)

## Data preparation
## Correlated Predictors

cor(subset(train, select = -c(sex, outcome)))

## Bilirubin
train %>% 
  ggplot(aes(x = tb, y = db)) + 
  geom_point() + 
  theme_minimal() +
  facet_grid(~ outcome)

## Drop direct bilirubin
train <- train %>% subset(select = -c(db))
test <- test %>% subset(select = -c(db))

## Aminotransferases
train %>% 
  ggplot(aes(x = sgot, y = sgpt, color = outcome, shape = outcome)) + 
  geom_point() + 
  theme_minimal()

## Drop Alanine Aminotransferase
train <- train %>% subset(select = -c(sgpt))
test <- test %>% subset(select = -c(sgpt))

## Protein
train %>% 
  ggplot(aes(x = tp, y = alb, color = outcome, shape = outcome)) + 
  geom_point() + 
  theme_minimal()

## Drop albumin
train <- train %>% subset(select = -c(alb))
test <- test %>% subset(select = -c(alb))

## Machine Learning Methods
results <- data.frame(Model = character(), 
                      Accuracy = double(), 
                      Sensitivity = double(), 
                      Specificity = double(), 
                      stringsAsFactors = FALSE)

## Naive Bayes
nb_model = train(outcome ~ ., data = train, method = "nb")
predictions = predict(nb_model, newdata = test)
confusionMatrix <- confusionMatrix(predictions, test$outcome)
results[nrow(results) + 1, ] <- c(as.character('Naive Bayes (nb)'), 
                                  confusionMatrix$overall['Accuracy'],  
                                  confusionMatrix$byClass['Sensitivity'], 
                                  confusionMatrix$byClass['Specificity'])
rm(nb_model, predictions)
confusionMatrix

## Linear Classifier
lc_model = train(outcome ~ ., data = train, method = "glmboost")
predictions = predict(lc_model, newdata = test)
confusionMatrix <- confusionMatrix(predictions, test$outcome)
results[nrow(results) + 1, ] <- c(as.character('Linear Classifier (glmboost)'), 
                                  confusionMatrix$overall['Accuracy'],  
                                  confusionMatrix$byClass['Sensitivity'], 
                                  confusionMatrix$byClass['Specificity'])
rm(lc_model, predictions)
confusionMatrix

## Calculating Prevalence using Naive Bayes Model and Linear Classifier Model
nb_model = train(outcome ~ ., data = train, method = "nb")
predictions = predict(nb_model, newdata = test)
confusionMatrix <- confusionMatrix(predictions, test$outcome, prevalence = 0.06)
rm(nb_model, predictions)
confusionMatrix

lc_model = train(outcome ~ ., data = train, method = "glmboost")
predictions = predict(lc_model, newdata = test)
confusionMatrix <- confusionMatrix(predictions, test$outcome, prevalence = 0.06)
rm(lc_model, predictions)
confusionMatrix

## Logistic Regression
lr_model = train(outcome ~ ., data = train, method = "bayesglm")
predictions = predict(lr_model, newdata = test)
confusionMatrix <- confusionMatrix(predictions, test$outcome, prevalence = 0.06)
results[nrow(results) + 1, ] <- c(as.character('Logistic Regression (bayesglm)'), 
                                  confusionMatrix$overall['Accuracy'],  
                                  confusionMatrix$byClass['Sensitivity'], 
                                  confusionMatrix$byClass['Specificity'])
rm(lr_model, predictions)
confusionMatrix

## K-Nearest Neighbours
## Value of K is a tuning parameter
knn_model = train(outcome ~ ., data = train, method = "knn", preProcess=c('knnImpute'))
knn_model

predictions = predict(knn_model, newdata = test)
confusionMatrix <- confusionMatrix(predictions, test$outcome, prevalence = 0.06)
results[nrow(results) + 1, ] <- c(as.character('K-nearest neighbours (knn)'), 
                                  confusionMatrix$overall['Accuracy'],  
                                  confusionMatrix$byClass['Sensitivity'], 
                                  confusionMatrix$byClass['Specificity'])
rm(knn_model, predictions)
confusionMatrix

## Random Forest
rf_model = train(outcome ~ ., data = train, method = "rf")
predictions = predict(rf_model, newdata = test)
confusionMatrix <- confusionMatrix(predictions, test$outcome, prevalence = 0.06)
results[nrow(results) + 1, ] <- c(as.character('Random Forest (rf)'), 
                                  confusionMatrix$overall['Accuracy'],  
                                  confusionMatrix$byClass['Sensitivity'], 
                                  confusionMatrix$byClass['Specificity'])
rm(rf_model, predictions)
confusionMatrix


rm(confusionMatrix)


## Filter the dataset with liver enzymes higher than 2 times of upper limit of normal range
## Creating another sets of training and test sets

set.seed(1, sample.kind = "Rounding") 
fullData2<-fullData%>%filter(alkphos>266&sgpt>112&sgot>70)
trainIndex2 <- createDataPartition(fullData2$outcome,times=1,p=0.7, list=FALSE)
train2 <- fullData2[trainIndex2,]
test2 <- fullData2[-trainIndex2,]
rm(trainIndex2)

head(fullData2)

nrow(fullData2)

## Machine Learning Methods
## Naive Bayes 

nb_model2 = train(outcome ~ ., data = train2, method = "nb")
predictions2 = predict(nb_model2, newdata = test2)
confusionMatrix2 <- confusionMatrix(predictions2, test2$outcome)
rm(nb_model2, predictions2)
confusionMatrix2

##Linear Classifier

lc_model2 = train(outcome ~ ., data = train2, method = "glmboost")
predictions2 = predict(lc_model2, newdata = test2)
confusionMatrix2 <- confusionMatrix(predictions2, test2$outcome)
rm(lc_model2, predictions2)
confusionMatrix2

##Logistic regression

lr_model2 = train(outcome ~ ., data = train2, method = "bayesglm")
predictions2 = predict(lr_model2, newdata = test2)
confusionMatrix2 <- confusionMatrix(predictions2, test2$outcome, prevalence = 0.06)
rm(lr_model2, predictions2)
confusionMatrix2

## knn model

knn_model2 = train(outcome ~ ., data = train2, method = "knn", preProcess=c('knnImpute'))
knn_model2

predictions2 = predict(knn_model2, newdata = test2)
confusionMatrix2 <- confusionMatrix(predictions2, test2$outcome, prevalence = 0.06)
rm(knn_model2, predictions2)
confusionMatrix2

## Random Forest

rf_model2 = train(outcome ~ ., data = train2, method = "rf")
predictions2 = predict(rf_model2, newdata = test2)
confusionMatrix2 <- confusionMatrix(predictions2, test2$outcome, prevalence = 0.06)
rm(rf_model2, predictions2)
confusionMatrix2

## Filter the dataset with liver enzymes raised but less than 2 times of upper limit of normal range
## Creating another sets of training and test sets

set.seed(1, sample.kind = "Rounding") 
fullData3<-fullData%>%filter(alkphos%in%(133:266)|sgpt%in%(56:112)|sgot%in%(35:70))
trainIndex3 <- createDataPartition(fullData3$outcome, times=1, p=0.7, list=FALSE)
train3 <- fullData3[trainIndex3,]
test3 <- fullData3[-trainIndex3,]
rm(trainIndex3)

head(fullData3)

nrow(fullData3)

## Machine Learning Methods
## Naive Bayes

nb_model3 = train(outcome ~ ., data = train3, method = "nb")
predictions3 = predict(nb_model3, newdata = test3)
confusionMatrix3 <- confusionMatrix(predictions3, test3$outcome)
rm(nb_model3, predictions3)
confusionMatrix3

##Linear Classifier

lc_model3 = train(outcome ~ ., data = train3, method = "glmboost")
predictions3 = predict(lc_model3, newdata = test3)
confusionMatrix3 <- confusionMatrix(predictions3, test3$outcome)
rm(lc_model3, predictions3)
confusionMatrix3

## Logistic regression

lr_model3 = train(outcome ~ ., data = train3, method = "bayesglm")
predictions3 = predict(lr_model3, newdata = test3)
confusionMatrix3 <- confusionMatrix(predictions3, test3$outcome, prevalence = 0.06)
rm(lr_model3, predictions3)
confusionMatrix3

## K-Nearest Neighbours

knn_model3 = train(outcome ~ ., data = train3, method = "knn", preProcess=c('knnImpute'))
knn_model3

predictions3 = predict(knn_model3, newdata = test3)
confusionMatrix3 <- confusionMatrix(predictions3, test3$outcome, prevalence = 0.06)
rm(knn_model3, predictions3)
confusionMatrix3

## Random Forest

rf_model3 = train(outcome ~ ., data = train3, method = "rf")
predictions3 = predict(rf_model3, newdata = test3)
confusionMatrix3 <- confusionMatrix(predictions3, test3$outcome, prevalence = 0.06)
rm(rf_model3, predictions3)
confusionMatrix3

## Filter the dataset with liver enzymes which are within normal range
## Creating another sets of training and test sets

set.seed(1, sample.kind = "Rounding") 
fullData4<-fullData%>%filter(alkphos%in%(41:133)&sgpt%in%(7:56)&sgot%in%(0:35))
trainIndex4 <- createDataPartition(fullData4$outcome,times=1,p=0.7, list=FALSE)
train4 <- fullData4[trainIndex4,]
test4 <- fullData4[-trainIndex4,]
rm(trainIndex4)

head(fullData4)

nrow(fullData4)

## Naive Bayes

nb_model4 = train(outcome ~ ., data = train4, method = "nb")
predictions4 = predict(nb_model4, newdata = test4)
confusionMatrix4 <- confusionMatrix(predictions4, test4$outcome)
rm(nb_model4, predictions4)
confusionMatrix4

## Linear regression

lc_model4 = train(outcome ~ ., data = train4, method = "glmboost")
predictions4 = predict(lc_model4, newdata = test4)
confusionMatrix3 <- confusionMatrix(predictions4, test4$outcome)
rm(lc_model4, predictions4)
confusionMatrix4

##Logistic regression

lr_model4 = train(outcome ~ ., data = train4, method = "bayesglm")
predictions4 = predict(lr_model4, newdata = test4)
confusionMatrix4 <- confusionMatrix(predictions4, test4$outcome, prevalence = 0.06)
rm(lr_model4, predictions4)
confusionMatrix4

## knn model
knn_model4 = train(outcome ~ ., data = train4, method = "knn", preProcess=c('knnImpute'))
knn_model4

predictions4 = predict(knn_model4, newdata = test4)
confusionMatrix4 <- confusionMatrix(predictions4, test4$outcome)
rm(knn_model4, predictions4)
confusionMatrix4

## Random Forest

rf_model4 = train(outcome ~ ., data = train4, method = "rf")
predictions4 = predict(rf_model4, newdata = test4)
confusionMatrix4 <- confusionMatrix(predictions4, test4$outcome, prevalence = 0.06)
rm(rf_model4, predictions4)
   
## Results
results %>% arrange(Accuracy)  
