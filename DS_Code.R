#Patrick Keppler and Matt Mietchen, 2019

library(dplyr)
library(ggplot2)
library(glmnet)

#rstudioapi is used to set working directory to
#the location of the r script in Rstudio. 
#if not installed, uncomment the next line:
#install.packages(rstudioapi)
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##########Load Model Data##########
model_data = read.csv("emulate_main.csv")


##########Exploratory Data Analysis##########

#plot(model_data$sigma, model_data$psi) #not very helpful
#plot(model_data$acquisitions, model_data$psi) #potentially helpful, too much noise
#get summary statistics
summary(model_data$acquisitions)
#show boxplot with discrete predictor (bedsize)
bedsize_factor = factor(model_data$bedsize)
ggplot(NULL,aes(x=bedsize_factor, y=model_data$acquisitions)) + geom_boxplot() + labs(title="Boxplot of Bed Size vs Acquisitions", x="Bed Size", y="Acquisitons")

#omit the lowest and highest 5% of acquisitions 
acquisition_hundredths = quantile(model_data$acquisitions, probs = seq(0,1, by=0.01))
reduced_data = model_data %>%
  filter(acquisitions < acquisition_hundredths[96])

summary(reduced_data$acquisitions)

#show boxplot with problematic values omitted
reduced_bedsize_factor = factor(reduced_data$bedsize)
ggplot(NULL,aes(x=reduced_bedsize_factor, y=reduced_data$acquisitions)) + geom_boxplot() + labs(title="Reduced Boxplot of Bed Size vs Acquisitions", x="Bed Size", y="Acquisitons")

##########Backwards Feature Selection##########

#create linear model with all predictors
full_lm = lm(acquisitions ~ ., data = reduced_data)
summary(full_lm)

#create a linear model with all predictors and all pairwise interactions
interaction_lm = lm(acquisitions~(.)^2, data = reduced_data)
summary(interaction_lm)

#remove insignificant predictors/interactions
reduced_lm = lm(acquisitions~(.-iota_D-tau_N-tau_D-gamma)^2 + sigma * (gamma + tau_N)
                + gamma * (theta + tau_D + bedsize) + theta * iota_D, data=reduced_data)
summary(reduced_lm)

##########LASSO##########

##"Interaction version", includes interactions between features##

#Feature matrix containing interactions must be generated before creating the model
feature_interactions = model.matrix(acquisitions~(.^2),reduced_data)[,-1]
response = reduced_data$acquisitions

#Start with cross-validation method to determine the best value for lambda
cv_lasso_model = cv.glmnet(feature_interactions,response, lambda = seq(0.01, 1, by=0.01), nfolds = 5)
lambda =  cv_lasso_model$lambda.min
lambda

#Construct model
lasso_model = glmnet(feature_interactions,response)
#show coeficients
beta_hat = coef(lasso_model,s = lambda)
beta_hat

###########Principal Components Analysis###########

#pca is somewhat inefficient to run on entire data set, run on smaller sample
#when run on full data set, the results were relatively unaffected
pca_sample = sample(nrow(reduced_data), 100)

#run pca
pca.out = prcomp(reduced_data[pca_sample,], scale=TRUE)

#produce a biplot to show the first two principal components
biplot(pca.out, scale=TRUE, cex=c(0.1,1))

#Compute Proportion of Variance Explained (PVE)
pca.out = prcomp(reduced_data, scale=TRUE)
pca.out$sdev
pca.var=pca.out$sdev^2
pca.var
pve=pca.var/sum(pca.var)
pve

#Show scree plots for PCA performance
plot(pve, xlab="Principal Component", ylab="Proportion of Variance Explained", ylim=c(0,1),type='b')
plot(cumsum(pve), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", ylim=c(0,1),type='b')


##########Model Evaluation##########

##R-Squared Comparision##

#Adjusted R-Squared values in summary
summary(full_lm)
summary(interaction_lm)
summary(reduced_lm)

##AIC scores##
#setup
full_glm = glm(acquisitions ~ ., data = reduced_data)
interaction_glm = glm(acquisitions~(.)^2, data = reduced_data)
reduced_glm = glm(acquisitions~(.-iota_D-tau_N-tau_D-gamma)^2 + sigma * (gamma + tau_N)
                  + gamma * (theta + tau_D + bedsize) + theta * iota_D, data=reduced_data)

#AIC in summary
summary(full_glm)
summary(interaction_glm)
summary(reduced_glm)


##Root Mean Squared Error (RMSE)##

#get data for testing from other file
test_data = read.csv("emulate_main2.csv")

#omit the lowest and highest 5% of acquisitions 
acquisition_hundredths = quantile(test_data$acquisitions, probs = seq(0,1, by=0.01))
reduced_test_data = model_data %>%
  filter(acquisitions < acquisition_hundredths[96])

#Full model
full_model_pred = predict(full_lm, newdata=reduced_test_data)
sqrt(mean((full_model_pred-reduced_test_data$acquisitions)^2))

#Interaction model
interaction_model_pred = predict(interaction_lm, newdata=reduced_test_data)
sqrt(mean((interaction_model_pred-reduced_test_data$acquisitions)^2))

#Reduced model
reduced_model_pred = predict(reduced_lm, newdata=reduced_test_data)
sqrt(mean((reduced_model_pred-reduced_test_data$acquisitions)^2))

#LASSO model
test_data_interactions = model.matrix(acquisitions~.*.,reduced_test_data)[,-1]
lasso_model_pred = predict(lasso_model, newx = test_data_interactions, s=0.01)
sqrt(mean((lasso_model_pred-reduced_test_data$acquisitions)^2))

###########Psi Approximation###############
#use interaction model
psi_lm = lm(psi~.*., data = reduced_data)
#note a very poor adjusted R-squared score
summary(psi_lm)

##LASSO##

#Feature matrix containing interactions must be generated before creating the model
psi_feature_interactions = model.matrix(psi~(.^2),reduced_data)[,-1]
psi_response = reduced_data$psi

#Start with cross-validation method to determine the best value for lambda
psi_cv_lasso_model = cv.glmnet(psi_feature_interactions, psi_response, lambda = seq(0.01, 1, by=0.01), nfolds = 5)
psi_lambda =  psi_cv_lasso_model$lambda.min
psi_lambda

#Construct lasso model
psi_lasso_model = glmnet(psi_feature_interactions,psi_response)
#show coeficients
psi_beta_hat = coef(psi_lasso_model,s = psi_lambda)


#Interaction model RMSE
psi_interaction_model_pred = predict(psi_lm, newdata=reduced_test_data)
sqrt(mean((psi_interaction_model_pred-reduced_test_data$psi)^2))

#LASSO RMSE
psi_test_data_interactions = model.matrix(psi~.*.,reduced_test_data)[,-1]
psi_lasso_model_pred = predict(psi_lasso_model, newx = psi_test_data_interactions, s=psi_lambda)
sqrt(mean((psi_lasso_model_pred-reduced_test_data$psi)^2))
