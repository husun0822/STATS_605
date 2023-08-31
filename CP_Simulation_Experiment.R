# ----- Conformal Prediction Simulation Experiment ----- #
library(dplyr)
library(tidyr)
library(EnvStats)
library(VGAM)

# ----- Study 1: Mis-specified Model (with Gaussian error) ----- #
full_result <- NULL
for (seed in 1:20){
  set.seed(seed)
  print(seed)
  n <- 1500
  x1 <- rnorm(n = n, mean = 0); x2 <- x1^2; x3 <- x1^3
  e_gauss <- rnorm(n = n); e_t <- rt(n = n, df = 2, ncp = 0)
  y1 <- 3 * x1 - 4 * x2 + 5 * x3 + e_gauss
  # y2 <- 3 * x1 - 4 * x2 + 5 * x3 + e_t
  
  # train-test split 
  n_train <- 1:1000; n_test <- 1001:n
  x1_train <- x1[n_train]; x2_train <- x2[n_train]; x3_train <- x3[n_train]; y1_train <- y1[n_train]; y2_train <- y2[n_train]
  x1_test <- x1[n_test]; x2_test <- x2[n_test]; x3_test <- x3[n_test]; y1_test <- y1[n_test]; y2_test <- y2[n_test]
  
  model1 <- lm(y1_train ~ x1_train)
  model1_pred <- predict(model1, data.frame(x1_train = x1_test), interval = "predict")
  model1_cov <- sum((y1_test >= model1_pred[,2]) & (y1_test <= model1_pred[,3])) / length(y1_test)
  model1_CI_bw <- mean(model1_pred[,3] - model1_pred[,2])
  
  model2 <- lm(y1_train ~ x1_train + x2_train)
  model2_pred <- predict(model2, data.frame(x1_train = x1_test, x2_train = x2_test), interval = "predict")
  model2_cov <- sum((y1_test >= model2_pred[,2]) & (y1_test <= model2_pred[,3])) / length(y1_test)
  model2_CI_bw <- mean(model2_pred[,3] - model2_pred[,2])
  
  model3 <- lm(y1_train ~ x1_train + x2_train + x3_train)
  model3_pred <- predict(model3, data.frame(x1_train = x1_test, x2_train = x2_test, x3_train = x3_test), interval = "predict")
  model3_cov <- sum((y1_test >= model3_pred[,2]) & (y1_test <= model3_pred[,3])) / length(y1_test)
  model3_CI_bw <- mean(model3_pred[,3] - model3_pred[,2])
  
  
  # ----- Study 1: Mis-specified Model with Split Conformal Prediction ----- #
  x1_train_cp <- x1_train[1:500]; x2_train_cp <- x2_train[1:500]; x3_train_cp <- x3_train[1:500]; y1_train_cp <- y1_train[1:500]
  x1_cali_cp <- x1_train[501:1000]; x2_cali_cp <- x2_train[501:1000]; x3_cali_cp <- x3_train[501:1000]; y1_cali_cp <- y1_train[501:1000]
  
  model1_cp <- lm(y1_train_cp ~ x1_train_cp)
  model1_cali_pred <- predict(model1_cp, data.frame(x1_train_cp = x1_cali_cp))
  model1_test_pred <- predict(model1_cp, data.frame(x1_train_cp = x1_test))
  ncs_cali <- abs(y1_cali_cp - model1_cali_pred)
  d <- sort(ncs_cali)[ceiling(501 * 0.95)]
  lb <- model1_test_pred - d; ub <- model1_test_pred + d
  model1_cp_cov <- sum((y1_test >= lb) & (y1_test <= ub)) / length(y1_test)
  model1_cp_CI_bw <- 2 * d
  
  model2_cp <- lm(y1_train_cp ~ x1_train_cp + x2_train_cp)
  model2_cali_pred <- predict(model2_cp, data.frame(x1_train_cp = x1_cali_cp, x2_train_cp = x2_cali_cp))
  model2_test_pred <- predict(model2_cp, data.frame(x1_train_cp = x1_test, x2_train_cp = x2_test))
  ncs_cali <- abs(y1_cali_cp - model2_cali_pred)
  d <- sort(ncs_cali)[ceiling(501 * 0.95)]
  lb <- model2_test_pred - d; ub <- model2_test_pred + d
  model2_cp_cov <- sum((y1_test >= lb) & (y1_test <= ub)) / length(y1_test)
  model2_cp_CI_bw <- 2 * d
  
  model3_cp <- lm(y1_train_cp ~ x1_train_cp + x2_train_cp + x3_train_cp)
  model3_cali_pred <- predict(model3_cp, data.frame(x1_train_cp = x1_cali_cp, x2_train_cp = x2_cali_cp, x3_train_cp = x3_cali_cp))
  model3_test_pred <- predict(model3_cp, data.frame(x1_train_cp = x1_test, x2_train_cp = x2_test, x3_train_cp = x3_test))
  ncs_cali <- abs(y1_cali_cp - model3_cali_pred)
  d <- sort(ncs_cali)[ceiling(501 * 0.95)]
  lb <- model3_test_pred - d; ub <- model3_test_pred + d
  model3_cp_cov <- sum((y1_test >= lb) & (y1_test <= ub)) / length(y1_test)
  model3_cp_CI_bw <- 2 * d
  
  full_result <- rbind(full_result, c(model1_cov, model2_cov, model3_cov, model1_cp_cov, model2_cp_cov, model3_cp_cov,
                                      model1_CI_bw, model2_CI_bw, model3_CI_bw, model1_cp_CI_bw, model2_cp_CI_bw, model3_cp_CI_bw))
}

full_result <- as.data.frame(full_result)
colnames(full_result) <- c("M1_Cov", "M2_Cov", "M3_Cov", "M1_CP_Cov", "M2_CP_Cov", "M3_CP_Cov",
                           "M1_CIwidth", "M2_CIwidth", "M3_CIwidth", "M1_CP_CIwidth", "M2_CP_CIwidth", "M3_CP_CIwidth")

full_result_mean <- full_result %>%
  summarise_all(mean)

full_result_sd <- full_result %>%
  summarise_all(sd)
full_result_sd <- full_result_sd * 1.96


# ----- Study 2: Non-Gaussian Error ----- #
full_result <- NULL
for (seed in 1:20){
  set.seed(seed)
  print(seed)
  n <- 1500
  x1 <- rnorm(n = n, mean = 0); x2 <- x1^2; x3 <- x1^3
  e_gauss <- rnorm(n = n); e_laplace <- rt(n = n, df = 2); e_pareto <- EnvStats::rpareto(n = n, location = 0.5, shape = 2)
  y1 <- 3 * x1 - 4 * x2 + 5 * x3 + e_gauss
  y2 <- 3 * x1 - 4 * x2 + 5 * x3 + e_laplace
  y3 <- 3 * x1 - 4 * x2 + 5 * x3 + e_pareto
  
  # train-test split 
  n_train <- 1:1000; n_test <- 1001:n
  x1_train <- x1[n_train]; x2_train <- x2[n_train]; x3_train <- x3[n_train]; y1_train <- y1[n_train]; y2_train <- y2[n_train]; y3_train <- y3[n_train]
  x1_test <- x1[n_test]; x2_test <- x2[n_test]; x3_test <- x3[n_test]; y1_test <- y1[n_test]; y2_test <- y2[n_test]; y3_test <- y3[n_test]
  
  # analytical solution
  model1 <- lm(y1_train ~ x1_train + x2_train + x3_train)
  model2 <- lm(y2_train ~ x1_train + x2_train + x3_train)
  model3 <- lm(y3_train ~ x1_train + x2_train + x3_train)
  
  model1_pred <- predict(model1, data.frame(x1_train = x1_test, x2_train = x2_test, x3_train = x3_test), interval = "predict")
  model2_pred <- predict(model2, data.frame(x1_train = x1_test, x2_train = x2_test, x3_train = x3_test), interval = "predict")
  model3_pred <- predict(model3, data.frame(x1_train = x1_test, x2_train = x2_test, x3_train = x3_test), interval = "predict")
  
  model1_cov <- sum((y1_test >= model1_pred[,2]) & (y1_test <= model1_pred[,3])) / length(y1_test)
  model2_cov <- sum((y2_test >= model2_pred[,2]) & (y2_test <= model2_pred[,3])) / length(y2_test)
  model3_cov <- sum((y3_test >= model3_pred[,2]) & (y3_test <= model3_pred[,3])) / length(y3_test)
  
  model1_CI_bw <- mean(model1_pred[,3] - model1_pred[,2])
  model2_CI_bw <- mean(model2_pred[,3] - model2_pred[,2])
  model3_CI_bw <- mean(model3_pred[,3] - model3_pred[,2])
  
  # inductive conformal prediction solution
  x1_train_cp <- x1_train[1:500]; x2_train_cp <- x2_train[1:500]; x3_train_cp <- x3_train[1:500];
  y1_train_cp <- y1_train[1:500]; y2_train_cp <- y2_train[1:500]; y3_train_cp <- y3_train[1:500]; 
  x1_cali_cp <- x1_train[501:1000]; x2_cali_cp <- x2_train[501:1000]; x3_cali_cp <- x3_train[501:1000]; 
  y1_cali_cp <- y1_train[501:1000]; y2_cali_cp <- y2_train[501:1000]; y3_cali_cp <- y3_train[501:1000]; 
  
  model1 <- lm(y1_train_cp ~ x1_train_cp + x2_train_cp + x3_train_cp)
  model2 <- lm(y2_train_cp ~ x1_train_cp + x2_train_cp + x3_train_cp)
  model3 <- lm(y3_train_cp ~ x1_train_cp + x2_train_cp + x3_train_cp)
  
  # calibration set prediction
  model1_pred <- predict(model1, data.frame(x1_train_cp = x1_cali_cp, x2_train_cp = x2_cali_cp, x3_train_cp = x3_cali_cp))
  model2_pred <- predict(model2, data.frame(x1_train_cp = x1_cali_cp, x2_train_cp = x2_cali_cp, x3_train_cp = x3_cali_cp))
  model3_pred <- predict(model3, data.frame(x1_train_cp = x1_cali_cp, x2_train_cp = x2_cali_cp, x3_train_cp = x3_cali_cp))
  d1 <- sort(abs(y1_cali_cp - model1_pred))[ceiling(501 * 0.95)]
  d2 <- sort(abs(y2_cali_cp - model2_pred))[ceiling(501 * 0.95)]
  d3 <- sort(abs(y3_cali_cp - model3_pred))[ceiling(501 * 0.95)]
  
  # test set prediction
  model1_pred <- predict(model1, data.frame(x1_train_cp = x1_test, x2_train_cp = x2_test, x3_train_cp = x3_test))
  model2_pred <- predict(model2, data.frame(x1_train_cp = x1_test, x2_train_cp = x2_test, x3_train_cp = x3_test))
  model3_pred <- predict(model3, data.frame(x1_train_cp = x1_test, x2_train_cp = x2_test, x3_train_cp = x3_test))
  
  model1_cp_cov <- sum((y1_test >= (model1_pred-d1)) & (y1_test <= (model1_pred+d1))) / length(y1_test)
  model2_cp_cov <- sum((y2_test >= (model2_pred-d2)) & (y2_test <= (model2_pred+d2))) / length(y2_test)
  model3_cp_cov <- sum((y3_test >= (model3_pred-d3)) & (y3_test <= (model3_pred+d3))) / length(y3_test)
  
  model1_cp_CI_bw <- 2 * d1
  model2_cp_CI_bw <- 2 * d2
  model3_cp_CI_bw <- 2 * d3
  
  full_result <- rbind(full_result, c(model1_cov, model2_cov, model3_cov, model1_cp_cov, model2_cp_cov, model3_cp_cov,
                                      model1_CI_bw, model2_CI_bw, model3_CI_bw, model1_cp_CI_bw, model2_cp_CI_bw, model3_cp_CI_bw))
}


full_result <- as.data.frame(full_result)
colnames(full_result) <- c("M1_Cov", "M2_Cov", "M3_Cov", "M1_CP_Cov", "M2_CP_Cov", "M3_CP_Cov",
                           "M1_CIwidth", "M2_CIwidth", "M3_CIwidth", "M1_CP_CIwidth", "M2_CP_CIwidth", "M3_CP_CIwidth")

full_result_mean <- full_result %>%
  summarise_all(mean)

full_result_sd <- full_result %>%
  summarise_all(sd)
full_result_sd <- full_result_sd * 1.96

