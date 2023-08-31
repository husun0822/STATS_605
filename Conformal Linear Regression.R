library(matlib)
library(datasets)
library(ggplot2)
library(latex2exp)
data(iris)


# ----- predict petal width with sepal length ----- #
x <- iris$Sepal.Length; y <- iris$Petal.Width; n <- length(x)

X <- matrix(1, n, 2); X[,2] <- x; # design matrix with intercept
Xp <- inv(t(X) %*% X) %*% t(X) # pre-compute (XTX)^(-1)XT


# ----- Full Conformal Prediction ----- #
# compute the a and b sequence
y_loo = matrix(y); y_loo[n,1] = 0
e = matrix(0,n,1); e[n,1] = 1               
a <- c(); b <- c();

for (i in 1:(n-1)){
  a_i = y[i] - matrix(c(1, x[i]), 1, 2) %*% Xp %*% y_loo
  b_i = -1 *(matrix(c(1, x[i]), 1, 2) %*% Xp %*% e)
  
  if (b_i < 0){
    a_i = -a_i
    b_i = -b_i
  }
  a <- c(a, a_i)
  b <- c(b, b_i)
}

a_n <- -1 * (matrix(c(1, x[n]), 1, 2) %*% Xp %*% y_loo)
b_n <- 1 - (matrix(c(1, x[n]), 1, 2) %*% Xp %*% e)
if (b_n < 0){
  a_n = -a_n
  b_n = -b_n
}

# ----- conformal prediction ----- # 
y_grid <- seq(from = 0.0, to = 3.0, by = 0.001)
pi_y <- c()

for (y_cand in y_grid){
  s <- 1
  
  for (i in 1:(n-1)){
    if (abs(a[i] + b[i] * y_cand) <= abs(a_n + b_n * y_cand)){
      s <- s+ 1
    }
  }
  pi_y <- c(pi_y, s)
}

cutoff <- range(y_grid[pi_y <= ceiling(n * 0.95)])

plt <- ggplot() + 
  geom_line(aes(x = y_grid, y = pi_y)) +
  geom_vline(aes(xintercept = cutoff), color="red") + 
  xlab("Candidate y") + 
  ylab(unname(TeX(c("$n\\times \\pi(y)$")))) +
  theme_bw() +
  ggtitle("95% Confidence Band for Sample #150 (true y=1.8)") +
  theme(plot.title = element_text(size = 12, face = "bold"))

plt
ggsave("./Single_Sample_Conformal_Function.png", width=15, height=7, units="cm", dpi=600)


# ----- conformal prediction for all 150 sample points ----- #
lb <- c(); ub <- c();
time_start <- Sys.time()
for (p in 1:n){
  # show progress of the CP
  if (p %% 10 == 0){
    print(paste(c("Progress: ", p, "/", n), collapse=""))
  }
  
  # compute the a and b sequence
  y_loo = matrix(y); y_loo[p,1] = 0
  e = matrix(0,n,1); e[p,1] = 1               
  a <- c(); b <- c();
  
  for (i in 1:n){
    if (i == p){
      a_i <- -1 * (matrix(c(1, x[p]), 1, 2) %*% Xp %*% y_loo)
      b_i <- 1 - (matrix(c(1, x[p]), 1, 2) %*% Xp %*% e)
    }else{
      a_i <- y[i] - matrix(c(1, x[i]), 1, 2) %*% Xp %*% y_loo
      b_i <- -1 *(matrix(c(1, x[i]), 1, 2) %*% Xp %*% e)
    }
    
    if (b_i < 0){
      a_i = -a_i
      b_i = -b_i
    }
    a <- c(a, a_i)
    b <- c(b, b_i)
  }
  
  # conformal prediction
  y_grid <- seq(from = -2.0, to = 5.0, by = 0.001)
  pi_y <- c()
  a_p <- a[p]; b_p <- b[p]
  
  for (y_cand in y_grid){
    s <- 0
    for (i in 1:n){
      if (abs(a[i] + b[i] * y_cand) <= abs(a_p + b_p * y_cand)){
        s <- s+ 1
      }
    }
    pi_y <- c(pi_y, s)
  }
  cutoff <- range(y_grid[pi_y <= ceiling(n * 0.95)])
  lb <- c(lb, cutoff[1]); ub <- c(ub, cutoff[2])
}

time.end <- Sys.time()
time_elapsed <- round(time.end - time_start, 2)
print(time_elapsed)

# ----- plot the prediction interval ----- #
plt2 <- ggplot() + 
  geom_point(aes(x = x, y = y)) +
  geom_errorbar(aes(x = x, ymin = lb, ymax = ub), color = "red") +
  theme_bw() + 
  ggtitle("Transductive Conformal Prediction Interval (Coverage 143/150, 95.3%)") +
  xlab("Petal Width") +
  ylab("Sepal Length") +
  theme(plot.title = element_text(size = 10, face = "bold"))

plt2

ggsave("./Transductive_Conformal_Interval.png", plt2, width=15, height=7, units="cm", dpi=600)


plt_save <- gridExtra::grid.arrange(plt, plt2, nrow = 1)
ggsave("./Conformalized_Linear_Model_Example.png", plt_save, width = 20, height = 6, units = "cm", dpi = 600)

# ----- Baseline Gaussian-Error Predictive Interval ----- #
idx = 1:n
lb_gauss <- c(); ub_gauss <- c();
for (i in 1:n){
  x1 <- x[!idx==i]
  y1 <- y[!idx==i]
  model <- lm(y1 ~ x1)
  res <- predict(model, data.frame(x1 = c(x[i])), interval="predict")
  lb_gauss <- c(lb_gauss, res[2])
  ub_gauss <- c(ub_gauss, res[3])
}

# ----- plot the prediction interval based on Gaussianity ----- #
plt3 <- ggplot() + 
  geom_point(aes(x = x, y = y)) +
  geom_errorbar(aes(x = x, ymin = lb_gauss, ymax = ub_gauss), color = "red") +
  theme_bw() + 
  ggtitle("Linear Regression Prediction Interval (Coverage 144/150, 96%)") +
  xlab("Petal Width") +
  ylab("Sepal Length") +
  theme(plot.title = element_text(size = 10, face = "bold"))

plt3
ggsave("./Prediction_Interval_Gaussian.png", plt3, width = 15, height = 7, units = "cm", dpi = 600)


# ----- Split Conformal Prediction ----- #
train_idx = 1:75
calibration_idx = 76:(n-1)
X_train = X[train_idx,]; y_train = matrix(y[train_idx], ncol = 1)
X_cali = X[calibration_idx,]; y_cali = matrix(y[calibration_idx], ncol = 1)

beta_fit = inv(t(X_train) %*% X_train) %*% t(X_train) %*% y_train
R = abs(y_cali - X_cali %*% beta_fit)
k = ceiling((length(R) + 1) * 0.95)
p1 <- ggplot() + 
  geom_histogram(aes(x = R), binwidth = 0.1, col = "black", fill = "blue") +
  geom_vline(aes(xintercept = sort(R)[k]), color = "red") +
  theme_bw() + 
  ggtitle("Non-conformity Measure Distribution (74 Calibration Samples)") +
  theme(plot.title = element_text(size = 10, face = "bold")) +
  xlab("Non-conformity Measure") + 
  ylab("Bin Count")
  
p1
ggsave("./Non-conformity_Measure_Distribution.png", p1, width=12, height = 6, units = "cm", dpi=600)

# ----- Apply Split Conformal to All 150 Samples ----- #
lb <- c(); ub <- c();
for (p in 1:n){
  idx = 1:n
  idx = idx[!idx == p]
  train_idx = idx[1:75]
  calibration_idx = idx[76:(n-1)]
  X_train = X[train_idx,]; y_train = matrix(y[train_idx], ncol = 1)
  X_cali = X[calibration_idx,]; y_cali = matrix(y[calibration_idx], ncol = 1)
  beta_fit = inv(t(X_train) %*% X_train) %*% t(X_train) %*% y_train
  R = abs(y_cali - X_cali %*% beta_fit)
  k = ceiling((length(R) + 1) * 0.95)
  
  l <- matrix(X[p,], 1, 2) %*% beta_fit - sort(R)[k]
  u <- matrix(X[p,], 1, 2) %*% beta_fit + sort(R)[k]
  
  lb <- c(lb, l)
  ub <- c(ub, u)
}

p2 <- ggplot() + 
  geom_point(aes(x = x, y = y)) +
  geom_errorbar(aes(x = x, ymin = lb, ymax = ub), color = "red") +
  theme_bw() + 
  ggtitle("Inductive Conformal Prediction Interval (Coverage 147/150, 98%)") +
  xlab("Petal Width") +
  ylab("Sepal Length") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p2
ggsave("./Inductive_Conformal_Prediction_Interval.png", p2, width=12, height = 6, units = "cm", dpi=600)

plt_save <- gridExtra::grid.arrange(p1, p2, nrow = 1)
ggsave("./Split_Conformalized_Linear_Model_Example.png", plt_save, width = 20, height = 6, units = "cm", dpi = 600)

