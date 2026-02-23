# This code evaluates LRT/QLS/DNN-LRT/DNN-QLS/DNN-ENS/DNN-LRT-RAW/DNN-QLS-RAW/DNN-LRT-AF/DNN-QLS-AF/DNN-ENS-AF under H0
# For Tab.3, Tab.4, and Tab.5 (with training/test division)
# source("C:/Users/admin-xw/Documents/CurrentProject/EL-CIMB/Code/Sim_Eva_All_H0.R")
library(mvtnorm)
library(keras3)
library(tensorflow)
library(abind)

workpath <- "C:/Users/admin-xw/Documents/CurrentProject/EL-CIMB"
source(paste0(workpath, "/Code/funEL.R"))
tensorflow::set_random_seed(123)

set.seed(123)

nfam <- 10                      # number of families
nind <- nfam * 10
nsnp <- 5000                    # number of snps
beta.vec <- c(1, 0.5, 0.8)      # covariate effects
gamma.vec <- rep(0, nsnp)       # genotype effects
p.vec <- runif(nsnp, 0.1, 0.5)  # allele frequencies
s2.vec <- c(1.5, 2.5)           # variance components

pedstr.pf.mat <- generatepedstr10()
PHI.pf.mat <- loadPHI.pf10()
s2e <- s2.vec[1]
s2a <- s2.vec[2]
nind.pf <- nrow(pedstr.pf.mat)
nind <- nfam * nind.pf

W1.vec <- runif(nind, 18, 80)
W2.vec <- rbinom(nind, 1, 0.7)
W.mat <- cbind(1, W1.vec, W2.vec)
pedstr.mat <- matrix(NA, nind, 3)
PHI.lst <- vector("list", nfam)
for (i in 1 : nfam) {
  PHI.lst[[i]] <- PHI.pf.mat
  ix.fi <- ((i - 1) * nind.pf + 1) : ((i - 1) * nind.pf + nind.pf)
  pedstr.mat[ix.fi, ] <- pedstr.pf.mat
}

response1.vec <- rep(NA, nsnp)
response2.vec <- rep(NA, nsnp)
feature.ary <- array(NA, dim = c(nsnp, nind, 5))  # 3D array
feature.mat <- matrix(NA, nsnp, 3)

cat("Running simulations...\n")
t <- 0
for (i in 1 : nsnp) {
  runt <- system.time({
    alfreqforfdr <- p.vec[i]
    X.vec <- generateG(nfam, pedstr.pf.mat, alfreqforfdr)
    phat.vec <- estimateAF(pedstr.mat, X.vec)
    Y.vec <- rep(NA, nind)
    mu.vec <- W.mat %*% beta.vec + X.vec * gamma.vec[i]
    SIGMA.pf.mat <- s2e * diag(nind.pf) + s2a * PHI.pf.mat
    for (j in 1 : nfam)
    {
      ix.fj <- ((j - 1) * nind.pf + 1) : ((j - 1) * nind.pf + nind.pf)
      mu.fj.vec <- mu.vec[ix.fj]
      Y.vec[ix.fj] <- rmvnorm(1, mu.fj.vec, SIGMA.pf.mat)
    }
    R.vec <- calculateR(Y.vec, W.mat, X.vec, PHI.lst)
    response1.vec[i] <- calculateZstatLRT(Y.vec, W.mat, X.vec, PHI.lst)
    response2.vec[i] <- calculateZstatMAS(Y.vec, W.mat, X.vec, PHI.lst, pedstr.mat)
    feature.ary[i, , ] <- cbind(X.vec, W1.vec, W2.vec, Y.vec, R.vec)
    feature.mat[i, ] <- c(phat.vec, sum(X.vec * R.vec))
  })
  t <- t + runt[3]
  if (i %% 500 == 0) cat("Completed", i, "SNPs, time =", t, "\n")
}
response.vec <- ifelse(abs(response1.vec) > abs(response2.vec), response1.vec, response2.vec)

# NORMALIZE PROPERLY - per feature, not per observation
normalize_design_matrix <- function(matrices) {
  # matrices shape: (n_sim, n_obs, n_vars)
  n_vars <- dim(matrices)[3]
  normalized <- array(0, dim = dim(matrices))
  
  for (var_idx in 1:n_vars) {
    # Normalize each variable across all simulations and observations
    var_data <- matrices[, , var_idx]
    var_mean <- mean(var_data)
    var_sd <- sd(var_data)
    normalized[, , var_idx] <- (var_data - var_mean) / var_sd
  }
  
  return(normalized)
}

normalize_features <- function(features) {
  n_features <- ncol(features)
  normalized <- matrix(0, nrow = nrow(features), ncol = n_features)
  
  for (i in 1 : n_features) {
    feature_vec <- features[, i]
    normalized[, i] <- (feature_vec - mean(feature_vec)) / sd(feature_vec)
  }
  normalized[is.na(normalized)] <- 0
  normalized[is.infinite(normalized)] <- 0
  
  return(normalized)
}

feature.norm <- normalize_design_matrix(feature.ary)
feature.AF.norm <- normalize_features(feature.mat)

# divide data into training and test sets
train.idx <- sample(1 : nsnp, floor(0.8 * nsnp))
test.idx <- setdiff(1 : nsnp, train.idx)
x.norm.train <- feature.norm[train.idx, , c(1, 5)]
x.norm.test <- feature.norm[test.idx, , c(1, 5)]
x.norm.RAW.train <- feature.norm[train.idx, , 1 : 4]
x.norm.RAW.test <- feature.norm[test.idx, , 1 : 4]
x.norm.AF.train <- feature.AF.norm[train.idx, ]
x.norm.AF.test <- feature.AF.norm[test.idx, ]

y.LRT.train <- response1.vec[train.idx]
y.LRT.train.mean <- mean(y.LRT.train)
y.LRT.train.sd <- sd(y.LRT.train)
y.LRT.train.scaled <- (y.LRT.train - y.LRT.train.mean) / y.LRT.train.sd
y.LRT.test <- response1.vec[test.idx]
y.LRT.test.mean <- mean(y.LRT.test)
y.LRT.test.sd <- sd(y.LRT.test)
y.LRT.test.scaled <- (y.LRT.test - y.LRT.test.mean) / y.LRT.test.sd
y.QLS.train <- response2.vec[train.idx]
y.QLS.train.mean <- mean(y.QLS.train)
y.QLS.train.sd <- sd(y.QLS.train)
y.QLS.train.scaled <- (y.QLS.train - y.QLS.train.mean) / y.QLS.train.sd
y.QLS.test <- response2.vec[test.idx]
y.QLS.test.mean <- mean(y.QLS.test)
y.QLS.test.sd <- sd(y.QLS.test)
y.QLS.test.scaled <- (y.QLS.test - y.QLS.test.mean) / y.QLS.test.sd
y.ENS.train <- response.vec[train.idx]
y.ENS.train.mean <- mean(y.ENS.train)
y.ENS.train.sd <- sd(y.ENS.train)
y.ENS.train.scaled <- (y.ENS.train - y.ENS.train.mean) / y.ENS.train.sd
y.ENS.test <- response.vec[test.idx]
y.ENS.test.mean <- mean(y.ENS.test)
y.ENS.test.sd <- sd(y.ENS.test)
y.ENS.test.scaled <- (y.ENS.test - y.ENS.test.mean) / y.ENS.test.sd
# may consider unscaled version

variance_preserving_loss <- function(y_true, y_pred) {
  mse_loss <- tf$reduce_mean(tf$square(y_true - y_pred))
  
  # Encourage predictions to have similar variance to true values
  pred_var <- tf$math$reduce_variance(y_pred)
  true_var <- tf$math$reduce_variance(y_true)
  var_loss <- tf$square(pred_var - true_var)
  
  # Weighted combination
  return(mse_loss + 0.02 * var_loss)
}

build_tabular_mlp_RAW <- function(input_shape) {
  
  model <- keras_model_sequential() %>%
    # First: Process each observation independently
    layer_reshape(target_shape = input_shape, input_shape = input_shape) %>%
    
    # Observation-wise processing
    layer_dense(units = 64, activation = "relu") %>%
    layer_batch_normalization() %>%
    
    # Aggregation across observations (critical!)
    layer_global_average_pooling_1d() %>%
    
    # Process aggregated information
    layer_dense(units = 128, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.4) %>%
    
    layer_dense(units = 64, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.3) %>%
    
    layer_dense(units = 32, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.2) %>%
    
    layer_dense(units = 16, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.1) %>%
    
    # Output
    layer_dense(units = 1, activation = "linear")
  
  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    # loss = "mse",
    loss = variance_preserving_loss,
    metrics = c("mae")
  )
  
  return(model)
}

build_tabular_mlp_AF <- function(input_shape) {
  
  model <- keras_model_sequential() %>%
    # Input layers for flattened features
    # layer_dense(units = 64, activation = "relu", input_shape = input_shape) %>%
    # layer_batch_normalization() %>%
    # layer_dropout(0.3) %>%
    # 
    # layer_dense(units = 32, activation = "relu") %>%
    # layer_batch_normalization() %>%
    # layer_dropout(0.2) %>%
    # 
    # layer_dense(units = 16, activation = "relu") %>%
    layer_dense(units = 8, activation = "relu", input_shape = input_shape) %>%
    layer_batch_normalization() %>%
    layer_dropout(0.1) %>%
    
    # Output
    layer_dense(units = 1, activation = "linear")
  
  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    # loss = "mse",
    loss = variance_preserving_loss,
    metrics = c("mae")
  )
  
  return(model)
}

build_tabular_mlp <- function(input_shape) {
  
  model <- keras_model_sequential() %>%
    # First: Process each observation independently
    layer_reshape(target_shape = input_shape, input_shape = input_shape) %>%
    
    # Observation-wise processing
    layer_dense(units = 64, activation = "relu") %>%
    layer_batch_normalization() %>%
    
    # Aggregation across observations (critical!)
    layer_global_average_pooling_1d() %>%
    
    # Process aggregated information
    # layer_dense(units = 128, activation = "relu") %>%
    # layer_batch_normalization() %>%
    # layer_dropout(0.4) %>%
    
    layer_dense(units = 64, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.3) %>%
    
    layer_dense(units = 32, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.2) %>%
    
    layer_dense(units = 16, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.1) %>%
    
    # Output
    layer_dense(units = 1, activation = "linear")
  
  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    # loss = "mse",
    loss = variance_preserving_loss,
    metrics = c("mae")
  )
  
  return(model)
}

input.shape <- c(nind, 2)  # 100 observations × 2 variables
input.shape.RAW <- c(nind, 4)  # 100 observations × 4 variables
input.shape.AF <- 3

mlp_model_LRT <- build_tabular_mlp(input.shape)
mlp_model_QLS <- build_tabular_mlp(input.shape)
mlp_model_ENS <- build_tabular_mlp(input.shape)
mlp_model_LRT_RAW <- build_tabular_mlp_RAW(input.shape.RAW)
mlp_model_QLS_RAW <- build_tabular_mlp_RAW(input.shape.RAW)
mlp_model_LRT_AF <- build_tabular_mlp_AF(input.shape.AF)
mlp_model_QLS_AF <- build_tabular_mlp_AF(input.shape.AF)
mlp_model_ENS_AF <- build_tabular_mlp_AF(input.shape.AF)

cat("Training MLP with LRT response...\n")
runt.LRT <- system.time({
history <- mlp_model_LRT %>% fit(
  x.norm.train, y.LRT.train.scaled,
  epochs = 200,
  batch_size = 16,  # Smaller batch size for complex data
  validation_split = 0.2,
  verbose = 1,
  callbacks = list(
    callback_early_stopping(patience = 20, restore_best_weights = TRUE),
    callback_reduce_lr_on_plateau(factor = 0.5, patience = 8),
    callback_model_checkpoint("best_model.h5", save_best_only = TRUE)
  )
)
})

cat("Training MLP with QLS response...\n")
runt.QLS <- system.time({
history <- mlp_model_QLS %>% fit(
  x.norm.train, y.QLS.train.scaled,
  epochs = 200,
  batch_size = 16,  # Smaller batch size for complex data
  validation_split = 0.2,
  verbose = 1,
  callbacks = list(
    callback_early_stopping(patience = 20, restore_best_weights = TRUE),
    callback_reduce_lr_on_plateau(factor = 0.5, patience = 8),
    callback_model_checkpoint("best_model.h5", save_best_only = TRUE)
  )
)
})

cat("Training MLP with ENS response...\n")
runt.ENS <- system.time({
  history <- mlp_model_ENS %>% fit(
    x.norm.train, y.ENS.train.scaled,
    epochs = 200,
    batch_size = 16,  # Smaller batch size for complex data
    validation_split = 0.2,
    verbose = 1,
    callbacks = list(
      callback_early_stopping(patience = 20, restore_best_weights = TRUE),
      callback_reduce_lr_on_plateau(factor = 0.5, patience = 8),
      callback_model_checkpoint("best_model.h5", save_best_only = TRUE)
    )
  )
})

cat("Training MLP-RAW with LRT response...\n")
runt.LRT.RAW <- system.time({
history <- mlp_model_LRT_RAW %>% fit(
  x.norm.RAW.train, y.LRT.train.scaled,
  epochs = 200,
  batch_size = 16,  # Smaller batch size for complex data
  validation_split = 0.2,
  verbose = 1,
  callbacks = list(
    callback_early_stopping(patience = 20, restore_best_weights = TRUE),
    callback_reduce_lr_on_plateau(factor = 0.5, patience = 8),
    callback_model_checkpoint("best_model.h5", save_best_only = TRUE)
  )
)
})

cat("Training MLP-RAW with QLS response...\n")
runt.QLS.RAW <- system.time({
history <- mlp_model_QLS_RAW %>% fit(
  x.norm.RAW.train, y.QLS.train.scaled,
  epochs = 200,
  batch_size = 16,  # Smaller batch size for complex data
  validation_split = 0.2,
  verbose = 1,
  callbacks = list(
    callback_early_stopping(patience = 20, restore_best_weights = TRUE),
    callback_reduce_lr_on_plateau(factor = 0.5, patience = 8),
    callback_model_checkpoint("best_model.h5", save_best_only = TRUE)
  )
)
})

cat("Training MLP-AF with LRT response...\n")
runt.LRT.AF <- system.time({
  history <- mlp_model_LRT_AF %>% fit(
    x.norm.AF.train, y.LRT.train.scaled,
    epochs = 200,
    batch_size = 16,  # Smaller batch size for complex data
    validation_split = 0.2,
    verbose = 1,
    callbacks = list(
      callback_early_stopping(patience = 20, restore_best_weights = TRUE),
      callback_reduce_lr_on_plateau(factor = 0.5, patience = 8),
      callback_model_checkpoint("best_model.h5", save_best_only = TRUE)
    )
  )
})

cat("Training MLP-AF with QLS response...\n")
runt.QLS.AF <- system.time({
  history <- mlp_model_QLS_AF %>% fit(
    x.norm.AF.train, y.QLS.train.scaled,
    epochs = 200,
    batch_size = 16,  # Smaller batch size for complex data
    validation_split = 0.2,
    verbose = 1,
    callbacks = list(
      callback_early_stopping(patience = 20, restore_best_weights = TRUE),
      callback_reduce_lr_on_plateau(factor = 0.5, patience = 8),
      callback_model_checkpoint("best_model.h5", save_best_only = TRUE)
    )
  )
})

cat("Training MLP-AF with ENS response...\n")
runt.ENS.AF <- system.time({
  history <- mlp_model_ENS_AF %>% fit(
    x.norm.AF.train, y.ENS.train.scaled,
    epochs = 200,
    batch_size = 16,  # Smaller batch size for complex data
    validation_split = 0.2,
    verbose = 1,
    callbacks = list(
      callback_early_stopping(patience = 20, restore_best_weights = TRUE),
      callback_reduce_lr_on_plateau(factor = 0.5, patience = 8),
      callback_model_checkpoint("best_model.h5", save_best_only = TRUE)
    )
  )
})

# Predictions
y.LRT.pred.scaled <- mlp_model_LRT %>% predict(x.norm.test)
y.LRT.pred <- y.LRT.pred.scaled * y.LRT.test.sd + y.LRT.test.mean
y.QLS.pred.scaled <- mlp_model_QLS %>% predict(x.norm.test)
y.QLS.pred <- y.QLS.pred.scaled * y.QLS.test.sd + y.QLS.test.mean
y.ENS.pred.scaled <- mlp_model_ENS %>% predict(x.norm.test)
y.ENS.pred <- y.ENS.pred.scaled * y.ENS.test.sd + y.ENS.test.mean
y.LRT.RAW.pred.scaled <- mlp_model_LRT_RAW %>% predict(x.norm.RAW.test)
y.LRT.RAW.pred <- y.LRT.RAW.pred.scaled * y.LRT.test.sd + y.LRT.test.mean
y.QLS.RAW.pred.scaled <- mlp_model_QLS_RAW %>% predict(x.norm.RAW.test)
y.QLS.RAW.pred <- y.QLS.RAW.pred.scaled * y.QLS.test.sd + y.QLS.test.mean
y.LRT.AF.pred.scaled <- mlp_model_LRT_AF %>% predict(x.norm.AF.test)
y.LRT.AF.pred <- y.LRT.AF.pred.scaled * y.LRT.test.sd + y.LRT.test.mean
y.QLS.AF.pred.scaled <- mlp_model_QLS_AF %>% predict(x.norm.AF.test)
y.QLS.AF.pred <- y.QLS.AF.pred.scaled * y.QLS.test.sd + y.QLS.test.mean
y.ENS.AF.pred.scaled <- mlp_model_ENS_AF %>% predict(x.norm.AF.test)
y.ENS.AF.pred <- y.ENS.AF.pred.scaled * y.ENS.test.sd + y.ENS.test.mean

# Calculate metrics
mae.LRT <- mean(abs(y.LRT.test - y.LRT.pred))
rmse.LRT <- sqrt(mean((y.LRT.test - y.LRT.pred)^2))
r2.LRT <- 1 - sum((y.LRT.test - y.LRT.pred)^2) / sum((y.LRT.test - mean(y.LRT.test))^2)
correlation.LRT <- cor(y.LRT.test, y.LRT.pred)
cat("\nPerformance Metrics for LRT Approximation:\n")
cat("  MAE:", round(mae.LRT, 4), "\n")
cat("  RMSE:", round(rmse.LRT, 4), "\n")
cat("  R-squared:", round(r2.LRT, 4), "\n")
cat("  Correlation:", round(correlation.LRT, 4), "\n")

mae.QLS <- mean(abs(y.QLS.test - y.QLS.pred))
rmse.QLS <- sqrt(mean((y.QLS.test - y.QLS.pred)^2))
r2.QLS <- 1 - sum((y.QLS.test - y.QLS.pred)^2) / sum((y.QLS.test - mean(y.QLS.test))^2)
correlation.QLS <- cor(y.QLS.test, y.QLS.pred)
cat("\nPerformance Metrics for QLS Approximation:\n")
cat("  MAE:", round(mae.QLS, 4), "\n")
cat("  RMSE:", round(rmse.QLS, 4), "\n")
cat("  R-squared:", round(r2.QLS, 4), "\n")
cat("  Correlation:", round(correlation.QLS, 4), "\n")

mae.ENS <- mean(abs(y.ENS.test - y.ENS.pred))
rmse.ENS <- sqrt(mean((y.ENS.test - y.ENS.pred)^2))
r2.ENS <- 1 - sum((y.ENS.test - y.ENS.pred)^2) / sum((y.ENS.test - mean(y.ENS.test))^2)
correlation.ENS <- cor(y.ENS.test, y.ENS.pred)
cat("\nPerformance Metrics for ENS Approximation:\n")
cat("  MAE:", round(mae.ENS, 4), "\n")
cat("  RMSE:", round(rmse.ENS, 4), "\n")
cat("  R-squared:", round(r2.ENS, 4), "\n")
cat("  Correlation:", round(correlation.ENS, 4), "\n")

mae.LRT.RAW <- mean(abs(y.LRT.test - y.LRT.RAW.pred))
rmse.LRT.RAW <- sqrt(mean((y.LRT.test - y.LRT.RAW.pred)^2))
r2.LRT.RAW <- 1 - sum((y.LRT.test - y.LRT.RAW.pred)^2) / sum((y.LRT.test - mean(y.LRT.test))^2)
correlation.LRT.RAW <- cor(y.LRT.test, y.LRT.RAW.pred)
cat("\nPerformance Metrics for LRT-RAW Approximation:\n")
cat("  MAE:", round(mae.LRT.RAW, 4), "\n")
cat("  RMSE:", round(rmse.LRT.RAW, 4), "\n")
cat("  R-squared:", round(r2.LRT.RAW, 4), "\n")
cat("  Correlation:", round(correlation.LRT.RAW, 4), "\n")

mae.QLS.RAW <- mean(abs(y.QLS.test - y.QLS.RAW.pred))
rmse.QLS.RAW <- sqrt(mean((y.QLS.test - y.QLS.RAW.pred)^2))
r2.QLS.RAW <- 1 - sum((y.QLS.test - y.QLS.RAW.pred)^2) / sum((y.QLS.test - mean(y.QLS.test))^2)
correlation.QLS.RAW <- cor(y.QLS.test, y.QLS.RAW.pred)
cat("\nPerformance Metrics for QLS-RAW Approximation:\n")
cat("  MAE:", round(mae.QLS.RAW, 4), "\n")
cat("  RMSE:", round(rmse.QLS.RAW, 4), "\n")
cat("  R-squared:", round(r2.QLS.RAW, 4), "\n")
cat("  Correlation:", round(correlation.QLS.RAW, 4), "\n")

mae.LRT.AF <- mean(abs(y.LRT.test - y.LRT.AF.pred))
rmse.LRT.AF <- sqrt(mean((y.LRT.test - y.LRT.AF.pred)^2))
r2.LRT.AF <- 1 - sum((y.LRT.test - y.LRT.AF.pred)^2) / sum((y.LRT.test - mean(y.LRT.test))^2)
correlation.LRT.AF <- cor(y.LRT.test, y.LRT.AF.pred)
cat("\nPerformance Metrics for LRT-AF Approximation:\n")
cat("  MAE:", round(mae.LRT.AF, 4), "\n")
cat("  RMSE:", round(rmse.LRT.AF, 4), "\n")
cat("  R-squared:", round(r2.LRT.AF, 4), "\n")
cat("  Correlation:", round(correlation.LRT.AF, 4), "\n")

mae.QLS.AF <- mean(abs(y.QLS.test - y.QLS.AF.pred))
rmse.QLS.AF <- sqrt(mean((y.QLS.test - y.QLS.AF.pred)^2))
r2.QLS.AF <- 1 - sum((y.QLS.test - y.QLS.AF.pred)^2) / sum((y.QLS.test - mean(y.QLS.test))^2)
correlation.QLS.AF <- cor(y.QLS.test, y.QLS.AF.pred)
cat("\nPerformance Metrics for QLS-AF Approximation:\n")
cat("  MAE:", round(mae.QLS.AF, 4), "\n")
cat("  RMSE:", round(rmse.QLS.AF, 4), "\n")
cat("  R-squared:", round(r2.QLS.AF, 4), "\n")
cat("  Correlation:", round(correlation.QLS.AF, 4), "\n")

mae.ENS.AF <- mean(abs(y.ENS.test - y.ENS.AF.pred))
rmse.ENS.AF <- sqrt(mean((y.ENS.test - y.ENS.AF.pred)^2))
r2.ENS.AF <- 1 - sum((y.ENS.test - y.ENS.AF.pred)^2) / sum((y.ENS.test - mean(y.ENS.test))^2)
correlation.ENS.AF <- cor(y.ENS.test, y.ENS.AF.pred)
cat("\nPerformance Metrics for ENS-AF Approximation:\n")
cat("  MAE:", round(mae.ENS.AF, 4), "\n")
cat("  RMSE:", round(rmse.ENS.AF, 4), "\n")
cat("  R-squared:", round(r2.ENS.AF, 4), "\n")
cat("  Correlation:", round(correlation.ENS.AF, 4), "\n")

# Check variance
# cat("\nVariance Analysis:\n")
# cat("  True LRT Zstats: mean =", round(mean(response1.vec), 4), 
#     "SD =", round(sd(response1.vec), 4), 
#     "range =", round(range(response1.vec), 3), "\n")
# cat("  Pred LRT Zstats: mean =", round(mean(y1.pred), 4),
#     "SD =", round(sd(y1.pred), 4),
#     "range =", round(range(y1.pred), 3), "\n")
# cat("  Variance ratio (pred/true):", round(var(y1.pred) / var(response1.vec), 4), "\n")
# 
# cat("  True QLS Zstats: mean =", round(mean(response2.vec), 4), 
#     "SD =", round(sd(response2.vec), 4), 
#     "range =", round(range(response2.vec), 3), "\n")
# cat("  Pred QLS Zstats: mean =", round(mean(y2.pred), 4),
#     "SD =", round(sd(y2.pred), 4),
#     "range =", round(range(y2.pred), 3), "\n")
# cat("  Variance ratio (pred/true):", round(var(y2.pred) / var(response2.vec), 4), "\n")

# ============================================
# 11. POST-TRAINING CALIBRATION (IF NEEDED)
# ============================================
# cat("\n=== Post-Training Calibration ===\n")
# 
# if (sd(y1.pred) < 0.8 * sd(response1.vec)) {
#   cat("Applying variance calibration...\n")
#   
#   # Scale predictions to match true variance
#   scaling_factor <- sd(response1.vec) / sd(y1.pred)
#   y_pred_cal <- y1.pred * scaling_factor
#   
#   # Adjust mean
#   pred_mean <- mean(y_pred_cal)
#   true_mean <- mean(response1.vec)
#   y_pred_cal <- y_pred_cal + (true_mean - pred_mean)
#   
#   # Update
#   y_pred_original <- y1.pred
#   y_pred <- y_pred_cal
#   
#   cat("  Scaling factor:", round(scaling_factor, 4), "\n")
#   cat("  Original pred SD:", round(sd(y_pred_original), 4), "\n")
#   cat("  Calibrated pred SD:", round(sd(y_pred), 4), "\n")
#   
#   # Recalculate metrics
#   r2_cal <- 1 - sum((response1.vec - y_pred)^2) / sum((response1.vec - mean(response1.vec))^2)
#   cat("  Calibrated R-squared:", round(r2_cal, 4), "\n")
# } else {
#   cat("Variance looks good - no calibration needed.\n")
# }
# 
# if (sd(y2.pred) < 0.8 * sd(response2.vec)) {
#   cat("Applying variance calibration...\n")
#   
#   # Scale predictions to match true variance
#   scaling_factor <- sd(response2.vec) / sd(y2.pred)
#   y_pred_cal <- y2.pred * scaling_factor
#   
#   # Adjust mean
#   pred_mean <- mean(y_pred_cal)
#   true_mean <- mean(response2.vec)
#   y_pred_cal <- y_pred_cal + (true_mean - pred_mean)
#   
#   # Update
#   y_pred_original <- y2.pred
#   y_pred <- y_pred_cal
#   
#   cat("  Scaling factor:", round(scaling_factor, 4), "\n")
#   cat("  Original pred SD:", round(sd(y_pred_original), 4), "\n")
#   cat("  Calibrated pred SD:", round(sd(y_pred), 4), "\n")
#   
#   # Recalculate metrics
#   r2_cal <- 1 - sum((response2.vec - y_pred)^2) / sum((response2.vec - mean(response2.vec))^2)
#   cat("  Calibrated R-squared:", round(r2_cal, 4), "\n")
# } else {
#   cat("Variance looks good - no calibration needed.\n")
# }

save.image(file = paste(workpath, "/Result/Tab345_H0.RData", sep = ""))

cat("\n=== Report Results ===\n")
pLRT.vec = pchisq(response1.vec ^ 2, 1, lower.tail = F)
pQLS.vec = pchisq(response2.vec ^ 2, 1, lower.tail = F)
pENS.vec = pchisq(response.vec ^ 2, 1, lower.tail = F)
cat("Type-I error for LRT:", mean(pLRT.vec < 0.05), "\n")
cat("Type-I error for QLS:", mean(pQLS.vec < 0.05), "\n")
cat("Type-I error for ENS:", mean(pENS.vec < 0.05), "\n")
pDNNLRT.vec = pchisq(y.LRT.pred ^ 2, 1, lower.tail = F)
pDNNQLS.vec = pchisq(y.QLS.pred ^ 2, 1, lower.tail = F)
pDNNENS.vec = pchisq(y.ENS.pred ^ 2, 1, lower.tail = F)
cat("Type-I error for DNN-LRT:", mean(pDNNLRT.vec < 0.05), "\n")
cat("Type-I error for DNN-QLS:", mean(pDNNQLS.vec < 0.05), "\n")
cat("Type-I error for DNN-ENS:", mean(pDNNENS.vec < 0.05), "\n")
pDNNLRTRAW.vec = pchisq(y.LRT.RAW.pred ^ 2, 1, lower.tail = F)
pDNNQLSRAW.vec = pchisq(y.QLS.RAW.pred ^ 2, 1, lower.tail = F)
cat("Type-I error for DNN-LRT-RAW:", mean(pDNNLRTRAW.vec < 0.05), "\n")
cat("Type-I error for DNN-QLS-RAW:", mean(pDNNQLSRAW.vec < 0.05), "\n")
pDNNLRTAF.vec = pchisq(y.LRT.AF.pred ^ 2, 1, lower.tail = F)
pDNNQLSAF.vec = pchisq(y.QLS.AF.pred ^ 2, 1, lower.tail = F)
pDNNENSAF.vec = pchisq(y.ENS.AF.pred ^ 2, 1, lower.tail = F)
cat("Type-I error for DNN-LRT-AF:", mean(pDNNLRTAF.vec < 0.05), "\n")
cat("Type-I error for DNN-QLS-AF:", mean(pDNNQLSAF.vec < 0.05), "\n")
cat("Type-I error for DNN-ENS-AF:", mean(pDNNENSAF.vec < 0.05), "\n")

cat("Training time for DNN-LRT:", runt.LRT[3], "\n")
cat("Training time for DNN-QLS:", runt.QLS[3], "\n")
cat("Training time for DNN-ENS:", runt.ENS[3], "\n")
cat("Training time for DNN-LRT-RAW:", runt.LRT.RAW[3], "\n")
cat("Training time for DNN-QLS-RAW:", runt.QLS.RAW[3], "\n")
cat("Training time for DNN-LRT-AF:", runt.LRT.AF[3], "\n")
cat("Training time for DNN-QLS-AF:", runt.QLS.AF[3], "\n")
cat("Training time for DNN-ENS-AF:", runt.ENS.AF[3], "\n")
