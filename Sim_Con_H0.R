# This code checks consistency among LRT/QLS/DNN-LRT/DNN-QLS/DNN-LRT-RAW/DNN-QLS-RAW under H0
# For Fig.4, Fig.S1, Fig.S2, Fig. S3, and Tab.2 (wo training/test division)
# source("C:/Users/admin-xw/Documents/CurrentProject/EL-CIMB/Code/Sim_Con_H0.R")
library(mvtnorm)
library(keras3)
library(tensorflow)
library(abind)

workpath <- "C:/Users/admin-xw/Documents/CurrentProject/EL-CIMB"
source(paste0(workpath, "/Code/funEL.R"))
tensorflow::set_random_seed(123)

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

cat("Running simulations...\n")
t <- 0
for (i in 1 : nsnp) {
  runt <- system.time({
    alfreqforfdr <- p.vec[i]
    X.vec <- generateG(nfam, pedstr.pf.mat, alfreqforfdr)
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
  })
  t <- t + runt[3]
  if (i %% 500 == 0) cat("Completed", i, "SNPs, time =", t, "\n")
}

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

feature.norm <- normalize_design_matrix(feature.ary)
x.norm <- feature.norm[, , c(1, 5)]
x.norm.RAW <- feature.norm[, , 1 : 4]
y.LRT.mean <- mean(response1.vec)
y.LRT.sd <- sd(response1.vec)
y.LRT.scaled <- (response1.vec - y.LRT.mean) / y.LRT.sd
y.QLS.mean <- mean(response2.vec)
y.QLS.sd <- sd(response2.vec)
y.QLS.scaled <- (response2.vec - y.QLS.mean) / y.QLS.sd
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

mlp_model_LRT <- build_tabular_mlp(input.shape)
mlp_model_QLS <- build_tabular_mlp(input.shape)
mlp_model_LRT_RAW <- build_tabular_mlp_RAW(input.shape.RAW)
mlp_model_QLS_RAW <- build_tabular_mlp_RAW(input.shape.RAW)

cat("Training MLP with LRT response...\n")
history <- mlp_model_LRT %>% fit(
  x.norm, y.LRT.scaled,
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

cat("Training MLP with QLS response...\n")
history <- mlp_model_QLS %>% fit(
  x.norm, y.QLS.scaled,
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

cat("Training MLP-RAW with LRT response...\n")
history <- mlp_model_LRT_RAW %>% fit(
  x.norm.RAW, y.LRT.scaled,
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

cat("Training MLP-RAW with QLS response...\n")
history <- mlp_model_QLS_RAW %>% fit(
  x.norm.RAW, y.QLS.scaled,
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

# Predictions
y.LRT.pred.scaled <- mlp_model_LRT %>% predict(x.norm)
y.LRT.pred <- y.LRT.pred.scaled * y.LRT.sd + y.LRT.mean
y.QLS.pred.scaled <- mlp_model_QLS %>% predict(x.norm)
y.QLS.pred <- y.QLS.pred.scaled * y.QLS.sd + y.QLS.mean
y.LRT.RAW.pred.scaled <- mlp_model_LRT_RAW %>% predict(x.norm.RAW)
y.LRT.RAW.pred <- y.LRT.RAW.pred.scaled * y.LRT.sd + y.LRT.mean
y.QLS.RAW.pred.scaled <- mlp_model_QLS_RAW %>% predict(x.norm.RAW)
y.QLS.RAW.pred <- y.QLS.RAW.pred.scaled * y.QLS.sd + y.QLS.mean

# Calculate metrics
mae.LRT <- mean(abs(response1.vec - y.LRT.pred))
rmse.LRT <- sqrt(mean((response1.vec - y.LRT.pred)^2))
r2.LRT <- 1 - sum((response1.vec - y.LRT.pred)^2) / sum((response1.vec - mean(response1.vec))^2)
correlation.LRT <- cor(response1.vec, y.LRT.pred)
cat("\nPerformance Metrics for LRT Approximation:\n")
cat("  MAE:", round(mae.LRT, 4), "\n")
cat("  RMSE:", round(rmse.LRT, 4), "\n")
cat("  R-squared:", round(r2.LRT, 4), "\n")
cat("  Correlation:", round(correlation.LRT, 4), "\n")

mae.QLS <- mean(abs(response2.vec - y.QLS.pred))
rmse.QLS <- sqrt(mean((response2.vec - y.QLS.pred)^2))
r2.QLS <- 1 - sum((response2.vec - y.QLS.pred)^2) / sum((response2.vec - mean(response2.vec))^2)
correlation.QLS <- cor(response2.vec, y.QLS.pred)
cat("\nPerformance Metrics for QLS Approximation:\n")
cat("  MAE:", round(mae.QLS, 4), "\n")
cat("  RMSE:", round(rmse.QLS, 4), "\n")
cat("  R-squared:", round(r2.QLS, 4), "\n")
cat("  Correlation:", round(correlation.QLS, 4), "\n")

mae.LRT.RAW <- mean(abs(response1.vec - y.LRT.RAW.pred))
rmse.LRT.RAW <- sqrt(mean((response1.vec - y.LRT.RAW.pred)^2))
r2.LRT.RAW <- 1 - sum((response1.vec - y.LRT.RAW.pred)^2) / sum((response1.vec - mean(response1.vec))^2)
correlation.LRT.RAW <- cor(response1.vec, y.LRT.RAW.pred)
cat("\nPerformance Metrics for LRT-RAW Approximation:\n")
cat("  MAE:", round(mae.LRT.RAW, 4), "\n")
cat("  RMSE:", round(rmse.LRT.RAW, 4), "\n")
cat("  R-squared:", round(r2.LRT.RAW, 4), "\n")
cat("  Correlation:", round(correlation.LRT.RAW, 4), "\n")

mae.QLS.RAW <- mean(abs(response2.vec - y.QLS.RAW.pred))
rmse.QLS.RAW <- sqrt(mean((response2.vec - y.QLS.RAW.pred)^2))
r2.QLS.RAW <- 1 - sum((response2.vec - y.QLS.RAW.pred)^2) / sum((response2.vec - mean(response2.vec))^2)
correlation.QLS.RAW <- cor(response2.vec, y.QLS.RAW.pred)
cat("\nPerformance Metrics for QLS-RAW Approximation:\n")
cat("  MAE:", round(mae.QLS.RAW, 4), "\n")
cat("  RMSE:", round(rmse.QLS.RAW, 4), "\n")
cat("  R-squared:", round(r2.QLS.RAW, 4), "\n")
cat("  Correlation:", round(correlation.QLS.RAW, 4), "\n")

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

# x11()
# setEPS()
# postscript(paste(workpath, "/Result/Fig4.eps", sep = ""), width = 7, height = 7)
# par(mfrow = c(2, 2))
# plot(response1.vec, y1.pred, type = "p", pch = 20, xlab = "Z_LRT", ylab = "Z_DNNLRT", main = "(A)")
# abline(0, 1, col = "red")
# plot(response2.vec, y2.pred, type = "p", pch = 20, xlab = "Z_QLS", ylab = "Z_DNNQLS", main = "(B)")
# abline(0, 1, col = "red")
# plot(response1.vec, response2.vec, type = "p", pch = 20, xlab = "Z_LRT", ylab = "Z_QLS", main = "(C)")
# abline(0, 1, col = "red")
# plot(y1.pred, y2.pred, type = "p", pch = 20, xlab = "Z_DNNLRT", ylab = "Z_DNNQLS", main = "(D)")
# abline(0, 1, col = "red")
# dev.off()

save.image(file = paste(workpath, "/Result/Fig4.RData", sep = ""))

# x11()
# par(mfrow = c(2, 2))
# plot.ecdf(response1.vec ^ 2, main = "S_LRT")
# xx = seq(0, max(response1.vec ^ 2), by = 0.1)
# lines(xx, pchisq(xx, 1), col = "red")
# plot.ecdf(response2.vec ^ 2, main = "S_QLS")
# xx = seq(0, max(response2.vec ^ 2), by = 0.1)
# lines(xx, pchisq(xx, 1), col = "red")
# plot.ecdf(y1.pred ^ 2, main = "S_DNNLRT")
# xx = seq(0, max(y1.pred ^ 2), by = 0.1)
# lines(xx, pchisq(xx, 1), col = "red")
# plot.ecdf(y2.pred ^ 2, main = "S_DNNQLS")
# xx = seq(0, max(y2.pred ^ 2), by = 0.1)
# lines(xx, pchisq(xx, 1), col = "red")

# x11()
# plot(-log10(pQLS.vec), -log10(pDNNQLS.vec), type = "p", pch = 20)
# abline(0, 1, col = "red")

cat("\n=== Report Results ===\n")
pLRT.vec = pchisq(response1.vec ^ 2, 1, lower.tail = F)
pQLS.vec = pchisq(response2.vec ^ 2, 1, lower.tail = F)
cat("Type-I error for LRT:", mean(pLRT.vec < 0.05), "\n")
cat("Type-I error for QLS:", mean(pQLS.vec < 0.05), "\n")
pDNNLRT.vec = pchisq(y.LRT.pred ^ 2, 1, lower.tail = F)
pDNNQLS.vec = pchisq(y.QLS.pred ^ 2, 1, lower.tail = F)
cat("Type-I error for DNN-LRT:", mean(pDNNLRT.vec < 0.05), "\n")
cat("Type-I error for DNN-QLS:", mean(pDNNQLS.vec < 0.05), "\n")
pDNNLRTRAW.vec = pchisq(y.LRT.RAW.pred ^ 2, 1, lower.tail = F)
pDNNQLSRAW.vec = pchisq(y.QLS.RAW.pred ^ 2, 1, lower.tail = F)
cat("Type-I error for DNN-LRT-RAW:", mean(pDNNLRTRAW.vec < 0.05), "\n")
cat("Type-I error for DNN-QLS-RAW:", mean(pDNNQLSRAW.vec < 0.05), "\n")

# x11()
setEPS()
postscript(paste(workpath, "/Result/Fig4.eps", sep = ""), width = 7, height = 7)
par(mfrow = c(3, 2))
plot(response1.vec, y.LRT.pred, type = "p", pch = 20, xlab = expression(italic("Z")[LS]), ylab = expression(italic("Z")[DNN-LS]), main = "(A)")
abline(0, 1, col = "red")
plot(response2.vec, y.QLS.pred, type = "p", pch = 20, xlab = expression(italic("Z")[QLS]), ylab = expression(italic("Z")[DNN-QLS]), main = "(B)")
abline(0, 1, col = "red")
plot(response1.vec, response2.vec, type = "p", pch = 20, xlab = expression(italic("Z")[LS]), ylab = expression(italic("Z")[QLS]), main = "(C)")
abline(0, 1, col = "red")
plot(y.LRT.pred, y.QLS.pred, type = "p", pch = 20, xlab = expression(italic("Z")[DNN-LS]), ylab = expression(italic("Z")[DNN-QLS]), main = "(D)")
abline(0, 1, col = "red")
plot(response1.vec, y.LRT.RAW.pred, type = "p", pch = 20, xlab = expression(italic("Z")[LS]), ylab = expression(italic("Z")[DNN-LS-RAW]), main = "(E)")
abline(0, 1, col = "red")
plot(response2.vec, y.QLS.RAW.pred, type = "p", pch = 20, xlab = expression(italic("Z")[QLS]), ylab = expression(italic("Z")[DNN-QLS-RAW]), main = "(F)")
abline(0, 1, col = "red")
dev.off()

setEPS()
postscript(paste(workpath, "/Result/FigS1.eps", sep = ""), width = 7, height = 7)
par(mfrow = c(3, 2))
plot(response1.vec ^ 2, y.LRT.pred ^ 2, type = "p", pch = 20, xlab = expression(italic("S")[LS]), ylab = expression(italic("S")[DNN-LS]), main = "(A)")
abline(0, 1, col = "red")
plot(response2.vec ^ 2, y.QLS.pred ^ 2, type = "p", pch = 20, xlab = expression(italic("S")[QLS]), ylab = expression(italic("S")[DNN-QLS]), main = "(B)")
abline(0, 1, col = "red")
plot(response1.vec ^ 2, response2.vec ^ 2, type = "p", pch = 20, xlab = expression(italic("S")[LS]), ylab = expression(italic("S")[QLS]), main = "(C)")
abline(0, 1, col = "red")
plot(y.LRT.pred ^ 2, y.QLS.pred ^ 2, type = "p", pch = 20, xlab = expression(italic("S")[DNN-LS]), ylab = expression(italic("S")[DNN-QLS]), main = "(D)")
abline(0, 1, col = "red")
plot(response1.vec ^ 2, y.LRT.RAW.pred ^ 2, type = "p", pch = 20, xlab = expression(italic("S")[LS]), ylab = expression(italic("S")[DNN-LS-RAW]), main = "(E)")
abline(0, 1, col = "red")
plot(response2.vec ^ 2, y.QLS.RAW.pred ^ 2, type = "p", pch = 20, xlab = expression(italic("S")[QLS]), ylab = expression(italic("S")[DNN-QLS-RAW]), main = "(F)")
abline(0, 1, col = "red")
dev.off()

setEPS()
postscript(paste(workpath, "/Result/FigS2.eps", sep = ""), width = 7, height = 7)
par(mfrow = c(3, 2))
plot(pLRT.vec, pDNNLRT.vec, type = "p", pch = 20, xlab = expression(Pval[LS]), ylab = expression(Pval[DNN-LS]), main = "(A)")
abline(0, 1, col = "red")
plot(pQLS.vec, pDNNQLS.vec, type = "p", pch = 20, xlab = expression(Pval[QLS]), ylab = expression(Pval[DNN-QLS]), main = "(B)")
abline(0, 1, col = "red")
plot(pLRT.vec, pQLS.vec, type = "p", pch = 20, xlab = expression(Pval[LS]), ylab = expression(Pval[QLS]), main = "(C)")
abline(0, 1, col = "red")
plot(pDNNLRT.vec, pDNNQLS.vec, type = "p", pch = 20, xlab = expression(Pval[DNN-LS]), ylab = expression(Pval[DNN-QLS]), main = "(D)")
abline(0, 1, col = "red")
plot(pLRT.vec, pDNNLRTRAW.vec, type = "p", pch = 20, xlab = expression(Pval[LS]), ylab = expression(Pval[DNN-LS-RAW]), main = "(E)")
abline(0, 1, col = "red")
plot(pQLS.vec, pDNNQLSRAW.vec, type = "p", pch = 20, xlab = expression(Pval[QLS]), ylab = expression(Pval[DNN-QLS-RAW]), main = "(F)")
abline(0, 1, col = "red")
dev.off()

setEPS()
postscript(paste(workpath, "/Result/FigS3.eps", sep = ""), width = 7, height = 7)
par(mfrow = c(3, 2))
qqnorm(response1.vec, main = expression(italic("Z")[LS]))
qqline(response1.vec, col = "red")
qqnorm(response2.vec, main = expression(italic("Z")[QLS]))
qqline(response2.vec, col = "red")
qqnorm(y.LRT.pred, main = expression(italic("Z")[DNN-LS]))
qqline(y.LRT.pred, col = "red")
qqnorm(y.QLS.pred, main = expression(italic("Z")[DNN-QLS]))
qqline(y.QLS.pred, col = "red")
qqnorm(y.LRT.RAW.pred, main = expression(italic("Z")[DNN-LS-RAW]))
qqline(y.LRT.RAW.pred, col = "red")
qqnorm(y.QLS.RAW.pred, main = expression(italic("Z")[DNN-QLS-RAW]))
qqline(y.QLS.RAW.pred, col = "red")
dev.off()
