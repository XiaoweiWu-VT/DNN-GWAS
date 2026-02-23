generatepedstr10 <- function() {
  # for particular simulation only
  # this function generates the size-10 pedigree structure (see Fig 1 in ELNN paper)
  # output: 
  #         a matrix of 10 * 3, 1st column: individual ID, 2nd column: father ID, 3rd column: mother ID
  
  pedstr.mat <- cbind(1 : 10, 0, 0)
  pedstr.mat[3, 2 : 3] <- c(1, 2)
  pedstr.mat[6, 2 : 3] <- c(1, 2)
  pedstr.mat[7, 2 : 3] <- c(3, 4)
  pedstr.mat[8, 2 : 3] <- c(5, 6)
  pedstr.mat[9, 2 : 3] <- c(5, 6)
  pedstr.mat[10, 2 : 3] <- c(5, 6)
  
  return(pedstr.mat)
}

loadPHI.pf10 <- function() {
  # for particular simulation only
  # this function generates kinship matrix for the size-10 pedigree (see Fig 1 in ELNN paper)
  # output: 
  #         a symmetric matrix of 10 * 10
  
  PHI.pf.mat <- rbind(c(   1,    0,  0.5,    0,    0,  0.5,  0.25,  0.25,  0.25,  0.25), 
                      c(   0,    1,  0.5,    0,    0,  0.5,  0.25,  0.25,  0.25,  0.25), 
                      c( 0.5,  0.5,    1,    0,    0,  0.5,   0.5,  0.25,  0.25,  0.25), 
                      c(   0,    0,    0,    1,    0,    0,   0.5,     0,     0,     0), 
                      c(   0,    0,    0,    0,    1,    0,     0,   0.5,   0.5,   0.5), 
                      c( 0.5,  0.5,  0.5,    0,    0,    1,  0.25,   0.5,   0.5,   0.5),
                      c(0.25, 0.25,  0.5,  0.5,    0, 0.25,     1, 0.125, 0.125, 0.125),
                      c(0.25, 0.25, 0.25,    0,  0.5,  0.5, 0.125,     1,   0.5,   0.5),
                      c(0.25, 0.25, 0.25,    0,  0.5,  0.5, 0.125,   0.5,     1,   0.5),
                      c(0.25, 0.25, 0.25,    0,  0.5,  0.5, 0.125,   0.5,   0.5,     1))
  
  return(PHI.pf.mat)
}

generatepedG <- function(pedstr.mat, G.fdr.vec) {
  # this function generates single-SNP genotype data for a predetermined pedigree by gene dropping
  # input:
  #         pedstr.mat: a matrix of size nind * 3
  #         G.fdr.vec: a vector of the genotypes of the founders in the pedigree
  # output: 
  #         a vector of genotypes for all individuals in the pedigree
  
  fdr.ix <- which((pedstr.mat[, 2] == 0) & (pedstr.mat[, 3] == 0))
  if (length(G.fdr.vec) != length(fdr.ix)) stop("Number of founders seems wrong!")
  
  nind <- nrow(pedstr.mat)
  G.vec <- rep(NA, nind)
  G.vec[fdr.ix] <- G.fdr.vec
  
  osp.ix <- setdiff(1 : nind, fdr.ix)
  nosp <- length(osp.ix)
  for (i in 1 : nosp)
  {
    id <- osp.ix[i]
    id.fa <- pedstr.mat[id, 2]
    id.mo <- pedstr.mat[id, 3]
    G.ospi.fa <- G.vec[id.fa]
    G.ospi.mo <- G.vec[id.mo]
    if ((G.ospi.fa == 0) & (G.ospi.mo == 0))
    {
      G.vec[id] <- 0
    }
    if (((G.ospi.fa == 0) & (G.ospi.mo == 1)) | ((G.ospi.fa == 1) & (G.ospi.mo == 0)))
    {
      G.vec[id] <- sample(c(0, 1), 1)
    }
    if (((G.ospi.fa == 0) & (G.ospi.mo == 2)) | ((G.ospi.fa == 2) & (G.ospi.mo == 0)))
    {
      G.vec[id] <- 1
    }
    if ((G.ospi.fa == 1) & (G.ospi.mo == 1))
    {
      G.vec[id] <- sample(c(0, 1, 1, 2), 1)
    }
    if (((G.ospi.fa == 1) & (G.ospi.mo == 2)) | ((G.ospi.fa == 2) & (G.ospi.mo == 1)))
    {
      G.vec[id] <- sample(c(1, 2), 1)
    }
    if ((G.ospi.fa == 2) & (G.ospi.mo == 2))
    {
      G.vec[id] <- 2
    }
  }
  
  return(G.vec)
}

generateG <- function(nfam, pedstr.pf.mat, alfreqforfdr) {
  # for particular simulation only
  # this function generates single-SNP genotype data for multiple families with the same pedigree structure
  # input:
  #         nfam: a scalar, the number of families
  #         pedstr.pf.mat: a matrix of size nind.pf * 3
  #         alfreqforfdr: a scalar, the allele frequency
  # output: 
  #         a vector of genotypes for the individuals in all nfam families
  
  nind.pf <- nrow(pedstr.pf.mat)
  fdr.ic <- (pedstr.pf.mat[, 2] == 0) & (pedstr.pf.mat[, 3] == 0)
  nfdr.pf <- sum(fdr.ic)
  nfdr <- nfdr.pf * nfam
  G.fdr.vec <- rep(NA, nfdr)
  nind <- nind.pf * nfam
  G.vec <- rep(NA, nind)
  for (i in 1 : nfam)
  {
    ix.fdrfi <- ((i - 1) * nfdr.pf + 1) : ((i - 1) * nfdr.pf + nfdr.pf)
    G.fdr.vec[ix.fdrfi] <- rbinom(nfdr.pf, 2, alfreqforfdr)
  }
  for (i in 1 : nfam)
  {
    ix.fi <- ((i - 1) * nind.pf + 1) : ((i - 1) * nind.pf + nind.pf)
    ix.fdrfi <- ((i - 1) * nfdr.pf + 1) : ((i - 1) * nfdr.pf + nfdr.pf)
    G.vec[ix.fi] <- generatepedG(pedstr.pf.mat, G.fdr.vec[ix.fdrfi])
  }
  
  return(G.vec)
}

estimatebeta <- function(s2.vec, Y.vec, W.mat, PHI.lst) {
  # this function estimates covariate effects from the null model
  # Y = W * beta + eps, eps ~ MVN(0, S), S = se2 * I + sa2 * PHI
  # input: 
  #         s2.vec: a vector of length 2, with the components se^2 and sa^2
  #         Y.vec: a vector of length nind, the trait
  #         W.mat: design matrix of size nind by q, including intercept
  #         PHI.lst: a list of length nfam, each containing the kinship matrix of a family
  # output: 
  #         a column vector of length q, the estimated covariate effects
  
  nind <- length(Y.vec)
  if (nrow(W.mat) != nind) stop("Dimension wrong!")
  q <- ncol(W.mat)
  nfam <- length(PHI.lst)
  nind.vec <- sapply(PHI.lst, nrow)
  if (prod(sapply(PHI.lst, ncol) == nind.vec) * (sum(nind.vec) == nind) == 0) stop("Dimension wrong!")
  
  WtSinvW <- matrix(0, q, q)
  WtSinvY <- matrix(0, q, 1)
  ed <- cumsum(nind.vec)
  ed1 <- ed + 1
  st <- c(1, ed1[-nfam])
  for (i in 1 : nfam)
  {
    Y.fi <- cbind(Y.vec[st[i] : ed[i]])
    W.fi <- matrix(W.mat[st[i] : ed[i], ], ncol = q)
    Wt.fi <- t(W.fi)
    nind.fi <- nind.vec[i]
    PHI.fi <- PHI.lst[[i]]
    S.fi <- s2.vec[1] * diag(nind.fi) + s2.vec[2] * PHI.fi
    S.fi.inv <- solve(S.fi)
    WtSinvW <- WtSinvW + Wt.fi %*% S.fi.inv %*% W.fi
    WtSinvY <- WtSinvY + Wt.fi %*% S.fi.inv %*% Y.fi
  }
  beta.hat.vec <- solve(WtSinvW) %*% WtSinvY
  
  return(beta.hat.vec)
}

nllnull <- function(s2.vec, Y.vec, W.mat, PHI.lst) {
  # this function calculates negative log-likelihood for the null model
  # Y = W * beta + eps, eps ~ MVN(0, S), S = se2 * I + sa2 * PHI
  # input: 
  #         s2.vec: a vector of length 2, containing two components: se^2 and sa^2
  #         Y.vec: a vector of length nind, the trait
  #         W.mat: design matrix of size nind by q, including intercept
  #         PHI.lst: a list of length nfam, each containing the kinship matrix of a family
  # output: 
  #         a scalar, the negative log-likelihood
  
  lb <- 1e7
  nfam <- length(PHI.lst)
  nind.vec <- sapply(PHI.lst, nrow)
  beta.hat.vec <- estimatebeta(s2.vec, Y.vec, W.mat, PHI.lst)
  mu.vec <- W.mat %*% beta.hat.vec
  
  ll <- 0
  ed <- cumsum(nind.vec)
  ed1 <- ed + 1
  st <- c(1, ed1[-nfam])
  for (i in 1 : nfam)
  {
    Y.fi <- Y.vec[st[i] : ed[i]]
    mu.fi <- mu.vec[st[i] : ed[i]]
    nind.fi <- nind.vec[i]
    PHI.fi <- PHI.lst[[i]]
    S.fi <- s2.vec[1] * diag(nind.fi) + s2.vec[2] * PHI.fi
    ll <- ll + dmvnorm(Y.fi, mu.fi, S.fi, log = T)
  }
  nll <- ifelse(is.infinite(-ll), lb, -ll)
  
  return(nll)
}

estimates2 <- function(Y.vec, W.mat, PHI.lst) {
  # this function estimates variance components from the null model
  # Y = W * beta + eps, eps ~ MVN(0, S), S = se2 * I + sa2 * PHI
  # input: 
  #         Y.vec: a vector of length nind, the trait
  #         W.mat: design matrix of size nind by q, including intercept
  #         PHI.lst: a list of length nfam, each containing the kinship matrix of a family
  # output: 
  #         a vector of length 2, (se2, sa2)
  
  nind <- length(Y.vec)
  if (nrow(W.mat) != nind) stop("Dimension wrong!")
  q <- ncol(W.mat)
  nfam <- length(PHI.lst)
  nind.vec <- sapply(PHI.lst, nrow)
  if (prod(sapply(PHI.lst, ncol) == nind.vec) * (sum(nind.vec) == nind) == 0) stop("Dimension wrong!")
  
  para0.vec <- c(1, 1)
  res <- optim(par = para0.vec, fn = nllnull, NULL, method = "L-BFGS-B", lower = c(0.01, 0.01), upper = c(10, 10), Y.vec = Y.vec, W.mat = W.mat, PHI.lst = PHI.lst)
  
  return(res$par)
}

calculateR <- function(Y.vec, W.mat, PHI.lst) {
  # this function calculates transformed phenotypic residual under H0
  # input: 
  #         Y.vec: a vector of length nind, the trait
  #         W.mat: design matrix of size nind by q, including intercept
  #         PHI.lst: a list of length nfam, each containing the kinship matrix of a family
  # output: 
  #         a vector of length nind
  nind <- length(Y.vec)
  if (nrow(W.mat) != nind) stop("Dimension wrong!")
  q <- ncol(W.mat)
  nfam <- length(PHI.lst)
  nind.vec <- sapply(PHI.lst, nrow)
  if (prod(sapply(PHI.lst, ncol) == nind.vec) * (sum(nind.vec) == nind) == 0) stop("Dimension wrong!")
  
  s2.hat.vec <- estimates2(Y.vec, W.mat, PHI.lst)
  beta.hat.vec <- estimatebeta(s2.hat.vec, Y.vec, W.mat, PHI.lst)
  
  ed <- cumsum(nind.vec)
  ed1 <- ed + 1
  st <- c(1, ed1[-nfam])
  R.vec <- rep(NA, nind)
  for (i in 1 : nfam)
  {
    Y.fi <- cbind(Y.vec[st[i] : ed[i]])
    W.fi <- W.mat[st[i] : ed[i], ]
    Wt.fi <- t(W.fi)
    nind.fi <- nind.vec[i]
    PHI.fi <- PHI.lst[[i]]
    S.fi <- s2.hat.vec[1] * diag(nind.fi) + s2.hat.vec[2] * PHI.fi
    S.fi.inv <- solve(S.fi)
    R.vec[st[i] : ed[i]] <- S.fi.inv %*% (Y.fi - W.fi %*% beta.hat.vec)
  }
  return(R.vec)
}

calculateZstatLRT <- function(Y.vec, W.mat, X.vec, PHI.lst) {
  # this function calculates p-value for the prospective LRT
  # input: 
  #         Y.vec: a vector of length nind, the trait
  #         W.mat: design matrix of size nind by q, including intercept
  #         X.vec: a vector of length nind, the genotype
  #         PHI.lst: a list of length nfam, each containing the kinship matrix of a family
  # output: 
  #         a scalar
  
  nind <- length(Y.vec)
  if (nrow(W.mat) != nind | length(X.vec) != nind) stop("Dimension wrong!")
  q <- ncol(W.mat)
  nfam <- length(PHI.lst)
  nind.vec <- sapply(PHI.lst, nrow)
  if (prod(sapply(PHI.lst, ncol) == nind.vec) * (sum(nind.vec) == nind) == 0) stop("Dimension wrong!")
  
  s2.hat.vec <- estimates2(Y.vec, W.mat, PHI.lst)
  beta.hat.vec <- estimatebeta(s2.hat.vec, Y.vec, W.mat, PHI.lst)
  
  WtSinvW <- matrix(0, q, q)
  # WtSinvY <- matrix(0, q, 1)
  WtSinvX <- matrix(0, q, 1)
  XtSinvX <- 0
  XtR <- 0
  ed <- cumsum(nind.vec)
  ed1 <- ed + 1
  st <- c(1, ed1[-nfam])
  for (i in 1 : nfam)
  {
    Y.fi <- cbind(Y.vec[st[i] : ed[i]])
    W.fi <- matrix(W.mat[st[i] : ed[i], ], ncol = q)
    Wt.fi <- t(W.fi)
    X.fi <- cbind(X.vec[st[i] : ed[i]])
    Xt.fi <- t(X.fi)
    nind.fi <- nind.vec[i]
    PHI.fi <- PHI.lst[[i]]
    S.fi <- s2.hat.vec[1] * diag(nind.fi) + s2.hat.vec[2] * PHI.fi
    S.fi.inv <- solve(S.fi)
    R.fi <- S.fi.inv %*% (Y.fi - W.fi %*% beta.hat.vec)
    XtR <- XtR + sum(X.fi * R.fi) # numerator
    WtSinvW <- WtSinvW + Wt.fi %*% S.fi.inv %*% W.fi
    # WtSinvY <- WtSinvY + Wt.fi %*% S.fi.inv %*% Y.fi
    WtSinvX <- WtSinvX + Wt.fi %*% S.fi.inv %*% X.fi
    XtSinvX <- XtSinvX + Xt.fi %*% S.fi.inv %*% X.fi
  }
  # num <- XtR ^ 2
  # den <- XtSinvX - t(WtSinvX) %*% solve(WtSinvW) %*% WtSinvX
  num <- XtR
  den <- sqrt(XtSinvX - t(WtSinvX) %*% solve(WtSinvW) %*% WtSinvX)
  
  # return(pchisq(num / den, 1, lower.tail = F))
  return(num / den)
}

estimateAF <- function(pedstr.mat, X.vec) {
  # this function estimates allele frequency, pf (BLUE) and pd
  # input: 
  #         pedstr.mat: a matrix of size nind * 3
  #         X.vec: a vector of length nind, the genotype
  # output: 
  #         a vector of length 2
  
  nind <- nrow(pedstr.mat)
  fdr.ix <- which((pedstr.mat[, 2] == 0) & (pedstr.mat[, 3] == 0))
  osp.ix <- setdiff(1 : nind, fdr.ix)
  
  return(c(mean(X.vec[fdr.ix]) / 2, mean(X.vec[osp.ix]) / 2))
}

calculateZstatMAS <- function(Y.vec, W.mat, X.vec, PHI.lst, pedstr.mat) {
  # this function calculates p-value for the retrospective MASTOR
  # input: 
  #         Y.vec: a vector of length nind, the trait
  #         W.mat: design matrix of size nind by q, including intercept
  #         X.vec: a vector of length nind, the genotype
  #         PHI.lst: a list of length nfam, each containing the kinship matrix of a family
  #         pedstr.mat: a matrix of size nind * 3
  # output: 
  #         a scalar
  
  nind <- length(Y.vec)
  if (nrow(W.mat) != nind | length(X.vec) != nind ) stop("Dimension wrong!")
  nfam <- length(PHI.lst)
  nind.vec <- sapply(PHI.lst, nrow)
  if (prod(sapply(PHI.lst, ncol) == nind.vec) * (sum(nind.vec) == nind) == 0) stop("Dimension wrong!")
  
  s2.hat.vec <- estimates2(Y.vec, W.mat, PHI.lst)
  beta.hat.vec <- estimatebeta(s2.hat.vec, Y.vec, W.mat, PHI.lst)
  
  XtR <- 0
  den <- 0
  ed <- cumsum(nind.vec)
  ed1 <- ed + 1
  st <- c(1, ed1[-nfam])
  for (i in 1 : nfam)
  {
    Y.fi <- cbind(Y.vec[st[i] : ed[i]])
    W.fi <- W.mat[st[i] : ed[i], ]
    Wt.fi <- t(W.fi)
    X.fi <- cbind(X.vec[st[i] : ed[i]])
    Xt.fi <- t(X.fi)
    nind.fi <- nind.vec[i]
    PHI.fi <- PHI.lst[[i]]
    S.fi <- s2.hat.vec[1] * diag(nind.fi) + s2.hat.vec[2] * PHI.fi
    S.fi.inv <- solve(S.fi)
    R.fi <- S.fi.inv %*% (Y.fi - W.fi %*% beta.hat.vec)
    XtR <- XtR + sum(X.fi * R.fi)
    den <- den + t(R.fi) %*% PHI.fi %*% R.fi
  }
  # num <- XtR ^ 2
  num <- XtR
  temp <- estimateAF(pedstr.mat, X.vec)
  p.hat <- temp[1]
  sx2.hat <- 2 * p.hat * (1 - p.hat)
  den <- den * sx2.hat
  den <- sqrt(den)
  
  # return(pchisq(num / den, 1, lower.tail = F))
  return(num / den)
}

# simulatedata <- function(nfam, beta.vec, gamma.vec, p.vec, s2.vec) {
#   # for particular simulation only
#   # this function simulates data from the prospective model
#   # Y = W * beta + X * gamma + eps, W = [1, W1, W2], eps ~ MVN(0, S), S = se2 * I + sa2 * PHI
#   # input: 
#   #         nfam: a scalar, the number of families
#   #         beta.vec: a vector of length q, the covariate effects
#   #         gamma.vec: a vector of length nsnp, the genotype effects
#   #         p.vec: a vector of length nsnp, the allele frequencies
#   #         s2.vec: a vector of length 2, (se2, sa2)
#   # output: 
#   #         a matrix of size nsnp * (4 * nind + 2)
#   
#   nsnp <- length(p.vec)
#   data.lst <- vector("list", nsnp)
#   pedstr.pf.mat <- generatepedstr10()
#   PHI.pf.mat <- loadPHI.pf10()
#   s2e <- s2.vec[1]
#   s2a <- s2.vec[2]
#   nind.pf <- nrow(pedstr.pf.mat)
#   nind <- nfam * nind.pf
#   data.mat <- matrix(NA, nsnp, 4 * nind + 2)
#   W1.vec <- scale(runif(nind, 18, 80))
#   W2.vec <- rbinom(nind, 1, 0.7)
#   W.mat <- cbind(1, W1.vec, W2.vec)
#   X.mat <- matrix(NA, nind, nsnp)
#   Y.mat <- matrix(NA, nind, nsnp)
#   PHI.lst <- vector("list", nfam)
#   pedstr.mat <- matrix(NA, nind, 3)
#   for (i in 1 : nfam) {
#     PHI.lst[[i]] <- PHI.pf.mat
#     ix.fi <- ((i - 1) * nind.pf + 1) : ((i - 1) * nind.pf + nind.pf)
#     pedstr.mat[ix.fi, ] <- pedstr.pf.mat
#   }
#   Y.vec <- rep(NA, nind)
#   mu.vec <- W.mat %*% beta.vec
#   SIGMA.pf.mat <- s2e * diag(nind.pf) + s2a * PHI.pf.mat
#   for (j in 1 : nfam)
#   {
#     ix.fj <- ((j - 1) * nind.pf + 1) : ((j - 1) * nind.pf + nind.pf)
#     mu.fj.vec <- mu.vec[ix.fj]
#     Y.vec[ix.fj] <- rmvnorm(1, mu.fj.vec, SIGMA.pf.mat)
#   }
#   
#   for (i in 1 : nsnp) {
#     alfreqforfdr <- p.vec[i]
#     X.tmp <- generateG(nfam, pedstr.pf.mat, alfreqforfdr)
#     X.vec <- scale(X.tmp)
#     X.mat[, i] <- X.vec
#     # mu.vec <- W.mat %*% beta.vec + X.vec * gamma.vec[i]
#     # SIGMA.pf.mat <- s2e * diag(nind.pf) + s2a * PHI.pf.mat
#     # Y.vec <- rep(NA, nind)
#     # for (j in 1 : nfam)
#     # {
#     #   ix.fj <- ((j - 1) * nind.pf + 1) : ((j - 1) * nind.pf + nind.pf)
#     #   mu.fj.vec <- mu.vec[ix.fj]
#     #   Y.vec[ix.fj] <- rmvnorm(1, mu.fj.vec, SIGMA.pf.mat)
#     # }
#     # Y.vec <- rnorm(nind, W.mat %*% beta.vec + X.vec * gamma.vec[i], s2e)  ###
#     # Y.vec <- rnorm(nind, W.mat %*% beta.vec, s2e)
#     # Y.mat[, i] <- Y.vec
#     pvalLRT <- calculatePvalLRT(Y.vec, W.mat, X.vec, PHI.lst)
#     pvalMAS <- calculatePvalMAS(Y.vec, W.mat, X.vec, PHI.lst, pedstr.mat)
#     # fit <- lm(Y.vec ~ W1.vec + W2.vec + X.vec) ###
#     # pval <- summary(fit)$coefficients[4, "Pr(>|t|)"] ###
#     
#     data.mat[i, ] <- c(W1.vec, W2.vec, X.vec, Y.vec, pvalLRT, pvalMAS)
#     # data.mat[i, ] <- c(W1.vec, W2.vec, X.vec, Y.vec, pval, NA)
#   }
#   colnames(data.mat) <- c(paste0("w1_", 1 : nind), paste0("w2_", 1 : nind), paste0("x_", 1 : nind), paste0("y_", 1 : nind), "pvalLRT", "pvalMAS")
#   
#   return(list(data.mat, W.mat, X.mat, Y.mat))
# }

normalize_features <- function(features) {
  # this function does column-wise normalization to a matrix
  # input: 
  #         features: a matrix, each column is a feature
  # output: 
  #         a matrix of the same size as features
  
  n_features <- ncol(features)
  normalized <- matrix(0, nrow = nrow(features), ncol = n_features)
  
  for (i in 1 : n_features) {
    feature_vec <- features[, i]
    normalized[, i] <- (feature_vec - mean(feature_vec)) / sd(feature_vec)
  }
  
  # Handle any NA/Inf values
  normalized[is.na(normalized)] <- 0
  normalized[is.infinite(normalized)] <- 0
  
  return(normalized)
}

build_tabular_mlp_4n1 <- function(input_shape) { # 400 obsolete
  
  model <- keras_model_sequential() %>%
    # Input layer for the flattened features
    # layer_dense(units = 512, activation = "relu", input_shape = input_shape) %>%
    # layer_batch_normalization() %>%
    # layer_dropout(0.6) %>%
    
    # Process aggregated information
    layer_dense(units = 1024, activation = "relu", input_shape = input_shape) %>%
    layer_batch_normalization() %>%
    layer_dropout(0.5) %>%
    
    layer_dense(units = 512, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.4) %>%
    
    layer_dense(units = 256, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.3) %>%
    
    layer_dense(units = 128, activation = "relu") %>%
    layer_dropout(0.2) %>%
    
    layer_dense(units = 64, activation = "relu") %>%
    layer_dropout(0.1) %>%
    
    # Output
    layer_dense(units = 1, activation = "linear")
  
  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    loss = "mse",
    metrics = c("mae")
  )
  
  return(model)
}

build_tabular_mlp_4n <- function(input_shape) {
  
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
    layer_dropout(0.2) %>%
    
    layer_dense(units = 16, activation = "relu") %>%
    layer_dropout(0.1) %>%
    
    # Output
    layer_dense(units = 1, activation = "linear")
  
  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    loss = "mse",
    metrics = c("mae")
  )
  
  return(model)
}

build_tabular_mlp_3np21 <- function(input_shape) { # 302
  
  model <- keras_model_sequential() %>%
    # Input layer for the flattened features
    layer_dense(units = 256, activation = "relu", input_shape = input_shape) %>%
    layer_batch_normalization() %>%
    layer_dropout(0.5) %>%
    
    # Process aggregated information
    layer_dense(units = 128, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.4) %>%
    
    layer_dense(units = 64, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.3) %>%
    
    layer_dense(units = 32, activation = "relu") %>%
    layer_dropout(0.2) %>%
    
    layer_dense(units = 16, activation = "relu") %>%
    layer_dropout(0.1) %>%
    
    # Output
    layer_dense(units = 1, activation = "linear")
  
  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    loss = "mse",
    metrics = c("mae")
  )
  
  return(model)
}

build_tabular_mlp_3np2 <- function(matrix_shape, scalar_dim) { # 302
  
  # Branch 1: Process matrix features (raw data)
  matrix_input <- layer_input(shape = matrix_shape, name = "matrix_input")
  
  matrix_branch <- matrix_input %>%
    # Process each observation
    layer_dense(units = 64, activation = "relu") %>%
    layer_batch_normalization() %>%
    
    # Aggregate across observations
    layer_global_average_pooling_1d() %>%
    
    # Further processing
    layer_dense(units = 32, activation = "relu") %>%
    layer_dropout(0.2) %>%
    
    layer_dense(units = 16, activation = "relu") %>%
    layer_dropout(0.1)
  
  # Branch 2: Process scalar features (summary stats)
  scalar_input <- layer_input(shape = scalar_dim, name = "scalar_input")
  
  scalar_branch <- scalar_input %>%
    layer_dense(units = 16, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.1) %>%
    
    layer_dense(units = 8, activation = "relu")
  
  # Combine branches
  combined <- layer_concatenate(c(matrix_branch, scalar_branch)) %>%
    layer_dense(units = 16, activation = "relu") %>%
    layer_dropout(0.1) %>%
    
    layer_dense(units = 8, activation = "relu") %>%
    
    layer_dense(units = 1, activation = "linear", name = "output")
  
  # Create model
  model <- keras_model(
    inputs = c(matrix_input, scalar_input),
    outputs = combined
  )

  model %>% compile(
    optimizer = optimizer_adam(0.001),
    loss = "mse",
    metrics = c("mae")
  )
  
  return(model)
}

build_tabular_mlp_1np2 <- function(input_shape) { # 102
  
  model <- keras_model_sequential() %>%
    # Input layer for the flattened features
    layer_dense(units = 64, activation = "relu", input_shape = input_shape) %>%
    layer_batch_normalization() %>%
    layer_dropout(0.3) %>%
    
    # Hidden layer 1
    layer_dense(units = 32, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.2) %>%
    
    # Hidden layer 2
    layer_dense(units = 16, activation = "relu") %>%
    layer_dropout(0.1) %>%
    
    # Output layer
    layer_dense(units = 1, activation = "linear")
  
  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    loss = "mse",
    metrics = c("mae")
  )
  
  return(model)
}

build_tabular_mlp_2n <- function(input_shape) {
  # input_shape: (nind, 2)
  
  model <- keras_model_sequential() %>%
    # Process each (x, residual) observation
    layer_dense(units = 8, activation = "relu", input_shape = input_shape) %>%
    layer_batch_normalization() %>%
    
    # Another layer for observation processing
    layer_dense(units = 4, activation = "relu") %>%
    layer_batch_normalization() %>%
    
    # Global pooling to aggregate across observations
    layer_global_average_pooling_1d() %>%
    
    # Process aggregated information
    layer_dense(units = 32, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.3) %>%
    
    layer_dense(units = 16, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.2) %>%
    
    layer_dense(units = 8, activation = "relu") %>%
    
    # Output layer - predict t-statistic
    layer_dense(units = 1, activation = "linear",
                kernel_initializer = initializer_random_normal(stddev = 0.3))
  
  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    loss = variance_preserving_loss,
    metrics = c("mae")
  )
  
  return(model)
}

build_real_data_mlp <- function(input_shape) {
  # input_shape: (n_obs, n_vars) = (8087, 8)
  
  # Create a model that processes each observation independently
  # then aggregates information across observations
  
  inputs <- layer_input(shape = input_shape)
  
  # Process each observation with shared weights (like a 1D convolution without spatial structure)
  # This is permutation-invariant!
  processed <- inputs %>%
    # First, process each observation (8 features -> 16 features)
    layer_dense(units = 16, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.2) %>%
    
    # Further processing
    layer_dense(units = 8, activation = "relu") %>%
    layer_batch_normalization() %>%
    
    # Aggregate across observations (global pooling makes it permutation-invariant)
    layer_global_average_pooling_1d() %>%
    
    # Process aggregated information
    layer_dense(units = 64, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.3) %>%
    
    layer_dense(units = 32, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(0.2) %>%
    
    layer_dense(units = 16, activation = "relu") %>%
    layer_dropout(0.1) %>%
    
    # Output layer
    layer_dense(units = 1, activation = "linear")
  
  model <- keras_model(inputs = inputs, outputs = processed)
  
  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    # loss = "mean_squared_error",
    loss = variance_preserving_loss,
    metrics = c("mean_absolute_error")
  )
  
  return(model)
}

build_tabular_mlp_8n <- function(input_shape) {
  
  model <- keras_model_sequential() %>%
    # First: Process each observation independently
    layer_reshape(target_shape = input_shape, input_shape = input_shape) %>%
    
    # Observation-wise processing
    layer_dense(units = 4096, activation = "relu") %>% ### need adjustment
    # layer_dense(units = 512, activation = "relu") %>% ### need adjustment
    layer_batch_normalization() %>%
    
    # Aggregation across observations (critical!)
    layer_global_average_pooling_1d() %>%
    
    # Process aggregated information
    layer_dense(units = 8192, activation = "relu") %>% ### need adjustment
    # layer_dense(units = 1024, activation = "relu") %>% ### need adjustment
    layer_batch_normalization() %>%
    layer_dropout(0.5) %>%
    # layer_dropout(0.4) %>%
    
    layer_dense(units = 2048, activation = "relu") %>% ### need adjustment
    # layer_dense(units = 512, activation = "relu") %>% ### need adjustment
    layer_batch_normalization() %>%
    layer_dropout(0.4) %>%
    # layer_dropout(0.3) %>%
    
    layer_dense(units = 512, activation = "relu") %>% ### need adjustment
    # layer_dense(units = 128, activation = "relu") %>% ### need adjustment
    layer_dropout(0.3) %>%
    # layer_dropout(0.2) %>%
    
    layer_dense(units = 128, activation = "relu") %>% ### need adjustment
    # layer_dense(units = 32, activation = "relu") %>% ### need adjustment
    layer_dropout(0.2) %>%
    # layer_dropout(0.1) %>%

    layer_dense(units = 32, activation = "relu") %>% ### need adjustment
    layer_dropout(0.1) %>%
    
    # Output
    layer_dense(units = 1, activation = "linear")

  model %>% compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    loss = "mse",
    metrics = c("mae")
  )
  
  return(model)
}

train_predict <- function(x.train.norm, x.test.norm, y.train.scaled, model) {
  
  cat(paste0("MLP-", model, "...\n"))
  if (model == "4n") {
    input.shape <- c(100, 4)
    mlp_model <- build_tabular_mlp_4n(input.shape)
  }
  if (model == "3np2") {
    matrix_shape <- c(100, 3)
    scalar_dim <- 2
    mlp_model <- build_tabular_mlp_3np2(matrix_shape, scalar_dim)
  }
  if (model == "1np2") {
    input.shape <- ncol(x.train.norm)
    mlp_model <- build_tabular_mlp_1np2(input.shape)
  }
  if (model == "8n") {
    input_shape <- c(8087, 8)  # 8087 observations × 8 variables
    mlp_model <- build_real_data_mlp(input_shape)
    # mlp_model <- build_tabular_mlp_8n(input.shape)
  }
  if (model == "2n") {
    input_shape <- c(8087, 2)  # 8087 observations × 2 variables
    mlp_model <- build_tabular_mlp_2n(input_shape)
    # mlp_model <- build_tabular_mlp_8n(input.shape)
  }
  runt <- system.time({
    # Training parameters
    epochs <- 100
    batch_size <- 16  # Smaller batch size due to large input size
    patience <- 15
    
    # Callbacks
    callbacks <- list(
      callback_early_stopping(
        monitor = "val_loss",
        patience = 20,
        restore_best_weights = TRUE
      ),
      callback_reduce_lr_on_plateau(
        monitor = "val_loss",
        factor = 0.5,
        patience = 10
      )
    )
    # callbacks <- list(
    #   callback_early_stopping(
    #     monitor = "val_loss",
    #     patience = patience,
    #     restore_best_weights = TRUE
    #   ),
    #   callback_reduce_lr_on_plateau(
    #     monitor = "val_loss",
    #     factor = 0.5,
    #     patience = patience %/% 2,
    #     min_lr = 1e-6
    #   ),
    #   callback_model_checkpoint(
    #     filepath = "best_mlp_model.keras",
    #     monitor = "val_loss",
    #     save_best_only = TRUE
    #   )
    # )
    
    # Train the model
    history <- mlp_model %>% fit(
      x.train.norm, y.train.scaled,
      # y = y_train_norm,
      # y = y_train,
      # validation_data = list(x_val_norm, y_val_norm),
      # validation_data = list(x_val_norm, y_val),
      validation_split = 0.2,
      epochs = 100,
      batch_size = 32,
      callbacks = callbacks,
      verbose = 1
      # epochs = epochs,
      # batch_size = batch_size,
      # callbacks = callbacks,
      # verbose = 1
    )

    # history <- mlp_model %>% fit(
    #   x.train.norm, y.train.scaled,
    #   epochs = 200,
    #   batch_size = 16,  # Smaller batch size for complex data
    #   validation_split = 0.2,
    #   verbose = 1,
    #   callbacks = list(
    #     callback_early_stopping(patience = 20, restore_best_weights = T),
    #     callback_reduce_lr_on_plateau(factor = 0.5, patience = 8),
    #     callback_model_checkpoint("best_model.h5", save_best_only = T)
    #   )
    # )
  })
  y.pred.scaled <- mlp_model %>% predict(x.test.norm)
  
  return(list(y.pred.scaled, runt[3]))
}

# Simple variance-preserving loss
variance_preserving_loss <- function(y_true, y_pred) {
  mse_loss <- tf$reduce_mean(tf$square(y_true - y_pred))
  
  # Encourage predictions to have similar variance to true values
  pred_var <- tf$math$reduce_variance(y_pred)
  true_var <- tf$math$reduce_variance(y_true)
  var_loss <- tf$square(pred_var - true_var)
  
  # Weighted combination
  return(mse_loss + 0.02 * var_loss)
}

# variance_preserving_loss <- function(y_true, y_pred) {
#   mse_loss <- tf$reduce_mean(tf$square(y_true - y_pred))
#   
#   # Encourage predictions to have similar variance to true values
#   pred_var <- tf$math$reduce_variance(y_pred)
#   true_var <- tf$math$reduce_variance(y_true)
#   var_loss <- tf$square(pred_var - true_var)
#   
#   # Weighted combination
#   return(mse_loss + 0.02 * var_loss)
# }

get_mode <- function(x) {
  # Create a frequency table
  freq_table <- table(x)
  
  # Sort the table in descending order of frequency
  sorted_table <- sort(freq_table, decreasing = T)
  
  # Extract the name of the most frequent value
  mode_value <- as.numeric(names(sorted_table[1]))
  
  return(mode_value)
}

# imputeg <- function(mat1) {
#   nind <- nrow(mat1)
#   nsnp <- ncol(mat1)
#   mat2 <- mat1
#   for (i in 1 : nsnp) {
#     g.vec <- mat1[, i]
#     # freq_table <- table(g.vec)
#     # sorted_table <- sort(freq_table, decreasing = T)
#     # mode_value <- as.numeric(names(sorted_table[1]))
#     ix <- is.na(g.vec)
#     if (sum(ix) > 0) {
#       g.vec[ix] <- get_mode(g.vec)
#     }
#     mat2[, i] <- g.vec
#   }
#   return(mat2)
# }

impute <- function(mat1) {
  nind <- nrow(mat1)
  nvar <- ncol(mat1)
  mat2 <- mat1
  for (i in 1 : nvar) {
    g.vec <- mat1[, i]
    # freq_table <- table(g.vec)
    # sorted_table <- sort(freq_table, decreasing = T)
    # mode_value <- as.numeric(names(sorted_table[1]))
    ix <- is.na(g.vec)
    if (sum(ix) > 0) {
      g.vec[ix] <- get_mode(g.vec)
    }
    mat2[, i] <- g.vec
  }
  return(mat2)
}