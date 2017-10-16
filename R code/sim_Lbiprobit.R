##### Lbiprobit simulations #####
M <- 5000
K <- 10
A.perc <- 0.1
r <- 8
beta0 <- c(-1, -2)
alpha <- c(0.2, 0.4)
rho <- 0.6

data <- generate_data(M, K, A.perc, r, beta0, alpha, rho)
result <- Lbiprobit(data$Pvalue, data$A)

# no X
M <- 5000
K <- 0
A.perc <- 0.1
r <- 8
beta0 <- c(-1, -2)
alpha <- c(0.2, 0.4)
rho <- 0.6

data_noX <- generate_data(M, K, A.perc, r, beta0, alpha, rho)
result_noX <- Lbiprobit.noX(data_noX$Pvalue)

library(fCopulae)
generate_data <- function(M, K, A.perc, r, beta0, alpha, rho){
  
  A         <- rep(0, M*K)
  indexA    <- sample(M*K, M*K*A.perc)
  A[indexA] <- 1
  A         <- matrix(A, M, K)
  
  beta    <- matrix(rnorm(2*K), 2, K)
  sigmae2 <- var(A %*% t(beta))/r # r is signal-noise ratio
  beta    <- beta/sqrt(diag(sigmae2))
  beta    <- cbind(as.matrix(beta0), beta)
  
  Z <- cbind(rep(1, M), A) %*% t(beta) + rnorm2d(M, rho)
  
  indexeta <- (Z > 0)
  eta      <- matrix(as.numeric(indexeta), M, 2)
  
  Pvalue <- matrix(runif(2*M), M, 2)
  Pvalue[indexeta[, 1], 1] <- rbeta(sum(indexeta[, 1]), alpha[1], 1)
  Pvalue[indexeta[, 2], 2] <- rbeta(sum(indexeta[, 2]), alpha[2], 1)
  
  return( list(Pvalue = Pvalue, A = A, beta = beta, eta = eta))
}
