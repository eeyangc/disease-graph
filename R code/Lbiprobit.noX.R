# Latent Bivariate Probit Model with no X
library(pbivnorm)

Lbiprobit.noX <- function(Pvalue){
  
  # Stage 1
  cat("Stage 1 starts... \n", " First GWAS ")
  result1.1 <- Update_alpha(Pvalue[, 1])
  cat("iterations=", result1.1$iter.times, "\n")
  cat(" Second GWAS ")
  result1.2 <- Update_alpha(Pvalue[, 2])
  cat("iterations=", result1.2$iter.times, "\n")
  
  # Stage 2
  cat("Stage 2 starts... \n", " First GWAS ")
  result2.1 <- Update_beta_noX(Pvalue[, 1], result1.1$alpha, 
                               beta0 = -qnorm(1 - result1.1$pi1))
  cat("iterations=", result2.1$iter.times, "\n")
  cat(" Second GWAS ")
  result2.2 <- Update_beta_noX(Pvalue[, 2], result1.2$alpha, 
                               beta0 = -qnorm(1 - result1.2$pi1))
  cat("iterations=", result2.2$iter.times, "\n")
  
  # Stage 3
  cat("Stage 3 starts... \n")
  alpha <- c(result2.1$alpha, result2.2$alpha)
  beta <- c(result2.1$beta, result2.2$beta)
  result3 <- Update_all_noX(Pvalue, alpha, beta, 0)
  cat(" iterations=", result3$iter.times, "\n")
  
  return(list(alpha = result3$alpha, alpha.stage1 = c(result1.1$alpha, result1.2$alpha), 
              alpha.stage2 = c(result2.1$alpha, result2.2$alpha), beta = result3$beta, 
              beta.stage2 = c(result2.1$beta, result2.2$beta), rho = result3$rho, 
              pi11 = result3$pi11, pi10 = result3$pi10, pi01 = result3$pi01, 
              pi00 = result3$pi00, L = result3$L, 
              iter.times.stage1 = c(result1.1$iter.times, result1.2$iter.times), 
              iter.times.stage2 = c(result2.1$iter.times, result2.2$iter.times), 
              iter.times.stage3 = result3$iter.times))
}

Update_alpha <- function(Pvalue, alpha = 0.1, pi1 = 0.1, maxiter = 1e4, tol = 1e-6){
  M <- length(Pvalue)
  
  L <- rep(0, maxiter)
  
  # EM
  for (iter in 1:maxiter){
    # E step
    comp.pos <- pi1*alpha*(Pvalue)^(alpha-1)
    eta.tilde <- comp.pos/(comp.pos + 1 - pi1)
    
    # compute incomplete log likelihood
    L[iter] <- sum(log(comp.pos + 1 - pi1))
    # print(L[iter])
    # cat("Stage 1 ", iter, "\n")
    
    # M step
    alpha <- -sum(eta.tilde)/sum(eta.tilde*log(Pvalue))
    pi1 <- sum(eta.tilde)/M
    
    # check convergence
    if (iter != 1){
      if (L[iter] < L[iter-1]){
        print("L is not increasing")
        break
      }
      if (abs((L[iter]-L[iter-1])/L[iter]) < tol){
        L <- L[1:iter]
        iter.times <- iter
        break
      }
      else
        iter.times <- maxiter
    }
  }
  
  return ( list (alpha = alpha, pi1 = pi1, pi_i1 = eta.tilde, L = L, iter.times = iter.times ))
  
}

Update_beta_noX <- function(Pvalue, alpha, beta0, maxiter = 1e4, tol = 1e-6){
  M <- length(Pvalue)
  L <- rep(0, maxiter)
  
  beta <- beta0
  
  for (iter in 1:maxiter){
    
    # E step
    Phi <- pnorm(beta)
    phi <- dnorm(beta)
    
    comp.pi1 <- Phi*alpha*Pvalue^(alpha - 1)
    comp.pi0 <- 1 - Phi
    comp.L <- comp.pi1 + comp.pi0
    E_eta <- comp.pi1/comp.L
    
    v <- phi*(alpha*Pvalue^(alpha - 1) - 1)/comp.L
    E_Z <- beta + v
    E_Z2 <- beta^2 + beta*v + 1
    
    # Log-likelihood
    L[iter] <- sum(log(comp.L))
    # print(L[iter])
    
    # M step
    gamma <- mean(E_Z) 
    sigma2 <- mean(E_Z2 - 2*E_Z*gamma + gamma^2)
    beta <- gamma/sqrt(sigma2)
    
    alpha <- -sum(E_eta)/sum(E_eta*log(Pvalue))
    
    # check convergence
    if (iter != 1){
      if (L[iter] < L[iter-1]){
        print("L is not increasing")
        break
      }
      if (abs((L[iter]-L[iter-1])/L[iter]) < tol){
        L <- L[1:iter]
        iter.times <- iter
        break
      }
      else
        iter.times <- maxiter
    }
    
  }
  
  return ( list (alpha = alpha, beta = beta, pi_i1 = E_eta, L = L, iter.times = iter.times ))
}

Update_all_noX <- function(Pvalue, alpha, beta, rho, maxiter = 1e4, tol = 1e-6){
  
  M <- nrow(Pvalue)
  L <- rep(0, maxiter)
  
  
  for (iter in 1:maxiter){
    # E step
    
    # E_eta
    Phi <- pnorm(beta)
    Phi11 <- pbivnorm(beta[1], beta[2], rho)
    Phi10 <- -Phi11 + Phi[1]
    Phi01 <- -Phi11 + Phi[2]
    Phi00 <- 1 + Phi11 - Phi[1] - Phi[2]
    
    comp.pi11 <- Phi11*alpha[1]*Pvalue[, 1]^(alpha[1] - 1)*alpha[2]*Pvalue[, 2]^(alpha[2] - 1)
    comp.pi10 <- Phi10*alpha[1]*Pvalue[, 1]^(alpha[1] - 1)
    comp.pi01 <- Phi01*alpha[2]*Pvalue[, 2]^(alpha[2] - 1)
    comp.pi00 <- Phi00
    comp.pi   <- comp.pi11 + comp.pi10 + comp.pi01 + comp.pi00
    
    pi11 <- comp.pi11/comp.pi
    pi10 <- comp.pi10/comp.pi
    pi01 <- comp.pi01/comp.pi
    pi00 <- comp.pi00/comp.pi
    
    E_eta1 <- pi11 + pi10
    E_eta2 <- pi11 + pi01
    
    # E_Z
    c <- 1/sqrt(1 - rho^2)
    
    phi <- dnorm(beta)
    
    red <- phi[1]*pnorm((beta[2] - rho*beta[1])*c)
    green <- phi[2]*pnorm((beta[1] - rho*beta[2])*c)
    
    comp1 <- red + rho*green
    comp2 <- green + rho*red
    comp3 <- -beta[1]*red - rho^2*beta[2]*green + 
      rho/c*phi[2]*dnorm((beta[1] - rho*beta[2])*c)
    comp4 <- -beta[2]*green - rho^2*beta[1]*red + 
      rho/c*phi[1]*dnorm((beta[2] - rho*beta[1])*c)
    comp5 <- -rho*beta[1]*red - rho*beta[2]*green +
      phi[1]/c*dnorm((beta[2] - rho*beta[1])*c)
    
    E_Z1 <- beta[1] + pi11*comp1/(Phi11 + (Phi11 == 0)) + 
      pi10*(phi[1] - comp1)/(Phi10 + (Phi10 == 0)) +
      pi01*(rho*phi[2] - comp1)/(Phi01 + (Phi01 == 0)) + 
      pi00*(-phi[1] - rho*phi[2] + comp1)/(Phi00 + (Phi00 == 0))
    E_Z2 <- beta[2] + pi11*comp2/(Phi11 + (Phi11 == 0)) + 
      pi10*(rho*phi[1] - comp2)/(Phi10 + (Phi10 == 0)) +
      pi01*(phi[2] - comp2)/(Phi01 + (Phi01 == 0)) + 
      pi00*(-phi[2] - rho*phi[1] + comp2)/(Phi00 + (Phi00 == 0))
    E_Z <- cbind(E_Z1, E_Z2)
    E_Z1_2 <- beta[1]*(2*E_Z1 - beta[1]) + 1 +
      pi11*comp3/(Phi11 + (Phi11 == 0)) + 
      pi10*(-beta[1]*phi[1] - comp3)/(Phi10 + (Phi10 == 0)) +
      pi01*(-rho^2*beta[2]*phi[2] - comp3)/(Phi01 + (Phi01 == 0)) +
      pi00*(beta[1]*phi[1] + rho^2*beta[2]*phi[2] + comp3)/(Phi00 + (Phi00 == 0))
    E_Z2_2 <- beta[2]*(2*E_Z2 - beta[2]) + 1 +
      pi11*comp4/(Phi11 + (Phi11 == 0)) + 
      pi10*(-rho^2*beta[1]*phi[1] - comp4)/(Phi10 + (Phi10 == 0)) +
      pi01*(-beta[2]*phi[2] - comp4)/(Phi01 + (Phi01 == 0)) +
      pi00*(beta[2]*phi[2] + rho^2*beta[1]*phi[1] + comp4)/(Phi00 + (Phi00 == 0))
    E_Z1_Z2 <- beta[1]*beta[2] + beta[1]*(E_Z2 - beta[2]) + 
      beta[2]*(E_Z1 - beta[1]) + rho + pi11*comp5/(Phi11 + (Phi11 == 0)) +
      pi10*(-rho*beta[1]*phi[1] - comp5)/(Phi10 + (Phi10 == 0)) +
      pi01*(-rho*beta[2]*phi[2] - comp5)/(Phi01 + (Phi01 == 0)) +
      pi00*(rho*(beta[1]*phi[1] + beta[2]*phi[2]) + comp5)/(Phi00 + (Phi00 == 0))
    
    # calculate L
    L[iter] <- sum(log(comp.pi))
    # print(L[iter])
    # print(alpha)
    # print(rho)
    # print(beta)
    
    # M step
    alpha[1] <- -sum(E_eta1)/sum(E_eta1*log(Pvalue[, 1]))
    alpha[2] <- -sum(E_eta2)/sum(E_eta2*log(Pvalue[, 2]))
    gamma <- as.matrix(as.numeric(colMeans(E_Z)))
    sum_EZZ <- matrix(c(sum(E_Z1_2), sum(E_Z1_Z2), sum(E_Z1_Z2), sum(E_Z2_2)), 2, 2)
    Sigma <- (sum_EZZ - 2 * gamma %*% colSums(E_Z) +
                M*gamma %*% t(gamma))/M
    beta <- gamma/sqrt(diag(Sigma))
    rho <- as.numeric(Sigma[1, 2]/sqrt(Sigma[1, 1])/sqrt(Sigma[2, 2]))
    
    # check convergence
    if (iter != 1){
      if (L[iter] < L[iter-1]){
        print("L is not increasing")
        break
      }
      if (abs((L[iter]-L[iter-1])/L[iter]) < tol){
        L <- L[1:iter]
        iter.times <- iter
        break
      }
      else
        iter.times <- maxiter
    }
    
  }
  
  return(list(alpha = alpha, beta = beta, rho = rho, pi11 = pi11, pi10 = pi10, pi01 = pi01,
              pi00 = pi00, L = L, iter.times = iter.times))
  
}
