# Preamble ----------------------------------------------------------------
require(Rcpp) # For C++ Gibbs sampler using RcppArmadillo
require(MASS) # For the mvrnorm function
require(plyr) # For its easy join function to join datasets (SQL-style)



# R Bayes Regression Function ---------------------------------------------------------

rbr <- function(iter = 999, y, X, prior.a, prior.b, prior.mean, prior.prec){
  # Sample size and number of regressors
  n <- length(y)
  k <- ncol(X)
  
  
  
  # Prior mean vector for betas
  if (length(prior.mean)<k){
    prior.mu <- rep(prior.mean[1], k)
  } else {
    prior.mu <- prior.mean
  }
  # Prior precision matrix for betas
  #prior.tau <- matrix(rep(prior.prec, k*k), ncol = k)
  if (length(prior.prec)<k){
    prior.tau <- diag(rep(prior.prec[1],k))
  } else {
    prior.tau <- diag(prior.prec)
  }
  
  # Posteriors for a and b
  v <- solve(t(X)%*%X)
  bhat <- v%*%t(X)%*%y
  posterior.b <- prior.b + t(y-X%*%bhat)%*%(y-X%*%bhat)/2
  posterior.a <- prior.a + n/2
  
  # Posterior for mu
  posterior.mu <- solve(t(X)%*%X+prior.tau)%*%(prior.tau%*%prior.mu + t(X)%*%y)
  
  # Monitor parameters
  var <- c(1, numeric(iter))
  beta <- matrix(numeric(k*(iter+1)), ncol = k, byrow = T)
  rsq <- c(1, numeric(iter))
  
  
  # Start the MCMC procedure
  for (i in 1:iter){
    # Betas from multivariate normal
    posterior.Sigma <- var[i] * solve(t(X)%*%X + prior.tau)
    beta[i+1,] <- mvrnorm(1, posterior.mu, posterior.Sigma)
    
  
    # Variance from inverse gamma
    var[i+1] <- 1/rgamma(1, posterior.a, posterior.b)
    
    # Calculate Rsquared
    # Residual Variance
    varres <- var(y-X%*%beta[i,])
    rsq[i+1] <- 1-varres/var(y)
    
  }
  
  result <- as.data.frame(cbind(var,beta,rsq))
  colnames(result) <- c("Variance",paste("b", rep(0:(k-1)), sep = ""), 
                        "Rsquared")
  
  return(result)
}


# C++ Bayes Regression Function -------------------------------------------------------

# Because of speed considerations, I wanted to try to write the Bayes Regression
# in C++ via the Rcpp package. The sampler code is in the file cppbr.cpp.
# I make a lot of use of the linear algebra library Armadillo, and I tried to 
# turn the different operations into matrix algebra as much as possible. The 
# code below compiles the C++ code and makes it available as a "regular"
# function in R for interfacing purposes.

# Compile the C++ Bayes Regression
Rcpp::sourceCpp("cppbr.cpp")


# R Gibbs Sampler ---------------------------------------------------------
rgs <- function(iter = 999, y, X, prior_a, prior_b, prior_mean, prior_prec,
                inits = NULL, verbose = T){
  # Sample size and number of regressors
  n <- length(y)
  k <- ncol(X)
  
  # Create prior mean vector
  if (!is.vector(prior_mean) || length(prior_mean) != k){
    prior_mu <- rep(prior_mean[1], k)
  } else {
    prior_mu <- prior_mean
  }
  
  # Create prior precision vector
  if (!is.vector(prior_prec) || length(prior_prec) != k){
    prior_tau <- rep(prior_prec[1], k)
  } else {
    prior_tau <- prior_prec
  }
  
  # Monitor parameters
  var <- numeric(iter+1)
  beta <- matrix(numeric(k*(iter+1)), ncol = k)
  rsq <- c(numeric(iter+1))
  
  
  # Set initial values
  if (is.null(inits) || length(inits) != k + 1 || !is.numeric(unlist(inits))
      || is.null(inits$var)){
    # Calculate initial values if inits have been incorrectly specified
    # or not specified at all.
    if (verbose == T) {cat("Beta inits assumed ols estimates, variance 1. \n")}
    var[1] <- 1
    beta[1,] <- solve(t(X)%*%X)%*%t(X)%*%y
    
  } else {
    # If inits have been correctly specified, use them.
    if (verbose == T) {
      cat("Inits: |", paste(names(unlist(inits)),": ",unlist(inits), " |", 
                            sep = ""),"\n")
    }
    var[1] <- inits$var # assign var init to var
    beta[1,] <- unlist(inits[- which(names(inits) == "var")]) # assign remaining
                                                              # ordered to beta
  }
  
  # Set bhat, the current estimate of the betas, to the initial values (for
  # separate updating of the betas)
  bhat <- beta[1,]
  
  # Start the MCMC procedure
  for (i in 1:iter){
    
    # sample conditional posterior for variance
    
    # Posterior alpha
    p_a <- prior_a + n/2
    
    # Posterior residuals
    p_res <- y - (X %*% beta[i,])
    
    # Posterior beta
    p_b <- prior_b + (t(p_res)%*%p_res)/2
    
    # Sample
    var[i+1] <- 1/rgamma(1, shape = p_a, rate = p_b)
    
    
    
    # sample conditional posterior for b's
    for (v in 1:k){
      
      # Posterior Variance
      p_tau <- 1/(t(X[,v])%*%X[,v]/var[i+1] + (prior_tau[v]^2))
      
      # Posterior Mean
      p_mu <- (prior_mu[v]*(prior_tau[v]^2) +
                 (y-bhat[-v]%*%t(X[,-v]))%*%X[,v]/var[i+1]) * p_tau
      
      # Sample conditional posterior for beta[v]
      bhat[v] <- rnorm(1, mean = p_mu, sd = sqrt(p_tau))
      
    }
    
    # Store betas from current bhat
    beta[i+1,] <- bhat
    
    # Calculate Rsquared
    # Residual Variance
    varres <- var(y-X%*%beta[i,])
    rsq[i+1] <- 1-varres/var(y)
    
  }
  
  result <- as.data.frame(cbind(var,beta,rsq))
  colnames(result) <- c("Variance",paste("b", rep(0:(k-1)), sep = ""), 
                        "Rsquared")
  
  return(result)
}


# R Metropolis-Hastings Sampler -------------------------------------------

rmhs <- function(iter = 999, y, X, prior_a, prior_b, prior_beta,
                 inits = NULL, dtuning = F, verbose = T){
  # Sample size and number of regressors
  n <- length(y)
  k <- ncol(X)
  
  # Monitor parameters
  var <- numeric(iter+1)
  beta <- matrix(numeric(k*(iter+1)), ncol = k)
  rsq <- c(numeric(iter+1))
  
  # Set scale parameter of proposal distribution of each predictor to 2.38
  tune <- rep(2.38,k)
  
  # For the gibbs-steps we need args_b and fun_prior_b
  fun_prior_b <- unlist(strsplit(prior_beta, "(", 
                                 fixed = T))[seq(1,2*length(prior_beta), 2)]
  args_b <- lapply(strsplit(gsub("[()[:alpha:]]", "", prior_beta), ","),
                   as.numeric)
  
  # Parse prior_beta to add x in the middle. This way we can create an 
  # expression which can be evaluated by eval() in the MH-step
  prior_beta <- unlist(lapply(prior_beta,function(x) {
    paste(
      append(
        unlist(
          strsplit(x, "(", fixed = T)
        ), 
        "(x,", after = 1), 
      collapse = "")
  }))
  
  
  # Set initial values
  if (is.null(inits) || length(inits) != k + 1 || !is.numeric(unlist(inits))
      || is.null(inits$var)){
    # Calculate initial values if inits have been incorrectly specified
    # or not specified at all.
    if (verbose == T) {cat("Beta inits assumed ols estimates, variance 1. \n")}
    var[1] <- 1
    beta[1,] <- solve(t(X)%*%X)%*%t(X)%*%y
    
  } else {
    # If inits have been correctly specified, use them.
    if (verbose == T) {
      cat("Inits: |", paste(names(unlist(inits)),": ",unlist(inits), " |", 
                            sep = ""),"\n")
    }
    var[1] <- inits$var # assign var init to var
    beta[1,] <- unlist(inits[- which(names(inits) == "var")]) # assign remaining
    # ordered to beta
  }
  
  # Set bhat, the current estimate of beta, to the inits. (for separate
  # updating of the betas)
  bhat <- beta[1,]
  
  # Start the MCMC procedure
  for (i in 1:iter){
    
    # sample conditional posterior for variance
    
    # Posterior alpha
    p_a <- prior_a + n/2
    
    # Posterior residuals
    p_res <- y - (X %*% beta[i,])
    
    # Posterior beta
    p_b <- prior_b + (t(p_res)%*%p_res)/2
    
    # Sample
    var[i+1] <- 1/rgamma(1, shape = p_a, rate = p_b)
    
    
    
    # sample conditional posterior for b's using MH OR Gibbs procedure
    for (v in 1:k){
      
      if (fun_prior_b[v] == "dnorm"){
        # Gibbs procedure if prior is conjugate
        # Posterior Variance
        p_tau <- 1/(t(X[,v])%*%X[,v]/var[i+1] + (1/args_b[[v]][2]^2))
        
        # Posterior Mean
        p_mu <- (args_b[[v]][1]/(args_b[[v]][2]^2) +
                   (y-bhat[-v]%*%t(X[,-v]))%*%X[,v]/var[i+1]) * p_tau
        
        # Sample conditional posterior for beta[v]
        bhat[v] <- rnorm(1, mean = p_mu, sd = sqrt(p_tau))
        
      } else {
        # Random Walk MH procedure if this is not the case
        
        # pdf for the target distribution
        # Likelihood Variance
        l_tau <- var[i+1]/t(X[,v])%*%X[,v]
        # Likelihood Mean
        l_mu <- (y-bhat[-v]%*%t(X[,-v]))%*%X[,v]/var[i+1] * l_tau
        # posterior <- prior * likelihood
        fx <- function(x){eval(parse(text=prior_beta[v])) *
            dnorm(x, mean = l_mu, sd = sqrt(l_tau))}
        
        # calculate the maximum of this function to serve as mean of proposal
        #propmean <- optimize(f = fx, maximum = T, lower = min(beta[,v])-100, 
        #upper = max(beta[,v])+100)$maximum
        
        #qx <- function(x) beta[i,v]+dnorm(x, 0, tune)
        
        # Sample candidate value from OLS proposal distribution
        bstar <- bhat[v]+rnorm(1, 0, tune[v])
        
        # Compute the acceptance ratio
        #r <- min(1, fx(bstar)/fx(beta[i,v])*qx(beta[i,v])/qx(bstar))
        r <- min(1,fx(bstar)/fx(beta[i,v]))
        
        # Assign when accepted
        if (runif(1) <= r){
          bhat[v] <- bstar
        }
        
        # Dynamic tuning procedure using information from
        # https://support.sas.com/documentation/cdl/en/statug/63033/HTML/default
        # /viewer.htm#statug_mcmc_sect022.htm
        
        # Dynamic tuning after every [dtuning]th iteration
        if (dtuning != F && i%%dtuning == 0){
          # Set the optimal acceptance rate
          popt <- 0.45
          
          # Calculate current acceptance rate
          pcur <- 1-mean(beta[1:i,v]==beta[2:(i+1),v])
          
          # If acceptance rate is too high, adjust scale parameter to be more 
          # wide. This will decrease the acceptance rate. If the acceptance rate
          # is too low, adjust the scale parameter downwards.
          if (pcur > 0.5 || pcur < 0.15){
            tune[v] <- (tune[v]*(1/pnorm(popt/2)))/(1/pnorm(pcur/2))
          }
        }
      }
    }
    
    # set beta to bhat
    beta[i+1,] <- bhat
    
    # Calculate Rsquared
    # Residual Variance
    varres <- var(y-X%*%beta[i,])
    rsq[i+1] <- 1-varres/var(y)
    
  }
  
  cat("Acceptance Rates:", round(apply(beta,2,function(x)1-mean(x[-iter]==x[-1])),3),"\n")
  cat("Tuning Parameter:", round(tune,3),"\n")
  
  result <- as.data.frame(cbind(var,beta,rsq))
  colnames(result) <- c("Variance",paste("b", rep(0:(k-1)), sep = ""), 
                        "Rsquared")
  
  return(result)
}
