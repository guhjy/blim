# Bayesian linear model function blm is a wrapper for the different functions
# I wrote: rbr (Bayesian Regression in R), cppbr (Bayesian Regression in C++),
# rgs (Gibbs Sampling in R), rmhs (Metropolis-Hastings Sampling in R). 
# It uses the same syntax as the lm and glm functions, using a formula
# to specify a model with at most 1 predicted value (univariate).
# It has options to specify method (rbr / cppbr / gs / mhs), 
# verbosity (T / F) to return information, priors in the form "beta(1,1)", as 
# well as several options such as burnin, inits, and iterations. Before running,
# check that the above mentioned functions are loaded in the environment.

# Erik-Jan van Kesteren

# Function Loading --------------------------------------------------------


blim <- function(formula, data, iter = 9999, burnin = 50, inits = NULL, thin = 0,
                prior_tau = "dgamma(0.001,0.001)", prior_b = "dnorm(0,10000)", 
                method = "cppbr", mtsprior = F, verbose = T, dtuning = F){
  
  # Extract response vector y and data matrix X
  mf <- model.frame(formula,data)
  y <- model.response(mf, "numeric") # assigns the response vector
  X <- model.matrix(formula, mf) # assigns the data matrix
  
  # If mtsprior is indicated, calculate a minimum-training-sample prior for b
  if (mtsprior == T){
    means <- tausq <- numeric(999)
    
    for (i in 1:999){
      v <- sample(y,2,replace = T)
      means[i] <- mean(v)
      tausq[i] <- ((v[1]-means[i])^2+(v[2]-means[i])^2)/4
    }
    
    prior_b <- rep(paste("dnorm(",mean(means),",",(mean(tausq)+var(means)),
                        ")",sep = ""), ncol(X))
    
    cat("MTS-prior calculated:", prior_b[1], "\n")
  }
  
  # Parse functions for priors and check if they exist
  fun_prior_tau <- strsplit(prior_tau, "(", fixed = T)[[1]][1]
  fun_prior_b <- unlist(strsplit(prior_b, "(", 
                                 fixed = T))[seq(1,2*length(prior_b), 2)]
  
  # If not all functions inside the priors for b exist then stop
  if (!all(apply(X = as.matrix(fun_prior_b),
                 MARGIN = 1, 
                 FUN = function(x) existsFunction(x))) || 
      !existsFunction(fun_prior_tau)){
    stop(paste(fun_prior_b, "or", fun_prior_tau, "does not exist!"))
  }
  
  # Parse arguments for the priors
  args_tau <- lapply(strsplit(gsub("[()]", "",
                                   unlist(strsplit(prior_tau, 
                                                   "(", fixed = T))[c(F,T)]), 
                              ","),
                     as.numeric)
  args_b <- lapply(strsplit(gsub("[()]", "", 
                                 unlist(strsplit(prior_b, 
                                                 "(", fixed = T))[c(F,T)]), 
                            ","),
                   as.numeric)
  
  if (length(fun_prior_b)<ncol(X)){
    fun_prior_b <- rep(fun_prior_b[1],ncol(X))
    args_b <- rep(args_b[1],ncol(X))
    prior_b <- rep(prior_b[1],ncol(X))
    cat("Using", prior_b[1], "as prior for all beta coefficients. \n")
  }
  
  if (toupper(method) == "RBR"){
      
    if (fun_prior_tau != "dgamma" || any(fun_prior_b != "dnorm")){
      stop("Method RBR only applicable with conjugate priors")
    }
  
    result <- rbr(iter = iter, y = y, X = X, 
                  prior.a = args_tau[[1]][1], prior.b = args_tau[[1]][2], 
                  prior.mean = unlist(args_b)[c(T,F)], 
                  prior.prec = 1/unlist(args_b)[c(F,T)])
    colnames(result) <- c("Variance", colnames(X), "Rsquared")
    algo <- "R-based Bayesian Linear Regression"
    
  } else if (toupper(method) == "CPPBR") {
    
    if (fun_prior_tau != "dgamma" || any(fun_prior_b != "dnorm")){
      stop("Method CPPBR only applicable with conjugate priors")
    }
    
    result <- cppbr(iter = iter, y = y, X = X, 
                  prior_a = args_tau[[1]][1], prior_b = args_tau[[1]][2], 
                  prior_mu = unlist(args_b)[c(T,F)], 
                  prior_prec = 1/unlist(args_b)[c(F,T)])
    colnames(result) <- c("Variance", colnames(X), "Rsquared")
    algo <- "C++ based Bayesian Linear Regression"
  
  } else if (toupper(method) == "RGS"){ 
    
    if (fun_prior_tau != "dgamma" || any(fun_prior_b != "dnorm")){
      stop("Method RGS only applicable with conjugate priors")
    }
    
    result <- rgs(iter = iter, y = y, X = X, 
                  prior_a = args_tau[[1]][1], prior_b = args_tau[[1]][2], 
                  prior_mean = unlist(args_b)[c(T,F)], 
                  prior_prec = 1/unlist(args_b)[c(F,T)],
                  verbose = verbose, inits = inits)
    colnames(result) <- c("Variance", colnames(X), "Rsquared")
    algo <- "R based Gibbs Sampling Regression"
    
  } else if (toupper(method) == "RMHS"){
    
    if (fun_prior_tau != "dgamma"){
      stop("Method RMHS only applicable with conjugate variance prior")
    }
    
    result <- rmhs(iter = iter, y = y, X = X, 
                   prior_a = args_tau[[1]][1], prior_b = args_tau[[1]][2], 
                   prior_beta = prior_b,verbose = verbose, inits = inits, 
                   dtuning = dtuning)
    
    
    algo <- "R based Metropolis-Hastings Sampling Regression"
  } else {
    stop("Unknown method: ", toupper(method))
  }
  
  
  
  if (verbose == T){cat(algo, "\n")}
  
  # Remove burnin and apply thinning to trace
  trace <- result[-(1:(burnin+1)),]
  if (thin > 0) {trace <- trace[(seq(from = 1, to = nrow(trace), by = thin)),]}
  
  # Create summary statistics
  summary <- cbind(mean = apply(X = trace, MARGIN = 2, 
                                FUN = mean), 
                   s.e. = apply(X = trace, MARGIN = 2, 
                                FUN = sd),
                   cred2.5 = apply(X = trace, MARGIN = 2,
                                 FUN = function(x) quantile(x, 0.025)),
                   cred97.5 = apply(X = trace, MARGIN = 2,
                                 FUN = function(x) quantile(x, 0.975)),
                   mcerror = apply(X = trace, MARGIN = 2,
                                   FUN = function(x) var(x)/sqrt(iter)))
  
  # Create trace autocorrelation matrix
  autocor <- matrix(numeric((ncol(X)+2)*40), ncol = 40)
  
  for (i in 1:40){
    autocor[,i] <- apply(X = trace, MARGIN = 2,
                     function(x) abs(cor(x[-(1:i)], 
                                         x[-((length(x)+1-i):length(x))]))
                     )
  }
  
  # Give the priors for the betas names before the output
  names(prior_b) <- paste("b",0:(ncol(X)-1),sep = "")
  
  blimfit <- list(algorithm = algo, 
                  summary = summary, 
                  auto = autocor, 
                  y = y, X = X, 
                  trace = trace, 
                  priors = prior_b)
  class(blimfit) <- "blimfit"
  
  return(blimfit)
}