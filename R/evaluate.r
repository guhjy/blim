# Evaluation functions that take blimfit-objects and perform operations on them
# Functions: Summary, Bayes Factor, Convergence Plots, Autocorrelation Plots

# Erik-Jan van Kesteren


# Summary -----------------------------------------------------------------

summary.blimfit <- function(x,...){
  print(x$summary)
}

# Bayes Factor ------------------------------------------------------------

BF <- function(blimfit, model = "par[1] < par[2]", complement = F){
  
  # Check if object has class blimfit
  if (class(blimfit) != "blimfit") stop("Please enter a blimfit object!")
  
  # Remove unnecessary information from the posterior distributions (the
  # variance and R-squared columns)
  trace <- blimfit$trace[,-c(1,ncol(blimfit$trace))]
  
  # Create evaluable expressions for priors
  p <- lapply(lapply(lapply(strsplit(blimfit$priors, "(", fixed = T),
                            function(x)append(x,"(nrow(trace),", after = 1)),
                     function(x)paste(x,collapse="")),
              function(x)paste("r",substr(x,2,nchar(x)), sep = ""))
  
  # Create a matrix of prior information NOTE that with few iterations this is
  # just a rough approximation as this random sequence is NOT taken from the 
  # Markov chain directly. Asymptotically, this becomes more precise.
  prior <- matrix(unlist(lapply(p,function(x)eval(parse(text=x)))),
                  nrow = nrow(trace))
  

  
  # Evaluate fit of the hypothesis
  fa <- mean(apply(trace, 1, function(par) eval(parse(text=model))))
  
  # Evaluate complexity of the hypothesis
  ca <- mean(apply(prior, 1, function(par) eval(parse(text=model))))
  
  if (complement == F){
    BFau <- fa/ca 
    # Bayes factor against the unconstrained hypothesis
    return(c(fit = fa, complexity = ca, BFau = BFau))
  } else {
    fc <- mean(apply(trace, 1, function(par){
      eval(parse(text=paste("!(",model,")", sep = "")))
    }))
    cc <- mean(apply(prior, 1, function(par){
      eval(parse(text=paste("!(",model,")", sep = "")))
    }))
    BFac <- (fa/ca)/(fc/cc)
    return(c(fit_a = fa, compl_a = ca, fit_c = fc, compl_c = cc, BFac = BFac))
  }
  
}


# Posterior Predictive Check ----------------------------------------------

DWPPC <- function(blimfit, fast = TRUE){
  
  # Check if object has class blimfit
  if (class(blimfit) != "blimfit") stop("Please enter a blimfit object!")
  
  # Remove unnecessary information from the posterior distributions (the
  # R-squared column)
  if (fast == F){
    trace <- blimfit$trace[,-c(ncol(blimfit$trace))]
    cat("Precise PPP calculation. This might take a while! \n")
  } else {
    trace <- blimfit$trace[sample(nrow(blimfit$trace),1000),
                           -c(ncol(blimfit$trace))]
    cat("Fast PPP calculation. Only 1000 iterations evaluated.
Turn off the option 'fast' to perform a precise PPP calculation. \n")
  }
  
  
  
  # Calculate the observed discrepancy measure for each column of the trace
  obsD <- apply(trace[,-1],1,function(b){
    # Calculate discrepancy measure
    resid <- blimfit$y-blimfit$X%*%b
    auto <- resid[-1]-resid[-length(resid)]
    DW <- (t(auto)%*%auto)/(t(resid)%*%resid)
    return(abs(DW-2))
    })
  
  # Calculate the replicated discrepancy measure for each column of the trace
  repD <- apply(trace,1,function(trace){
    # Assign proper parameters
    var <- trace[1]
    b <- trace[-1]
    # Generate data under model
    y <- apply(blimfit$X%*%b,1,function(mean){rnorm(1,mean,sqrt(var))}) 
    
    # Calculate discrepancy measure
    resid <- y-blimfit$X%*%b
    auto <- resid[-1]-resid[-length(resid)]
    DW <- (t(auto)%*%auto)/(t(resid)%*%resid)
    return(abs(DW-2))
  })
  
  # PPP is how often the observed D is smaller than replicated D
  return(mean(obsD < repD))

}


# Convergence Plots -------------------------------------------------------


cplots <- function(blimfit, params = NULL, limit = TRUE){
  
  # Check if object has class blimfit
  if (class(blimfit) != "blimfit") stop("Please enter a blimfit object!")
  
  # Get the trace
  trace <- blimfit$trace
  
  if (limit == T && nrow(trace)>10000) {
    thin <- round(nrow(trace)/10000)
    select <- (1:nrow(trace))[seq(thin,nrow(trace), length.out = 10000)]
    cat("Trace has",nrow(trace),"samples. Only 10000 shown.
Turn off the option 'limit' to show all observations (slow!).")
  } else {
    select <- 1:nrow(trace)
  }
  
  
  # If parameters not given assume all params.
  if (!is.numeric(params)){
    params <- 1:ncol(trace)
  }
  
  # Store and set margins
  margins <- par("mar")
  par(mar = c(2.1,4.1,0.6,1.1))
  
  # Store and set columns & rows
  if (length(params) < 5){
    par(mfrow = c(length(params),2))
  } else {
    par(mfrow = c(4,2))
  }
  
  for (i in params){
    plot(select, trace[select,i], type = "l", col = "blue", main = "", 
         xlab = "", ylab = colnames(trace)[i])
    hist(trace[select,i], breaks = "FD", xlab = "", ylab = "frequency", main = "", 
         col = "dark green", border = "dark green")
    abline(v = mean(trace[,i]), col = "white")
  }
  
  par(mfrow = c(1,1), mar = margins)
  
}


# Autocorrelation Plots ---------------------------------------------------

aplots <- function(blimfit, params = NULL, limit = TRUE){
  
  # Check if object has class blimfit
  if (class(blimfit) != "blimfit") stop("Please enter a blimfit object!")
  
  # MH warning message
  if (blimfit$algorithm == "R based Metropolis-Hastings Sampling Regression"){
    cat("Note that Metropolis-Hastings sampling leads to high autocorrelation!
Check that acceptance rates are between 0.15 and 0.5 in console at blim output.",
        "\n")
  }
  
  # Pick up trace for the column names
  trace <- blimfit$trace[1:2,]
  
  # Pick up autocorrelation
  auto <- blimfit$auto
  
  
  # If parameters not given assume all params.
  if (!is.numeric(params)){
    params <- 1:ncol(trace)
  }
  
  # Store and set margins
  margins <- par("mar")
  par(mar = c(2.1,4.1,0.6,1.1))
  
  # Store and set columns & rows
  if (length(params) < 5){
    par(mfrow = c(length(params),1))
  } else {
    par(mfrow = c(4,1))
  }
  
  for (i in params){
    plot(auto[i,], type = "h", col = "blue", main = "", 
         xlab = "", ylab = colnames(trace)[i], lwd = 3, ylim = c(0,1))
  }
  
  par(mfrow = c(1,1), mar = margins)
  
}
