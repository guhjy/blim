\name{blim-package}
\alias{blim-package}
\docType{package}
\title{
\packageTitle{blim}
}
\description{
\packageDescription{blim}
}


\author{
\packageAuthor{blim}

Maintainer: \packageMaintainer{blim}
}
\references{
Dynamic tuning procedure using information from
https://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_mcmc_sect022.htm
}

\keyword{ package, blim, bayesian, linear, modelling, univariate }

\seealso{

\code{\link[<blim.Rd>]{blim}}

}
\examples{
  # fit the model
  fit1 <- blim(dist~speed, data=cars)

  # check convergence and autocorrelation
  cplots(fit1)
  aplots(fit1)

  # look at results
  summary(fit1)

  # calculate bayes factor of informative hypothesis
  BF(fit1, model = "par[1] > par[2]", complement = TRUE)

  # perform posterior predictive check for autocorrelation of residuals
  DWPPC(fit1)
}
