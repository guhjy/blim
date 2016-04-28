\name{blim}
\alias{blim}
\title{Bayesian Linear Model}
\description{Fit a univariate Bayesian linear regression model.}
\usage{
blim(formula, data, iter = 9999, burnin = 50, inits = NULL, thin = 0
     prior_tau = "dgamma(0.001,0.001)", prior_b = "dnorm(0,10000)",
     method = "cppbr", mtsprior = F, verbose = T, dtuning = F)

\method{summary}{blim}(blimfit, ...)
}
\arguments{
\item{formula}{ a symbolic description of the model to be fit. }
\item{data}{ an optional data frame containing the variables in the model. }
\item{iter}{ amount of iterations for the sampler to perform and return. }
\item{burnin}{ burn-in amount of iterations ignored in the trace. }
\item{inits}{ initial values for variance and betas in a \code{"list"} format. }
\item{thin}{ amount of thinning to be performed }
\item{prior_tau}{ prior for the residual variance in the linear model, of the
                  form dgamma(shape, scale). Only dgamma allowed! }
\item{prior_b}{ (vector of) priors for beta coefficients. With
                \code{method = "rmhs"}, nonconjugate priors are allowed. Any
                density function works, as long as it has a d and an r method. }
\item{method}{ cppbr = Bayes Regression in C++, rbr = Bayes Regression in R
               rgs = Bayesian Regression under Gibbs Sampling,
               rmhs = Bayesian Regression under Gibbs Sampling with Metropolis-
               Hastings step for Nonconjugate Priors }
\item{mtsprior}{ Boolean value indicating whether a Minimum Training Sample
                 prior should be calculated. For use in conjuction with ANOVA &
                 Bayes Factors for informative hypotheses }
\item{vebose}{ print all info in the console or not }
\item{dtuning}{ dynamic tuning of the random walk Metropolis-Hastings step }
}
\value{
An object of class \code{blimfit}, including all necessary information such as
summary, trace and original data.
}
\author{Erik-Jan van Kesteren}
\examples{
# Most simple example with full functionality
fit1 <- blim(dist~speed, data=cars)
summary(fit1)

# More complex example with nonconjugate priors
fit2 <- blim(formula = dist~speed, data = cars, iter = 50000, burnin = 2000,
             inits = list(var = 200, b0 = -15, b1 = 3), thin = 0,
             prior_tau = "dgamma(0.001,0.001)",
             prior_b = c("dcauchy(0,2000)", "dunif(3,4.5)"),
             method = "rmhs", dtuning = 1000)

summary(fit2)

}
\keyword{regression, bayesian, linear, modelling}
