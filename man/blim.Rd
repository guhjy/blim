\name{blim}
\alias{blim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Linear Model
}
\description{
Fit a univariate Bayesian linear regression model.
}
\usage{
blim(formula, data, iter = 9999, burnin = 50, inits = NULL, thin = 0
     prior_tau = "dgamma(0.001,0.001)", prior_b = "dnorm(0,10000)",
     method = "cppbr", mtsprior = F, verbose = T, dtuning = F)

\method{summary}{blim}(blimfit, ...)
}
%- maybe also 'usage' for other objects documented here.
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
\item{method}{ cppbr = Bayes Regression in C++, rbr = Bayes Regression in R,
               rgs = Bayesian Regression under Gibbs Sampling,
               rmhs = Bayesian Regression under Gibbs Sampling with 
               Metropolis-Hastings step for Nonconjugate Priors }
\item{mtsprior}{ Boolean value indicating whether a Minimum Training Sample
                 prior should be calculated. For use in conjuction with ANOVA &
                 Bayes Factors for informative hypotheses }
\item{verbose}{ print all info in the console or not (in development!) }
\item{dtuning}{ dynamic tuning of the random walk Metropolis-Hastings step }
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{blimfit}, including all necessary information such as
summary, trace and original data.
}
\references{
Beland, S., Klugkist, I., Raiche, G. & Magis, D. (2012). A short introduction 
into Bayesian evaluation of informative hypotheses as an alternative to 
exploratory comparisons of multiple group means. Tutorials in Quantitative 
Methods for Psychology 8(2), p. 122-126.

Raftery, A. E. (1996). Hypothesis testing and model selection. Markov chain 
Monte Carlo in practice, 163-188.
}
\author{
Erik-Jan van Kesteren
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cplots}}, \code{\link{aplots}}, \code{\link{BF}}, 
\code{\link{PPC}}, \code{\link{ICBF}}, \code{\link{DIC}}
}
\examples{
# Most simple example with full functionality
fit1 <- blim(dist~speed, data=cars)

summary(fit1)

# More complex example with nonconjugate priors
fit2 <- blim(formula = dist~speed, data = cars, iter = 9999, burnin = 100,
             inits = list(var = 200, b0 = -15, b1 = 3), thin = 0,
             prior_tau = "dgamma(0.001,0.001)",
             prior_b = c("dcauchy(0,2000)", "dunif(3,4.5)"),
             method = "rmhs", dtuning = 1000)

summary(fit2)

# check convergence and autocorrelation
cplots(fit1)
aplots(fit1)

# calculate Bayes Factor of model 2 versus model 1
BF(fit2, fit1)

# perform posterior predictive check for normality of residuals
PPC(fit2)

# evaluate informative hypothesis in model 1
ICBF(fit1, model = "par[1] < par[2]", complement = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ bayesian }
\keyword{ linear }
\keyword{ modelling }% __ONLY ONE__ keyword per line
