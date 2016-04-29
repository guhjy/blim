\name{cplots}
\alias{cplots}

\title{
Convergence Plots Function
}
\description{
Display Convergence plots in the graphics device. Both trace and density are 
shown.
}
\usage{
cplots(blimfit, params = NULL, limit = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blimfit}{ a blimfit object from blim }
  \item{params}{ a vector of indices for the parameters to be displayed }
  \item{limit}{ a boolean value indicating whether limiter should be applied. 
  for performance, a standard limit of 10000 random samples from the trace is 
  used for the plots.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Plot output in the graphics device.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Erik-Jan van Kesteren
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{blim}}, \code{\link{aplots}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ convergence }
\keyword{ plot }% __ONLY ONE__ keyword per line
