\name{bootstrap.geoGAM}
\alias{bootstrap.geoGAM}
\alias{bootstrap}
\alias{bootstrap.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrapped predictive distribution
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Method for class \code{geoGAM} to compute model based bootstrap for point predictions. Returns complete predictive distribution of which prediction intervals can be computed.
}
\usage{
\method{bootstrap}{default}(object, ...)

\method{bootstrap}{geoGAM}(object, newdata, R = 100,
          back.transform = c("none", "log", "sqrt"),
          seed = NULL, cores = detectCores(), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{geoGAM object}
  \item{newdata}{data frame in which to look for covariates with which to predict.}
  \item{R}{number of bootstrap replicates, single positive integer.}
  \item{back.transform}{sould to \code{log} or \code{sqrt} transformed responses unbiased back transformation be applied? Default is \code{none}.}
  \item{seed}{seed for simulation of new response. Set seed for reproducible results.}
  \item{cores}{number of cores to be used for parallel computing.}
  \item{...}{further arguments.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~

To obtain predictive distribution for continuous responses \code{bootstrap} implements a model based bootstrap approach  (\cite{Davison-Hinkley-2008, pp. 262,	285}) for \code{geoGAM} models. Errors are simulated from Gaussian distribution with \code{R} repetitions involving the following steps:

	\enumerate{
		\item simulating new response \eqn{Y(s)^{*}}{Y(s)*} from Gaussian
		distribution \eqn{\mathcal{N} (\hat f(\mathbf{x},s), \hat \sigma^2)}{NN(f(x,s), sigma^2)} with
		fitted values of the final model \eqn{\hat f(\mathbf{x}, s)}{f(x,s)} and residual
		variance \eqn{\hat \sigma^2}{sigma^2},

		\item selecting geoadditive model with \code{\link{geoGAM}} for
		\eqn{Y(s)^{*}}{Y(s)*},

		\item computing prediction error

		\deqn{\sigma_{i}^{*} = \hat
		f(\mathbf{x},s)^{*} - ( \hat f(\mathbf{x},s) + \mathcal{N} (0, \hat
		\sigma^{*2}))}{sigma_i* = f(x,s)* - f(x,s) + NN(0, sigma*^2)}


		with fitted values \eqn{\hat f(\mathbf{x},s)^{*}}{f(x,s)*} and
		residual variance \eqn{\hat \sigma^{*2}}{sigma*^2} of the model built on the
		simulated response \eqn{Y(s)^{*}}{Y(s)*}.
    %
		Simulated point predictions are then computed by \eqn{\hat f(\mathbf{x},s)^{*} - \sigma_{i}^{*}}{f(x,s)* - sigma_i*}.
}

This results in a predictive distribution for each site. For \code{back.transform} = \code{log} or \code{sqrt} the simulated predictions are backtransformed after the steps above (see \code{\link{predict.geoGAM}} for more information).
Prediction intervals can be achieved by computing the desired percentiles (e.g. 2.5 \% and 97.5 \% percentiles to get lower and upper limits of 95 \% prediction intervals.)

}
\value{

Data frame of \code{nrows(newdata)} rows and \code{R + 2} columns with \code{x} and \code{y} indicating coordinates of the location and \code{P1} to \code{P...R} the prediction at this location from \code{1...R} replications.
}
\references{
%% ~put references to the literature/web site here ~
Nussbaum, M., Papritz, A., and others (in prep). Mapping of soil properties at high resolution by
geo-additive models.

Davison, A. C. and Hinkley, D. V. 2008. Bootstrap Methods and Their Applications. Cambridge University Press.
}
\author{
%%  ~~who you are~~
Madlene Nussbaum, \email{madlene.nussbaum@env.ethz.ch}
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
To create geoGAM objects see \code{\link{geoGAM}} and to predict without simulation of the predictive distribution see \code{\link{predict.geoGAM}}.
}
\examples{

\dontrun{
data(quakes)

# group stations to ensure min 20 observations per factor level
# and reduce number of levels for speed
quakes$stations <- factor( cut( quakes$stations, breaks = c(0,15,19,23,30,39,132)) )

# Artificially split data to create prediction data set
set.seed(1)
quakes.pred <- quakes[ ss <- sample(1:nrow(quakes), 500), ]
quakes <- quakes[ -ss, ]

quakes.geogam <- geoGAM(response = "mag",
                        covariates = c("stations", "depth"),
                        coords = c("lat", "long"),
                        data = quakes,
                        max.stop = 20)


## compute model based bootstrap with 100 repetitions
quakes.boot <- bootstrap(quakes.geogam,
                         newdata = quakes.pred,
                         R = 100)


# plot predictive distribution for site in row 9
hist( as.numeric( quakes.boot[ 9, -c(1:2)] ), col = "grey",
      main = paste("Predictive distribution at", paste( quakes.boot[9, 1:2], collapse = "/" )),
      xlab = "predicted magnitude")

# compute 95 \% prediction interval and add to plot
quant95 <- quantile( as.numeric( quakes.boot[ 9, -c(1:2)] ), probs = c(0.025, 0.975) )
abline(v = quant95[1], lty = "dashed")
abline(v = quant95[2], lty = "dashed")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }
\keyword{ models & regression & nonlinear }% __ONLY ONE__ keyword per line
