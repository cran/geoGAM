\name{bootstrap.geoGAM}
\alias{bootstrap.geoGAM}
\alias{bootstrap}
\alias{bootstrap.default}

\title{Bootstrapped predictive distribution

}
\description{

Method for class \code{geoGAM} to compute model based bootstrap for point predictions. Returns complete predictive distribution of which prediction intervals can be computed.
}
\usage{
\method{bootstrap}{default}(object, ...)

\method{bootstrap}{geoGAM}(object, newdata, R = 100,
          back.transform = c("none", "log", "sqrt"),
          seed = NULL, cores = detectCores(), ...)

}

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

Soil properties are predicted for new locations \eqn{\mathbf{s_{+}}}{s+} from the final \code{\link{geoGAM}} fit by \eqn{\tilde{Y}(\mathbf{s_+})=\hat f(\mathbf{x}\mathbf{(s_+)})}{Y(s+) = f(x(s+))}, see function \code{\link{predict.geoGAM}}.
%
	% bootstraped intervals
	To model the predictive distributions for continuous responses \code{\link{bootstrap.geoGAM}} uses a
	non-parametric, model-based bootstrapping approach (\cite{Davison and Hinkley 1997}, pp. 262, 285) as follows:
	%
\enumerate{

		\item New values of the response are simulated according to \eqn{Y(\mathbf{s})^{*} = \hat f(\mathbf{x}(\mathbf{s}))+\mathbf{\epsilon}}{Y(s)* = f(x(s)) + epsilon}, where \eqn{\hat f(\mathbf{x}(\mathbf{s}))}{f(x(s))} are the fitted values of the final model
		and \eqn{\epsilon}{epsilon} are errors randomly sampled with replacement from the centred, homoscedastic residuals of the final model \cite{Wood 2006}, p. 129).

		\item \code{\link{geoGAM}} is fitted to \eqn{Y(\mathbf{s})^{*}}{Y(s)*}.

		\item Prediction errors are computed according to \eqn{\delta_{+}^{*} =  \hat f(\mathbf{x}(\mathbf{s_{+}}))^{*} - (\, \hat f(\mathbf{x}(\mathbf{s_{+}})) + \mathbf{\epsilon} \,)}{delta+* = f(x(s+))* - (f(x(s+))+epsilon)},
		where \eqn{\hat f(\mathbf{x}(\mathbf{s_{+}}))^{*}}{f(x(s+))*} are predicted values at new locations \eqn{\mathbf{s_{+}}}{s+} of the
		model built with the simulated response \eqn{Y(\mathbf{s})^{*}}{Y(s)*} in step B above, and the errors \eqn{\epsilon}{epsilon} are again randomly sampled from the centred, homoscedastic residuals of the final model (see step A).
}

	Prediction intervals are computed according to
	%
	\deqn{
	[\hat f(\mathbf{x}(\mathbf{s_{+})}) - \delta_{+\,(1-\alpha)}^{*}\,; \hat f(\mathbf{x}\mathbf{(s_{+}})) - \delta_{+\,(\alpha)}^{*}]}{
	[f(x(s+)) - delta+*_(1-alpha); f(x(s+)) - delata+*_(alpha)]
	}
	%
	where \eqn{\delta_{+\,(\alpha)}^{*}}{delta+*_(alpha)} and \eqn{\delta_{+\,(1-\alpha)}^{*}}{delta+*_(1-alpha)} are the \eqn{\alpha}{alpha}- and \eqn{(1-\alpha)}{(1-alpha)}-quantiles of \eqn{\delta_{+}^{*}}{delta+*}, pooled over all 1000 bootstrap repetitions.

	Predictive distributions for binary and ordinal responses are
	directly obtained from a final \code{\link{geoGAM}} fit by predicting probabilities
	of occurrence \eqn{\mathrm{\widetilde{Prob}}(Y(\mathbf{s})=r\,|\,\mathbf{x}(\mathbf{s)})}{Prob((s) = r|x(s))}
	(\cite{Davison and Hinkley 1997}, p. 358).

}
\value{

Data frame of \code{nrows(newdata)} rows and \code{R + 2} columns with \code{x} and \code{y} indicating coordinates of the location and \code{P1} to \code{P...R} the prediction at this location from \code{1...R} replications.
}
\references{

Nussbaum, M., Walthert, L., Fraefel, M., Greiner, L., and Papritz, A.: Mapping of soil properties at high resolution in Switzerland using boosted geoadditive models, SOIL, 3, 191-210, doi:10.5194/soil-3-191-2017, 2017.

Davison, A. C. and Hinkley, D. V., 2008. Bootstrap Methods and Their Applications. Cambridge University Press.
}
\author{

M. Nussbaum
}
\seealso{

To create geoGAM objects see \code{\link{geoGAM}} and to predict without simulation of the predictive distribution see \code{\link{predict.geoGAM}}.
}
\examples{

\donttest{
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
                        max.stop = 20,
                        cores = 1)


## compute model based bootstrap with 10 repetitions (use at least 100)
quakes.boot <- bootstrap(quakes.geogam,
                         newdata = quakes.pred,
                         R = 10, cores = 1)


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
\keyword{ spatial }
\keyword{ models & regression & nonlinear }
