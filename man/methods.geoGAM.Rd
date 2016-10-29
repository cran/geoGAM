\name{methods}
\alias{summary.geoGAM}
\alias{print.geoGAM}
\alias{plot.geoGAM}
\alias{summary}
\alias{print}
\alias{plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Methods for \code{geoGAM} objects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Methods for models fitted by \code{geoGAM()}.
}
\usage{
\method{summary}{geoGAM}(object, \dots, what = c("final", "path"))

\method{print}{geoGAM}(x, \dots)

\method{plot}{geoGAM}(x, \dots, what = c("final", "path"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{geoGAM}}
  \item{x}{an object of class \code{geoGAM}}
  \item{\dots}{other arguments passed to \code{summary.gam}, \code{plot.gam} or \code{plot.mboost}}
  \item{what}{print summary or plot partial effects of \code{final} selected model or print summary or plot gradient boosting path of model selection \code{path}.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
\code{summary} with \code{what = "final"} calls \code{\link{summary.gam}} to display a summary of the final (geo)additive model. \code{plot} with  \code{what = "final"} calls \code{\link{plot.gam}} to plot partial residual plots of the final model.

\code{summary} with \code{what = "path"} give a summary of covariates selected in each step of model building.
\code{plot} with \code{what = "path"} calls \code{\link{plot.mboost}} to plot the path of the gradient boosting algorithm.

}
\value{

For \code{what == "final"} summary returns a list of 3:
\item{summary.gam}{containing the values of \code{\link{summary.gam}}.}

\item{summary.validation$cv}{cross validation statistics.}
\item{summary.validation$validation}{validation set statistics.}

For \code{what == "path"} summary returns a list of 13:
\item{response}{name of response.}
\item{family}{family used for \code{geoGAM} fit.}
\item{n.obs}{number of observations used for model fitting.}
\item{n.obs.val}{number of observations used for model validation.}
\item{n.covariates}{number of initial covariates including factors.}
\item{n.cov.chosen}{number of covariates in final model.}
\item{list.factors}{list of factors chosen as offset.}
\item{mstop}{number of optimal iterations of gradient boosting.}
\item{list.baselearners}{list of covariate names selected by gradient boosting.}
\item{list.effect.size}{list of covariate names after cross validation of effect size in gradient boosting.}
\item{list.backward}{list of covariate names after backward selection.}
\item{list.aggregation}{list of aggregated factor levels.}
\item{list.gam.final}{list of covariate names in final model.}
}
\references{
%% ~put references to the literature/web site here ~
Nussbaum, M., Papritz, A., and others (in prep). Mapping of soil properties at high resolution by
geo-additive models.
}
\author{
%%  ~~who you are~~
Madlene Nussbaum, \email{madlene.nussbaum@env.ethz.ch}
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{geoGAM}}, \code{\link{gam}}, \code{\link{predict.gam}}
}
\examples{

### small example with earthquake data

data(quakes)
set.seed(2)

quakes <- quakes[ sample(1:nrow(quakes), 50), ]

quakes.geogam <- geoGAM(response = "mag",
                        covariates = c("depth", "stations"),
                        data = quakes,
                        seed = 2,
                        max.stop = 5,
                        cores = 1)

summary(quakes.geogam)
summary(quakes.geogam, what = "path")

plot(quakes.geogam)
plot(quakes.geogam, what = "path")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }
\keyword{ models & regression & nonlinear }% __ONLY ONE__ keyword per line
