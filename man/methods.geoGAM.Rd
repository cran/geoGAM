\name{methods}
\alias{summary.geoGAM}
\alias{print.geoGAM}
\alias{plot.geoGAM}
\alias{summary}
\alias{print}
\alias{plot}
\title{

Methods for \code{geoGAM} objects
}
\description{

Methods for models fitted by \code{geoGAM()}.
}
\usage{
\method{summary}{geoGAM}(object, \dots, what = c("final", "path"))

\method{print}{geoGAM}(x, \dots)

\method{plot}{geoGAM}(x, \dots, what = c("final", "path"))
}

\arguments{
  \item{object}{an object of class \code{geoGAM}}
  \item{x}{an object of class \code{geoGAM}}
  \item{\dots}{other arguments passed to \code{summary.gam}, \code{plot.gam} or \code{plot.mboost}}
  \item{what}{print summary or plot partial effects of \code{final} selected model or print summary or plot gradient boosting path of model selection \code{path}.}
}

\details{

\code{summary} with \code{what = "final"} calls \code{\link{summary.gam}} to display a summary of the final (geo)additive model. \code{plot} with  \code{what = "final"} calls \code{\link{plot.gam}} to plot partial residual plots of the final model.

\code{summary} with \code{what = "path"} give a summary of covariates selected in each step of model building.
\code{plot} with \code{what = "path"} calls \code{plot.mboost} to plot the path of the gradient boosting algorithm.

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

Nussbaum, M., Walthert, L., Fraefel, M., Greiner, L., and Papritz, A.: Mapping of soil properties at high resolution in Switzerland using boosted geoadditive models, SOIL, 3, 191-210, doi:10.5194/soil-3-191-2017, 2017.
}
\author{

M. Nussbaum
}
\seealso{

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
\keyword{ spatial }
\keyword{ models & regression & nonlinear }
