\name{predict.geoGAM}
\alias{predict}
\alias{predict.default}
\alias{predict.geoGAM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Prediction from fitted geoGAM model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes a fitted \code{\link{geoGAM}} object and produces point predictions for a new set of covariate values. If no new data is provided fitted values are returned. Centering and scaling is applied with the same parameters as for the calibration \code{data} set given to \code{\link{geoGAM}}. Factor levels are aggregated according to the final model fit.
}
\usage{
\method{predict}{geoGAM}(object, newdata,
        type = c("response", "link", "probs", "class"),
        back.transform = c("none", "log", "sqrt"),
        threshold = 0.5, se.fit = F, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{geoGAM}}
  \item{newdata}{An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used. If newdata is provided then it should contain all the variables needed for prediction: a warning is generated if not.}
  \item{type}{Type of prediction.}
  \item{back.transform}{Should to \code{log} or \code{sqrt} transformed responses unbiased back transformation be applied? Default is \code{none}. Ignored for categorical responses.}
  \item{threshold}{ Ignored for \code{type = c("response", "link", "probs")} and for \code{type = "class"} for responses with more than two levels.}
  \item{se.fit}{logical. Default is FALSE.}
  \item{\dots}{further arguments to \code{predict()}.}

}

\details{
%%  ~~ If necessary, more details than the description above ~~
	Returns point predictions for new locations \eqn{s} from linear and smooth trends \eqn{\hat f(\mathbf{x},s)}{f(x,s)} estimated
	by penalized least squares geoGAM by calling the function \code{\link{predict.gam}}.
	For lognormal responses (\code{back.transform = 'log'}) unbiased
	back transformation (\cite{Cressie 1993},p 135) was computed by

  \deqn{ \tilde Y(s) = exp( \hat f(\mathbf{x}, s) - \frac{1}{2} \hat \sigma^2 ) }{
  Y(s) = exp( f(x,s) - 1/2*sigma^2) }

	with \eqn{\hat\sigma^2}{sigma^2} being the estimated variance of \eqn{\hat
	f(\mathbf{x}, s)}{f(x,s)} (see \code{\link{predict.gam}} with \code{se.fit=TRUE}). For responses with square root transformation (\code{back.transform = 'sqrt'}) unbiased backtransform was computed by

	\deqn{
	\tilde Y(s) = \hat f(\mathbf{x}, s)^2 - \hat \sigma^2 }{
  Y(s) = f(x,s)^2 + sigma^2	}

  .\cr
	For binary and ordered responses predictions yield
	predicted occurrence probabilities \eqn{\tilde
	P(Y(\mathbf{s})=\mathbf{r}|\mathbf{x},s)}{P(Y(s)=r|x,s)} for response classes \eqn{\mathbf{r}}{r}.

	To obtain binary class predictions a \code{threshold} can be given. A threshold of 0.5 (default) maximizes percentage correct of predicted classes. For binary responses of rare events this threshold may not be optimal. Maximizing on e.g. Gilbert Skill Score (GSS, Wilks, 2011, chap. 8) on cross-validation predictions of the final geoGAM might be a better strategy. GSS is excluding the correct predictions of the more abundant class 	and is preferably used in case of unequal distribution of binary responses
	(direct implementation of such a cross validation procedure planed.)


	For ordered responses \code{predict} with \code{type = 'class'} selects the class to which the median of the
	probability distribution over the ordered categories is assigned (\cite{Tutz 2012, p. 475}).

}
\value{

Vector of point predictions for the sites in \code{newdata} is returned, with unbiased back transformation applied according to option \code{back.transform}.

If \code{se.fit = TRUE} then a 2 item list is returned with items \code{fit} and \code{se.fit} containing predictions and associated standard error estimates as computed by \code{\link{predict.gam}}.

}
\references{
%% ~put references to the literature/web site here ~
Cressie, N. A. C. 1993. Statistics for Spatial Data, John Wiley \& Sons.

Nussbaum, M., Papritz, A., and others (in prep). Mapping of soil properties at high resolution by
geo-additive models.

Tutz, G. 2012. Regression for Categorical Data, Cambridge University Press.

Wilks, D. S. 2011. Statistical Methods in the Atmospheric Sciences, Academic Press.
}
\author{
%%  ~~who you are~~
Madlene Nussbaum, \email{madlene.nussbaum@env.ethz.ch}
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{geoGAM}}, \code{\link{gam}}, \code{\link{predict.gam}}, \code{\link{summary.geoGAM}}, \code{\link{plot.geoGAM}}
}
\examples{

data(quakes)
set.seed(2)

quakes <- quakes[ ss <- sample(1:nrow(quakes), 50), ]

# Artificially split data to create prediction data set
quakes.pred <- quakes[ -ss, ]

quakes.geogam <- geoGAM(response = "mag",
                        covariates = c("depth", "stations"),
                        data = quakes,
                        max.stop = 5,
                        cores = 1)

predicted <- predict(quakes.geogam, newdata = quakes.pred, type = "response" )


\dontrun{


## Use soil data set of soil mapping study area near Berne

library(raster)

data(berne)
data(berne.grid)

# Split data sets and
# remove rows with missing values in response and covariates

d.cal <- berne[ berne$dataset == "calibration" & complete.cases(berne), ]

### Model selection for binary response
ph10.geogam <- geoGAM(response = "ph.0.10",
                             covariates = names(d.cal)[14:ncol(d.cal)],
                             coords = c("x", "y"),
                             data = d.cal,
                             seed = 1)

# Create GRID output with predictions
sp.grid <- berne.grid[, c("x", "y")]

sp.grid$pred.ph.0.10 <- predict(ph10.geogam, newdata = berne.grid)

# transform to sp object
coordinates(sp.grid) <- ~ x + y

# assign Swiss CH1903 / LV03 projection
proj4string(sp.grid) <- CRS("+init=epsg:21781")

# transform to grid
gridded(sp.grid) <- TRUE

plot(sp.grid)

# optionally save result to GeoTiff
# writeRaster(raster(sp.grid, layer = "pred.ph.0.10"),
#             filename= "raspH10.tif", datatype = "FLT4S", format ="GTiff")

}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }
\keyword{ models & regression & nonlinear }% __ONLY ONE__ keyword per line
