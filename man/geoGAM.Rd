\name{geoGAM}
\alias{geoGAM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Select sparse geoadditive model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Selects a sparse geoadditive model from a large set of covariates with the aim of spatial prediction.
As covariates categoric and continuous are handeled and smooth spatial surfaces can be included for model selection.
}
\usage{
geoGAM(response, covariates = names(data)[!(names(data) \%in\% c(response,coords))],
       data, coords = NULL, weights = rep(1, nrow(data)),
       offset = T, max.stop = 300, non.stationary = F,
       sets = NULL, seed = NULL, validation.data = NULL,
       verbose = 0, cores = min(detectCores(),10))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response}{ name of response as character. Responses currently supported: gaussian, binary, ordered.}
  \item{covariates}{ character vector of all covariates (factor, continuous). If not given, all columns of \code{data} are used.}
  \item{data}{ data frame containing response, coordinates and covariates. }
  \item{coords}{character vector of column names indicating spatial coordinates.}
  \item{weights}{ weights used for model fitting.}
  \item{offset}{ logical, use offset for component wise gradient boosting algorithm.}
  \item{max.stop}{ maximal number of boosting iterations.}
  \item{non.stationary}{logical, include non-stationary effects in model selection. This allows for spatial varying coefficients for continuous covariates, but increases computational effort.}
  \item{sets}{ give predefined cross validation sets. }
  \item{seed}{ set random seed for splitting of the cross validation sets, if no \code{sets} are given.}
  \item{validation.data}{data frame containing response, coordinates and covariates to compute independent validation statistics. This data set is used to calculate predictive performance at the end of model selection only.}
  \item{verbose}{Should screen output be generated? 0 = none, >0 create output.}
  \item{cores}{ number of cores to be used for parallel computing }
}

\details{

\strong{Generic model representation}

	% continuous response
	Response \eqn{\mathit{Y}}{Y} at location \eqn{\mathit{s}}{s} measured
	on interval scale is modeled by the regression function
	\eqn{f} from the environmental covariates \eqn{\mathbf{x}}{x} as

  \deqn{
	Y(\mathit{\mathbf{x},s})=f(\mathbf{x},s) + \epsilon}{
	Y(x,s) = f(x,s) + epsilon
	}

	with \eqn{\epsilon}{epsilon} being a zero mean spatially uncorrelated error that
	follows a given distribution (e.g. Gaussian).
	For positively skewed \eqn{Y(\mathit{\mathbf{x},s})}{Y(x,s)} transformed by natural
	logarithm \eqn{\epsilon}{epsilon} follows a lognormal distribution.
	%
	% Binary response
	The probability of presence of a event or an observed class \eqn{r} is modeled
	similarly using the inverse logistic transformation (e.g. \cite{Tutz 2012}, p. 37) by

  \deqn{
	P(Y(\mathbf{s}) = r| \mathbf{x},s) = logit^{-1} \Big(~f( \mathbf{x},s)~
	\Big) =
	\frac{ exp{(~f( \mathbf{x},s)~)}  }{  1 + exp{(~f( \mathbf{x},s)~)}}.
  }{
  P(Y(s)=r|x,s) = logit^-1 f(x,s) = exp(f(x,s)) / ( 1 + exp(f(x,s)) )
  }

	Analogous ordered classes at location \eqn{\mathit{s}}{s} are modeled
	by cumulative probabilities (\cite{Tutz-2012}, pp. 244, cumulative logit model, proportional odds model)
	of the class being equal or smaller than \eqn{r} with

  \deqn{
	P(Y(\mathbf{s}) \leq r| \mathbf{x},s) = \frac{ exp{(~f(\nu_r,
			\mathbf{x},s)~)}  }{  1 + exp{(~f(\nu_r, \mathbf{x},s)~)}},
  }{
  P(Y(s)>=r|x,s) = exp(f(nu_r,x,s)) / ( 1 + exp(f(nu_r,x,s)))
  }

	where \eqn{\nu_1 ... \nu_{r-1}}{nu_1 ... nu_(r-1)} are class specific thresholds or model
	intercepts. This specification follows the proportional odds
	assumption by estimating equal covariate effects independent from the
	target category. \\


  .\cr
\strong{Model selection}

Sparse models are desirable if models should be open to feasibility checks in regard of meaningful relationships between response and (environmental) covariates. A fully automated model selection procedure is required for large number of responses or simulation of predictive distributions (\code{\link{bootstrap.geoGAM}}). The following steps are performed in \code{\link{geoGAM}}:

\enumerate{

		% centering
		\item To ensure stability of the algorithm continuous covariates are
		centered and scaled by their range.

		% Offset Lasso, polr
		\item (optional with \code{offset = TRUE}) The group Lasso (least absolute shrinkage and selection
		operator, \cite{Breheny and Huang 2015}, \code{\link{grpreg}})  -- an algorithm that likely excludes
		non-relevant covariates and treats categorical covariates as groups --
		is used to select relevant categorical covariates.
		%
		For ordered responses
		step wise proportional odds logistic regression (\code{\link{polr}}) in both directions with
		BIC is used, because no Lasso is available.

		% Boosting, baseleaner definition
		\item By component wise gradient boosting a subset of
		relevant continuous and categoric covariates and spatial effects is selected.
		Boosting is a slow stage-wise additive learning algorithm. It expands
		\eqn{f(\mathbf{x},s)}{f(x,s)} in a set of base procedures (baselearners) and
		approximates \eqn{f( \mathbf{x},s)}{f(x,s)} by a finite sum of them. Best fitting
		baselearners update the model only in small steps.
		% step size
		The default step length (\eqn{\nu = 0.1}{v = 0.1}) is used as this is not a
		sensitive parameter as long as it is small (\cite{Hofner et al. 2012}).
		Early stopping of the algorithm (\code{mstop}) is achieved by
		minimizing cross validation error.
		For continuous covariates non-parametric penalized smoothing
		spline baselearners (\cite{Kneib et al. 2009}) are used. Categorical covariates are
		assigned to linear baselearners. To capture the spatial structure in
		the data smooth spatial surface established by a bivariate
		tensor-product P-spline of spatial coordinates are added.
		Spatially varying effects are estimated by baselearners formed by the
		product continuous covariates and the smooth spatial surface (\code{non.stationary = TRUE}).
		%
		To equalize the inclusion probability of baselearners each should
		have the same flexibility to fit the data
		(\cite{Hofner et al. 2011}). Unbiased baselearner
		selection is controlled by assigning 5 degrees of freedom to each.
		Where linear categorical covariates did not reach 6 levels (5 degrees of freedom) we
		aggregated them to grouped baselearners.

		% offset
		With option \code{offset = TRUE} fitted values of group Lasso and proportional odds regression
		respectively (step 1) are used as an offset to the boosting
		algorithm. Setting an offset accelerates model building and ensures
		inclusion of relevant categorical covariates due to the limit of
		equal degrees of freedom of 5.


		% Select relevant baselarners by magnitude
		\item To remove baselearners with strongly shrunken coefficients the effect size is evaluated.
		As effect size the range of coefficients for categorical covariates is computed and
		the magnitude of smooth baselearners is obtained after removal of extreme
		values from the partial effects (outlier < 1. quantile -
		1.5*interquantile distance and > 3. quantile + 1.5* interquantile
		distance, \cite{Frigge et al. 1989}, \code{\link{boxplot}}).

		The optimal effect size is found by fitting geoadditive models by
		retaining the same degrees of freedom as in gradient boosting for the
		non-parametric smooth terms. The categorical covariates used to compute
		the offset (step 1) are included as linear parametric effects.

		% Backward selection
		\item The model is further reduced by stepwise removal of covariates.
		The covariate to drop is chosen by largest \code{p} value of \code{F} tests for categorical terms and approximate \code{F} test
		(\cite{Wood 2011}) for smooth terms.

		% Factor aggregation
		\item Similar levels of categorical covariates are merged stepwise
		based on highest \code{p} values from two sample \code{t} tests (\code{\link{t.test}}) of partial
		residuals. For interactions by smooth terms and categorical covariates
		the t-test of the interaction coefficient obtained from linear
		regression computed on the partial residuals are used.

}

	The optimal model building parameters (\eqn{\lambda}{lambda} in group Lasso,
	number of boosting iterations (\eqn{\nu}{v}), magnitude of baselearners, number of
	remaining covariates and degree of aggregation of factor levels) is done by
	minimizing root mean squared error (RMSE) for interval scaled, Brier Score
	for binary and ranked probability score (\cite{Wilks 2011, chap. 8}) for ordered
	responses computed from 10fold cross-validation done with the same subsets.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use

Object of class \code{geoGAM}:

\item{offset.grplasso}{Cross validation for grouped LASSO, object of class \code{\link{cv.grpreg}} of package \code{\link{grpreg}}). Empty for \code{offset = FALSE}.}
\item{offset.factors}{Character vector of factor names chosen for the offset computation. Empty for \code{offset = FALSE}.}
\item{gamboost}{Gradient boosting with smooth components, object of class \code{\link{gamboost}} of package \code{\link{mboost}}.}
\item{gamboost.cv}{Cross validation for gradient boosting, object of class \code{\link{cvrisk}} of package \code{\link{mboost}}.}
\item{gamboost.mstop}{Mstop used for gamboost.}
\item{gamback.cv}{List of cross validation error for tuning parameter magnitude.}
\item{gamback.backward}{List of cross validation error path for backward selection of \code{\link{gam}} fit.}
\item{gamback.aggregation}{List(s) of cross validation error path for aggregation of factor levels.}
\item{gam.final}{Final selected geoadditive model fit, object of class \code{\link{gam}}.}
\item{gam.final.cv}{Data frame with original response and cross validation predictions.}
\item{gam.final.extern}{Data frame with original response data and predictions of \code{gam.final}.}
\item{data}{Original data frame for model calibration.}
\item{parameters}{List of parameters handed to geoGAM (used for subsequent bootstrap of prediction intervals).}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Breheny, P. and Huang, J., 2015. Group descent algorithms for nonconvex penalized linear
and logistic regression models with grouped predictors. Statistics and Computing, 25, 173--187.

Frigge, M., Hoaglin, D. C., and Iglewicz, B., 1989. Some implementations of the boxplot.
The American Statistician, 43(1), 50--54.

Hofner, B., Hothorn, T., Kneib, T., and Schmid, M., 2011. A framework for unbiased
model selection based on boosting. Journal of Computational and Graphical Statistics,
20(4), 956--971.

Hofner, B., Mayr, A., Robinzonov, N., and Schmid, M., 2012. Model-based boosting in R,
a hands-on tutorial using the R package mboost. Technical report no. 120,
Ludwig-Maximilians-Universitaet Muenchen, Department of Statistics.
\url{http://epub.ub.uni-muenchen.de/12754/}

Kneib, T., Hothorn, T., and Tutz, G., 2009. Variable selection and model choice in geoad-
ditive regression models. Biometrics, 65(2), 626--634.

Nussbaum, M., Papritz, A., and others, in prep. Mapping of soil properties at high resolution by
geo-additive models.

Wood, S. N., 2011. Fast stable restricted maximum likelihood and marginal likelihood
estimation of semiparametric generalized linear models. Journal of the Royal Statistical
Society (B), 73(1), 3--36.
}
\author{
%%  ~~who you are~~
Madlene Nussbaum, \email{madlene.nussbaum@env.ethz.ch}
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
The model selection is based on packages \code{\link{grpreg}} (function \code{\link{cv.grpreg}}), \code{MASS} (function \code{\link{polr}}), \code{\link{mboost}} (functions \code{\link{gamboost}}, \code{\link{cv}}, \code{\link{cvrisk}}) and \code{\link{mgcv}} (function \code{\link{gam}}). For further information please see documentation and vignettes for these packages.
}
\examples{
### small examples with earthquake data

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

\donttest{
data(quakes)

# create grouped factor with reduced number of levels
quakes$stations <- factor( cut( quakes$stations, breaks = c(0,15,19,23,30,39,132)) )

quakes.geogam <- geoGAM(response = "mag",
                        covariates = c("stations", "depth"),
                        coords = c("lat", "long"),
                        data = quakes,
                        max.stop = 10,
                        cores = 1)

summary(quakes.geogam)
summary(quakes.geogam, what = "path")
}

\dontrun{

## Use soil data set of soil mapping study area near Berne

data(berne)
set.seed(1)

# Split data sets and
# remove rows with missing values in response and covariates

d.cal <- berne[ berne$dataset == "calibration" & complete.cases(berne), ]
d.val <- berne[ berne$dataset == "validation" & complete.cases(berne), ]


### Model selection for continuous response
ph10.geogam <- geoGAM(response = "ph.0.10",
                      covariates = names(d.cal)[14:ncol(d.cal)],
                      coords = c("x", "y"),
                      data = d.cal,
                      offset = T,
                      sets = mboost::cv(rep(1, nrow(d.cal)), type = "kfold"),
                      validation.data = d.val)
summary(ph10.geogam)
summary(ph10.geogam, what = "path")


### Model selection for binary response
waterlog100.geogam <- geoGAM(response = "waterlog.100",
                             covariates = names(d.cal)[c(14:54, 56:ncol(d.cal))],
                             coords = c("x", "y"),
                             data = d.cal,
                             offset = F,
                             sets = sample( cut(seq(1,nrow(d.cal)),breaks=10,labels=FALSE) ),
                             validation.data = d.val)
summary(waterlog100.geogam)
summary(waterlog100.geogam, what = "path")


### Model selection for ordered response
dclass.geogam <- geoGAM(response = "dclass",
                        covariates = names(d.cal)[14:ncol(d.cal)],
                        coords = c("x", "y"),
                        data = d.cal,
                        offset = T,
                        non.stationary = T,
                        seed = 1,
                        validation.data = d.val)
summary(dclass.geogam)
summary(dclass.geogam, what = "path")

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }
\keyword{ models & regression & nonlinear }% __ONLY ONE__ keyword per line
