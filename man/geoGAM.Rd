\name{geoGAM}
\alias{geoGAM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Select sparse geoadditive model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Selects a parsimonious geoadditive model from a large set of covariates with the aim of (spatial) prediction.
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

\strong{Summary}

\code{\link{geoGAM}} models smooth nonlinear relations between responses and single covariates and combines these model terms additively. Residual spatial autocorrelation is captured by a smooth function of spatial coordinates and nonstationary effects are included by interactions between covariates and smooth spatial functions. The core of fully automated model building for geoGAM is componentwise gradient boosting. The model selection procedures aims at obtaining sparse models that are open to check feasibilty of modelled relationships (\cite{Nussbaum et al. 2017a}).

\code{\link{geoGAM}} to date models continuous, binary and ordinal responses. It is able to cope with numerous continuous and categorical covariates.

\strong{Generic model representation}


  GAM expand the (possibly transformed) conditional
	expectation of a response at given covariates \eqn{s} as an additive series
  %
 \deqn{
    g\left(\rule{0pt}{14pt}\mathrm{E}[Y(\ensuremath{\mathbf{s}})\,|\,\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})]\right)
        = \nu + f(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}))
    = \nu +  \sum_{j} f_{j}(x_{j}(\ensuremath{\mathbf{s}})),
  }{
  g(E[Y(s)|x(s)]) = nu + f(x(s)) = nu + sum( f_j(x_j(s)) ),
  }
  %
	where \eqn{\nu}{nu} is a constant and \eqn{f_{j}(x_{j}(\ensuremath{\mathbf{s}}))}{f_j(x_j(s))} are linear
	terms or unspecified ``smooth'' nonlinear functions of single covariates \eqn{x_{j}(\ensuremath{\mathbf{s}})}{x_j(s)}
	(e.g. smoothing spline, kernel or
	any other scatterplot smoother) and \eqn{g(\cdot)}{g(.)} is again a link function. A generalized additive model (GAM) is based on the following components (\cite{Hastie and Tibshirani 1990, Chapt. 6}):

\enumerate{
  %
  \item \emph{Response distribution}: Given \eqn{\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}) = x_1(\ensuremath{\mathbf{s}}), x_2(\ensuremath{\mathbf{s}}), \\ ..., x_p(\ensuremath{\mathbf{s}})}{x(s)=x_1(s), x_2(s), ..., x_p(s)}, the
  \eqn{Y(\ensuremath{\mathbf{s}})}{Y(s)} are conditionally independent observations from simple
  exponential family distributions.
  %
  \item \emph{Link function}: \eqn{g(\cdot)}{g(.)} relates the expectation
  \eqn{\mu(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})) = \mathrm{E}[Y(\ensuremath{\mathbf{s}})|\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})]}{mu(x) = E[Y(s)|x(s)]} of the response
  distribution to
  %
  \item the \emph{additive predictor}  \eqn{\sum_{j} f_{j}(x_{j}(\ensuremath{\mathbf{s}}))}{sum( f_j(x_j(s)) )}.
}
  geoGAM extend GAM by allowing a more complex form of the additive
  predictor (\cite{Kneib et al. 2009, Hothorn et al. 2011}): First, one can
  add a smooth function \eqn{f_{\ensuremath{\scriptstyle \mathbf{s}}}(\ensuremath{\mathbf{s}})}{f_s(s)} of the spatial
  coordinates (smooth spatial surface) to the additive predictor to account for residual
  autocorrelation.
 	%
  More complex relationships between \eqn{Y} and \eqn{\ensuremath{\mathbf{x}}}{x} can be modelled
  by adding terms like \eqn{f_{j}(x_{j}(\ensuremath{\mathbf{s}})) \cdot
  f_{k}(x_{k}(\ensuremath{\mathbf{s}}))}{f_j(x_j(s))*f_k(x_k(s))} -- capturing the effect of interactions
  between covariates -- and \eqn{f_{\ensuremath{\scriptstyle \mathbf{s}}}(\ensuremath{\mathbf{s}}) \cdot
  f_{j}(x_{k}(\ensuremath{\mathbf{s}}))}{f_s(s)*f_j(x_k(s))} -- accounting for spatially changing
  dependence between \eqn{Y} and \eqn{\ensuremath{\mathbf{x}}}{x}.  Hence, in its full generality,
  a generalized additive model for spatial data is represented by
  %
  \deqn{
    g(\mu(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}))) = \nu + f(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}))   =  \nonumber }{}
  \deqn{
    \nu +
    \underbrace{
    \sum_{u} f_{j_{u}}(x_{j_{u}}(\ensuremath{\mathbf{s}})) + \sum_{v}
    f_{j_{v}}(x_{j_{v}}(\ensuremath{\mathbf{s}})) \cdot f_{k_{v}}(x_{k_{v}}(\ensuremath{\mathbf{s}}))
    }_{\mbox{global marginal and interaction effects}} \nonumber }{}
  \deqn{
     +
    \underbrace{ \sum_{w} f_{\ensuremath{\scriptstyle \mathbf{{s}}_{w}}}(\ensuremath{\mathbf{s}}) \cdot
    f_{j_{w}}(x_{j_{w}}(\ensuremath{\mathbf{s}})) }_{\mbox{nonstationary effects}} +
    \underbrace{\hspace{5mm} f_{\ensuremath{\scriptstyle \mathbf{s} }}(\ensuremath{\mathbf{s}})
    \hspace{5mm}}_{\mbox{autocorrelation}}.
    }{
    g(mu(x(s))) = nu + f(x(s)) = nu + sum( f_j(x_j(s)) ) + sum( f_j(x_j(s))*f_k(x_k(s)) ) + sum( f_s(s)*f_j(x_j(s)) ) + f_s(s).
    }
  %
  \cite{Kneib et al. (2009)} called the above equation a geoadditive model,
  a name coined before by \cite{Kammann and Wand 2003} for a combination
  of additive models with a geostatistical error model.
  It remains to specify what response distributions and link functions
  should be used for the various response types: For (possibly
  transformed) \emph{continuous} responses one uses often a normal
  response distribution combined with the identity link
  \eqn{g\left(\mu(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}))\right) = \mu(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}))}{
  g(mu(x(s))) = mu(x(s))}.
  For binary data (coded as 0 and 1), one assumes a Bernoulli
  distribution and uses often a logit link
  %
  \deqn{
     g\left(\mu(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}))\right) =\log\left(
     \frac{\mu(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}))}{1-\mu(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}))} \right),
  }{
  g(mu(x(s))) =log( mu(x(s)) / ( 1-mu(x(s)) ) )
  }
  %
  where
  %
  \deqn{
    \mu(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})) =
    \mathrm{Prob}[Y(\ensuremath{\mathbf{s}})=1\,|\,\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})] =
    \frac{\exp(\nu +f(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})))}{1+\exp(\nu +f(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})))}.
  }{
  mu(x(s)) = Prob[Y(s) = 1|x(s)] = exp(nu + f (x(s))) / (1 + exp(nu + f (x(s))))
  }
  %
  For ordinal data, with ordered response levels, \eqn{1, 2, \ldots, k}{1, 2, ..., k},
  the cumulative logit or proportional odds model
  (\cite{Tutz 2012}, Sect. 9.1) is used.  For any given level \eqn{r \in (1, 2,
  \ldots, k)}{r in(1, 2, ..., k)}, the logarithm of the odds of the event \eqn{Y(\ensuremath{\mathbf{s}}) \leq
  r \, | \, \ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})}{Y(s) <= r|x(s)} is then modelled by
  %
  \deqn{
    \log\left(
    \frac{\mathrm{Prob}[Y(\ensuremath{\mathbf{s}}) \leq
    r \, | \, \ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}))]}{\mathrm{Prob}[Y(\ensuremath{\mathbf{s}}) > r \, | \,
    \ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}))]}\right) = \nu_{r} + f(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})),
  }{
  log( Prob[Y(s) <= r|x(s)] / Prob[Y(s) > r|x(s)] ) = nu_r + f(x(s))
  }
  with \eqn{\nu_{r}}{nu_r} a sequence of level-specific constants satisfying
  \eqn{\nu_{1} \leq \nu_{2} \leq \ldots \leq \nu_{r}}{nu_1 <= nu_2 <= ... <= nu_r}. Conversely,
  %
  \deqn{
    \mathrm{Prob}[Y(\ensuremath{\mathbf{s}})\leq r\,|\,\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})] =
    \frac{\exp(\nu_{r} + f(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})))}{1+\exp(\nu_{r} + f(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})))}. }{
    Prob[Y(s) <= r|x(s)] = exp(nu_r + f(x(s))) / (1 + exp(nu_r + f(x(s)))).
    }
  %
  Note that \eqn{\mathrm{Prob}[Y(\ensuremath{\mathbf{s}})\leq r\,|\,\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})]}{Prob[Y(s) <= r|x(s)]}
  depends on \eqn{r} only through the constant \eqn{\nu_{r}}{nu_r}.  Hence, the ratio
  of the odds of two events \eqn{Y(\ensuremath{\mathbf{s}}) \leq r \, | \,
  \ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}})}{Y(s) <= r|x(s)} and \eqn{(\ensuremath{\mathbf{s}}) \leq r \, | \,
  \tilde{\ensuremath{\mathbf{x}}}(\ensuremath{\mathbf{s}})}{Y(s) <= r|x(s)} is the same for all \eqn{r}
  (\cite{Tutz 2012}, p. 245).




\strong{Model building (selection of covariates)}


	To build parsimonious models that can readily be checked for agreement understanding in regards to the analized subject. The following steps 1--6 are implemented in \code{\link{geoGAM}} toa achieve sparse models in a fully automated way.
	In several of these steps tuning parameters are optimized by 10-fold cross-validation with fixed subsets using either root mean squared error (RMSE), continuous responses), Brier score (BS), binary responses) or ranked probability score (RPS), ordinal responses) as optimization criteria (see \cite{Wilks, 2011}).
	To improve the stability of the algorithm continuous covariates are first scaled (by difference of maximum and minimum value) and centred.

	\enumerate{

		% Offset Lasso, polr
		\item Boosting (see step 2 below) is more stable and converges more quickly when the effects of categorical covariates (factors) are accounted for as model offset. Therefore, the group lasso (least absolute shrinkage and selection operator, \cite{Breheny and Huang 2015}, \code{\link{grpreg}}))  -- an algorithm that likely
		excludes non-relevant covariates and treats factors
		as groups -- is used to select important factors for the offset. %
		For ordinal responses
		stepwise proportional odds logistic regression in both directions
		with BIC (e. g. \cite{Faraway 2005}, p. 126) is used to select the offset covariates because lasso cannot be
		used for such responses.

		% Boosting, baseleaner definition
		\item Next, a subset of
		relevant factors, continuous covariates and spatial effects is selected by componentwise gradient boosting.
		Boosting is a slow stagewise additive learning algorithm. It
		expands \eqn{f(\ensuremath{\mathbf{x}}(\ensuremath{\mathbf{s}}))}{f(x(s))} in a set of base procedures (baselearners)
		and approximates the additive predictor by a finite sum of
		them as follows (\cite{Buehlmann and Hothorn 2007}):

		\enumerate{
			\item Initialize \eqn{\hat{f}( \mathbf{x}(\ensuremath{\mathbf{s}}))^{[m]}}{f(x(s))^[m]}
			with offset of step 1 above and set \eqn{m=0}.
			\item Increase \eqn{m} by 1. Compute
			the negative gradient vector \eqn{\ensuremath{\mathbf{U}}^{[m]}}{U^[m]} (e.g. residuals) for a loss
			function \eqn{l(\cdot)}{l(.)}.
			\item Fit all baselearners \eqn{g(
			\mathbf{x}(\ensuremath{\mathbf{s}}))_{1..p}}{g(x(s))_(1..p)} to \eqn{\ensuremath{\mathbf{U}}^{[m]}}{U^[m]} and select the baselearner, say
			\eqn{g(\mathbf{x}(\ensuremath{\mathbf{s}}))_{j}^{[m]}}{g(x(s))_j^[m]} that minimizes \eqn{l(\cdot)}{l(.)}.
			\item Update
			\eqn{\hat{f}( \mathbf{x}(\ensuremath{\mathbf{s}}))^{[m]} = \hat{f}( \mathbf{x}(\ensuremath{\mathbf{s}}))^{[m-1]} +
			v\cdot g( \mathbf{x}(\ensuremath{\mathbf{s}}))_{j}^{[m]}}{f(x(s))^[m] = f(x(s))^[m-1] + v*g(x(s))_j^[m]} with step size \eqn{v\leq1}{v<=1}.
			\item Iterate steps (b) to (d) until \eqn{m = m_{stop}}{m = m_stop} (main tuning
			parameter).
		}

		%
		The following settings are used in above algorithm: % loss
		As loss functions \eqn{l(\cdot)}{l(.)} \eqn{L_2} is used for continuous,
		negative binomial likelihood for binary
		(\cite{Buehlmann and Hothorn 2007}) and proportional odds likelihood
		for ordinal responses (\cite{Schmid et al. 2011}). %
		% mstop
		Early stopping of the boosting algorithm is achieved by
		determining optimal \eqn{m_{stop}}{m_stop} by cross-validation. % step size
		Default step length (\eqn{\upsilon = 0.1}{upsilon = 0.1}) is used. This is not a
		sensitive parameter as long as it is clearly below 1 (\cite{Hofner et al. 2014}).
		%
		For continuous covariates penalized smoothing
		spline baselearners (\cite{Kneib et al. 2009}) are used. Factors
		are treated as linear baselearners. To capture residual autocorrelation
		a bivariate tensor-product P-spline of spatial coordinates
		(\cite{Wood 2006, pp. 162}) is added to the additive predictor. Spatially varying effects
		are modelled by baselearners formed by multiplication of
		continuous covariates with tensor-product P-splines of spatial coordinates
		(\cite{Wood 2006, pp. 168}). Uneven degree of freedom of baselearners biases
		baselearner selection (\cite{Hofner et al. 2011b}). Therefore, each baselearner is penalized to 5 degrees of
		freedom (\eqn{df}). Factors with less than 6 levels (\eqn{df<5})
		are aggregated to grouped baselearners. By using an offset, effects of important factors with more than 6 levels
		are implicitly accounted for without penalization.

		% Select relevant baselarners by magnitude
		\item At \eqn{m_{stop}}{m_stop} (see step 2 above), many included baselearners may have very small effects only. To remove these
		the effect size \eqn{e_j} of each baselearner \eqn{f_j(x_j(\ensuremath{\mathbf{s}}))}{f_j(x_j(s)))} is computed. 	For factors the effect size \eqn{e_j} is the largest difference between effects of two levels and for continuous covariates it is equal to the maximum contrast of estimated partial effects (after removal of extreme values as in boxplots, \cite{Frigge et al. 1989}). Generalized additive models (GAM, \cite{Wood 2011}) are fitted including smooth and factor effects depending on the effect size \eqn{e_j} of the corresponding baselearner \eqn{j}. The procedure iterates through \eqn{e_j} and excludes covariates with \eqn{e_j} smaller than a threshold effect size \eqn{e_t}. Optimal \eqn{e_t} is determined by 10-fold cross-validation of GAM. In these GAM fits smooth effects are penalized to 5 degrees of freedom as imposed by componentwise gradient boosting (step 2 above). The factors selected as offset in step 1 are included in the main GAM, that is now fitted without offset.

		% Backward selection
		\item The GAM is further reduced by stepwise removal of
		covariates by cross-validation. The candidate covariate to drop is chosen by largest \eqn{p} value
		of \eqn{F} tests for linear factors and approximate \eqn{F} test
		(\cite{Wood 2011}) for smooth terms. %	For smooth + spatial
		% approximate F-Test (outuput from anova.gam,
		%for single model, drop1)

		% Factor aggregation
		\item Factor levels with similar estimated effects are merged stepwise again by cross-validation
		based on largest \eqn{p} values from two sample \eqn{t}-tests of partial
		residuals.

		% finales modell
		\item The final model (used to compute spatial predictions) results ideally in a parsimonious
		GAM. Because of step 5, factors have possibly a reduced number of coefficients.
		Effects of continuous covariates are modelled by smooth functions and  -- if at all present --
		spatially structured residual variation (autocorrelation) is represented by a smooth spatial surface.
		To avoid over-fitting both types of smooth effects are penalized to 5 degrees of freedom (as imposed by step 2).

}
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
Breheny, P. and Huang, J., 2015. Group descent algorithms for nonconvex penalized linear and logistic regression models with grouped predictors. Statistics and Computing, 25, 173--187.

Buehlmann, P. and Hothorn, T., 2007. Boosting algorithms: Regularization, prediction and model fitting, Stat Sci, 22, 477--505, doi:10.1214/07-sts242.

Faraway, J. J., 2005. Linear Models with R, vol. 63 of Texts in Statistical Science, Chapman & Hall/CRC, Boca Raton.

Frigge, M., Hoaglin, D. C., and Iglewicz, B., 1989. Some implementations of the boxplot. The American Statistician, 43(1), 50--54.

Hastie, T. J. and Tibshirani, R. J., 1990. Generalized Additive Models, vol. 43 of Monographs on Statistics and Applied Probability, Chapman and Hall, London.

Hofner, B., Hothorn, T., Kneib, T., and Schmid, M., 2011. A framework for unbiased model selection based on boosting. Journal of Computational and Graphical Statistics, 20(4), 956--971.

Hofner, B., Mayr, A., Robinzonov, N., and Schmid, M., 2014. Model-based boosting in R: A hands-on tutorial using the R package mboost, Computation Stat, 29, 3--35, doi:10.1007/s00180-012-0382-5.

Hothorn, T., Mueller, J., Schroder, B., Kneib, T., and Brandl, R., 2011. Decomposing environmental, spatial, and spatiotemporal components of species distributions, Ecol Monogr, 81, 329--347.

Kneib, T., Hothorn, T., and Tutz, G., 2009. Variable selection and model choice in geoadditive regression models. Biometrics, 65(2), 626--634.

Nussbaum, M., Walthert, L., Fraefel, M., Greiner, L., and Papritz, A., 2017a. Mapping of soil properties at high resolution in Switzerland using boosted geoadditive models, SOIL Discuss., https://www.soil-discuss.net/soil-2017-13/, in review.

Schmid, M., Hothorn, T., Maloney, K. O., Weller, D. E., and Potapov, S., 2011. Geoadditive regression modeling of stream biological condition, Environ Ecol Stat, 18, 709--733, doi:10.1007/s10651-010-0158-4.

Tutz, G., 2012, Regression for Categorical Data, Cambridge University Press,
doi:10.1017/cbo9780511842061.

Wilks, D. S., 2011. Statistical Methods in the Atmospheric Sciences, Academic Press, 3 edn.

Wood, S. N., 2006. Generalized Additive Models: An Introduction with R, Chapman and Hall/CRC.

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
