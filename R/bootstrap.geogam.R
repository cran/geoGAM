


bootstrap <- function(object, ...) { UseMethod("bootstrap") }

bootstrap.default <- function(object, ...) { UseMethod("bootstrap") }

### 7. Prediction Errors by Bootstrap Methods --------------

# Bootstrapped prediction errors for the full GAM model selection process
bootstrap.geoGAM <- function( object,
                              newdata,
                              R = 100,
                              back.transform = c("none", "log", "sqrt"),
                              seed = NULL,
                              cores = detectCores(), ...)  {


  if( object$parameters$family[[2]] != "Gaussian") stop( "No bootstrap available for categoric response.")

  back.transform <- match.arg(back.transform)


  # if windows set cores to 1
  # if ( get_os() == "windows") n.cores <- 1


  # simulate new response
  f.simulate.response <- function(jj, object, t.sigma){

    # 0) Take predictions from Boosting with Original Data
    l.fit <- fitted( object$gam.final )

    # 0) Standardize residuals according to Hinkely
    mod.resid <- object$gam.final$residuals / sqrt( 1 - object$gam.final$hat )
    mod.resid <- mod.resid - mean(mod.resid)

    # 1) Simulate new Observations, with f(x)_n, page 262, Davison + Hinkley
    return( l.fit + sample(mod.resid, length(l.fit), replace = TRUE) )

  }

  # model fit
  f.gam.model.fit.on.err <- function(jj, object, sim.resp, pred.original){

    cat( "Bootstrap iteration: ", jj, "\n\n")

    # 0) Take predictions from Boosting with Original Data
    d.boot <- object$data
    l.fit <- fitted( object$gam.final )

    # 1) Simulate new Observations, with f(x)_n
    d.boot[, object$parameters$response] <- sim.resp[[jj]] #rnorm( l.fit*0, mean = l.fit, sd = t.sigma )

    # 2) Fit Model to new Observations with full model selection process
    t <- capture.output( boot.fit <- geoGAM( response = object$parameters$response,
                                             covariates = object$parameters$covariates,
                                             sets = object$parameters$sets,
                                             offset = object$parameters$off.choice,
                                             coords = object$parameters$coords,
                                             data = d.boot,
                                             weights = object$parameters$weights,
                                             non.stationary = object$parameters$non.stationary,
                                             max.stop = object$parameters$max.stop,
                                             verbose = 0, cores = 1) )

    # 3) Prediction
    if( back.transform == "none"){

      new.dat.pred <- predict(boot.fit, newdata = newdata)

    } else {

      pred.obj <- predict(boot.fit, newdata = newdata, se.fit = TRUE, back.transform = back.transform)
      new.dat.pred <- as.numeric(pred.obj$pred)
      new.dat.pred.se <- as.numeric(pred.obj$pred.se)

    }


    # Compute prediction error, according to Hinkley, p. 285, Alg. 6.4
    # Notizen im Heft 2013-2014 (klein blau)
    #
    # Compute sigma
    df <- round( sum( summary( boot.fit$gam.final )$edf ) +
                   sum( summary(boot.fit$gam.final)$pTerms.df) )
    boot.res <- resid( boot.fit$gam.final )
    boot.sigma <- sqrt( sum( ( boot.res - mean( boot.res ) )^2 ) * 1/(length(boot.res)-df) )

    # Compute prediction + prediction error
    p.err <- new.dat.pred - pred.original +
      rnorm(new.dat.pred*0, mean = 0, sd = boot.sigma)

    new.dat.perr <- new.dat.pred - p.err # Hinkley, p. 286

    # unbiased backtransform
    if( back.transform == "sqrt"){
      new.dat.perr <- as.numeric( new.dat.perr^2 - new.dat.pred.se  )
    }

    if( back.transform == "log"){
      new.dat.perr <- as.numeric( exp( new.dat.perr - new.dat.pred.se*0.5 ) )
    }


    return( new.dat.perr )
  }


  pred.original <- predict( object, newdata = newdata, type = "response", se.fit = FALSE, back.transform = "none" )


  ## compute residual standard error sigma
  # Degrees of Freedom from GAM model (approximate)
  df <- round( sum( summary(object$gam.final )$edf ) +
                 sum( summary( object$gam.final )$pTerms.df) )
  t.res <- resid( object$gam.final)
  t.sigma <- sqrt( sum( ( t.res - mean( t.res ) )^2 ) * 1/(length(t.res)-df) )


  # Apply to number of bootstraps
  # ensure repeatable random number generation
  #RNGkind( "L'Ecuyer-CMRG" )
  if( !is.null(seed)) set.seed(seed)
  #mc.reset.stream()

  # 1. simulate new responses
  sim.resp <- mclapply(1:R, f.simulate.response, mc.allow.recursive = FALSE,
                      mc.set.seed = TRUE, mc.cores = cores,
                      object = object, t.sigma = t.sigma)

  # 2. fit + predict
  m.p.err <- mclapply(1:R, f.gam.model.fit.on.err, mc.allow.recursive = FALSE,
                      mc.set.seed = FALSE, mc.cores = cores,
                      object = object, sim.resp = sim.resp, pred.original = pred.original)

  d.dist <- matrix( unlist( m.p.err ), byrow = FALSE, ncol = R )

  d.dist  <- data.frame( cbind( newdata[, object$parameters$coords ], d.dist  ) )
  names(d.dist) <- c( object$parameters$coords, paste0("P", 1:R) )


  return(d.dist)

}


