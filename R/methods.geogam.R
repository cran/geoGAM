

# Function to create summary of final model
# Redirect to summary.gam on gam object and plot validation statistics

print.geoGAM <- function(x, ...) { UseMethod("summary") }



summary.geoGAM <- function(object, ..., what = c("final", "path")){

  what = match.arg(what)

  switch( what,
          final = {

            ss <-  summary(object$gam.final)

            stat <- ifelse( object$parameters$family[[4]] == "rmse", "st", object$parameters$family[[4]] )
            v.cv <- f.validate.predictions( object$gam.final.cv$data,
                                            object$gam.final.cv[, grepl("pred", names(object$gam.final.cv)) ], stat )


            if( !is.null(object$gam.final.extern) ){

              v.t <- f.validate.predictions(object$gam.final.extern$dat,
                                            object$gam.final.extern[, grepl("pred", names(object$gam.final.extern)) ], stat)


            } else {
              v.t <- NULL}


            l.s <- list( summary.gam = ss, summary.validation = list( cv = v.cv, validation = v.t))
            class(l.s) <- "summary.geoGAM"

            return( l.s )

          },

          path = {

            t.te <- ifelse( sum(grepl( "^te\\(", names(object$gam.final$coefficients))), 1, 0)
            p.lin <- sum( !grepl( "s\\(|te\\(", strsplit( as.character( formula(object$gam.final)[3]), split = "\\+")[[1]]))
            p.smooth <- length(object$gam.final$smooth) - t.te # mit Interaktionen einzeln

            list.aggr <- list()
            if( length(object$gamback.aggregation) > 0){
            for( ll in 1:length(object$gamback.aggregation)){
              list.aggr[[names(object$gamback.aggregation[1])]] <- object$gamback.aggregation[[ll]][[2]]
            }}

            t.final <- names(object$gam.final$model)[
              -match(c(object$parameters$response, "(weights)"), names(object$gam.final$model) )]
            t.f <- unlist(lapply(object$gam.final$model[, t.final, drop = FALSE], is.factor))
            t.final[t.f] <- gsub("ag$", "", t.final[t.f])

            tt <- capture.output(
              t.boost.names <- names(coef(object$gamboost[object$gamboost.mstop] )), type = "message" )

            l.s <- list(

              response = object$parameters$response,
              family = object$parameters$family[[2]],
              n.obs = nrow(object$data),
              n.obs.val = ifelse( is.null(object$validation.data), 0, nrow(object$validation.data)) ,
              n.covariates = length(object$parameters$covariates),

              n.cov.chosen = sum(p.lin,p.smooth,t.te),

              list.factors = object$offset.factors,

              mstop = object$gamboost.mstop,
              list.baselearners = t.boost.names,

              list.effect.size = names(object$gamback.cv[[3]]$model)[
                -match(c(object$parameters$response, "(weights)"), names(object$gamback.cv[[3]]$model) )],

              list.backward = names(object$gamback.backward[[3]]$model)[
                -match(c(object$parameters$response, "(weights)"),
                       names(object$gamback.backward[[3]]$model) )],

              list.aggregation = list.aggr,

              list.gam.final = t.final


            )

            class(l.s) <- "summary.path.geoGAM"
            return( l.s )

          } )
}

print.summary.path.geoGAM <- function(x, ... ){

  object <- x
  cat("\nGeoGAM for:", object$response )
  cat("\nFamily:", object$family)
  cat("\n---\nObservations for fitting:", object$n.obs )
  cat("\nObservations used for model validation only:", object$n.obs.val )

  cat("\n---\nPotential covariates:", object$n.covariates)
  cat("\nFinal covariates:", object$n.cov.chosen)

  cat("\n\n\n***Model selection***\n\n")


  if( "" %in% object$list.factors ){
    cat("No offset used for gradient boosting.")

  } else if( length(object$list.factors) < 1 ){
    cat("No offset chosen for gradient boosting.")

  } else {
    cat("Factor(s) chosen as offset: \n" )
    cat( paste( object$list.factors, sep = ", " ) )

  }

  cat("\n---\nCovariates chosen by gradient boosting after", object$mstop,"iterations: \n")
  cat( paste( object$list.baselearners, collapse = "\n" ) )

  cat("\n---\nCovariates chosen after cross validation of effect size:\n")
  cat( paste( object$list.effect.size, collapse = "\n" ) )


  cat("\n---\nCovariates chosen after backward selection:\n")
  cat( paste( object$list.backward, collapse = "\n" ) )

  if( length(object$list.aggregation) > 0){
  cat("\n---\nAggregated factors:")
  for( ll in 1:length( object$list.aggregation )){
    cat("\n* For factor '", names(object$list.aggregation[ll]), "' ",
        t.n <- sum(grepl("--",  object$list.aggregation[[ll]]) ), " level(s) aggregated")
    if( t.n > 0) cat(": \n ", grep("--",  object$list.aggregation[[ll]], value = TRUE)  )
  }}

  cat("\n---\nCovariates chosen for final model:\n")
  cat( paste( object$list.gam.final, collapse = "\n" ) )
  cat("\n\n")


}


# print method for summary
print.summary.geoGAM <- function(x, ... ){

  print( x$summary.gam )

  cat("---\nCross validation of final geoadditive model:\n")
  print(  round( x$summary.validation$cv, 5) )

  if( !is.null(x$summary.validation$test) ){
    cat("\n\nExternal validation of final geoadditive model:\n")
    print( round(  x$summary.validation$validation, 5 ) )
  }
}


# Plot method
plot.geoGAM <- function(x, ..., what = c("final", "path")){

  what = match.arg(what)

  switch( what,
          final = {
            # NextMethod("plot")
            plot(x$gam.final)
          },
          path = {
            # NextMethod("plot")
            plot(x$gamboost.cv)
          }
  )

}

