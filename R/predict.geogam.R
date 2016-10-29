

## Create batch predictions from a GAM fit

### 2. Create Predictions with GAM models ------


# Function to predict and to save GIS files
predict.geoGAM <- function( object,
                            newdata,
                            type = c("response", "link", "probs", "class"),
                            back.transform = c("none", "log", "sqrt"),
                            threshold = 0.5,
                            se.fit = F, ...) {



  type <- match.arg(type)
  back.transform <- match.arg(back.transform)

  if( object$parameters$family[[2]] == "Gaussian" & type %in% c("probs", "class") ){
    warning( cat( "Predictions of type = '",
                  type,"' not possible for continuous response. \nDoing predictions for type = 'response' ..", sep = "") ) }


  # Transformed calibration data
  daten <-  object$gam.final$model[, -match("(weights)", names(object$gam.final$model)), drop = F ]

  l.fact <- gsub("ag$", "", names(daten[,-1,drop = F])[ unlist( lapply( daten[,-1, drop = F], is.factor) )  ] )
  l.no.fact <- names(daten[,-1,drop = F])[ !unlist( lapply( daten[,-1, drop = F], is.factor) )  ]

  newdata <- newdata[ , c(object$parameters$coords, l.fact, l.no.fact) , drop = F]

  ## Transform the data like the original data set
  # Scale numerical data to [0;1], then center (- mean), but not the response
  # Use means / max-min as in calibration data set!
  l.sub.fact <- l.sub.fact.1 <- names(newdata) %in% l.fact
  names(l.sub.fact) <- names(newdata)
  l.sub.fact[ names(l.sub.fact) %in% c( object$parameters$response, object$parameters$coords, "(weights)") ] <- TRUE

  fun.center.v <- function(col, vali, kali, daten ){

    x.c <- ( vali[, col] - min(kali[, col]) ) /
      ( max(kali[, col]) - min(kali[, col]) )

    # Mean from [0,1]-column from calibration-data set
    x.kali.m <- ( kali[, col] - min(kali[, col]) ) /
      ( max(kali[, col]) - min(kali[, col]))

    # center with this mean
    x.sc <- scale(x.c, center = mean( x.kali.m ), scale = F)[, 1]
    return(x.sc)
  }

  if( length(l.sub.fact) > 0){
    newdata[, !l.sub.fact]  <- data.frame(
      lapply( names(newdata[, !l.sub.fact, drop = F]), FUN = fun.center.v,
              daten = daten, vali = newdata, kali = object$data) )
  }

  # Center coordinates with mean from calibration data
  if( !is.null(object$parameters$coords)){
    coord.means <-  colMeans( object$data[ , object$parameters$coords])
    newdata[, object$parameters$coords] <- data.frame(
      scale( newdata[, object$parameters$coords], center = coord.means, scale = F))
  }

  # Add intercept (and an alibi SE)
  newdata$int <- newdata$se <- rep(1, nrow(newdata) )

  # Add renamed factors for main effect
  if( length(l.sub.fact.1) > 0){
    newdata[ , paste( names( l.sub.fact.1 )[ l.sub.fact.1 ], "ag", sep = "") ] <-
      newdata[ , paste( names( l.sub.fact.1 )[ l.sub.fact.1 ] ) ]
  }

  ##  Aggregate factors accordingly
  l.factors <- names( daten[,-1, drop = F] )[ unlist( lapply( daten[, -1, drop = F], is.factor ) ) ]

  for( fact in l.factors ){

    if( fact %in% names(newdata) ){
      ii.new.fact <- which( names(newdata) == fact )
      fact.nag <- fact
    } else {
      fact.nag <- gsub("ag$", "", fact)
      ii.new.fact <- which( names(newdata) == fact.nag )
      # Add, if renamed factor
      newdata[, fact ] <- newdata[ , ii.new.fact ]
    }

    n.lev.cal <- length( lev.cal <- levels( daten[, fact] ) )
    n.lev.val <- length( lev.val <- levels( newdata[, fact.nag ]) )

    if( sum( !lev.val %in% lev.cal ) > 1){
      # Aggregate each level stepwise
      for( ii in 1:n.lev.val ){
        levels(newdata[, fact ])[ match( c( strsplit( lev.cal[ii], "--" )[[1]] ),
                                         levels( newdata[, fact ]) )]  <- lev.cal[ii]
      }
    }
  }


  se.fit <- ifelse( back.transform != "none", T, se.fit)
  type.pr <- ifelse( type %in% c("probs", "class"), "response", type)

  # normal untranformed Gaussian prediction
  pred <- predict( object$gam.final, newdata = newdata, type = type.pr, se.fit = se.fit)


  # Back transform LOG
  # approximation of unbiased backtransform using Var(XBeta) from GAM model
  # (true var(XBeta) probably larger, because model selection is ignored)
  if( back.transform == "log" ){

    pred <- exp(pred$fit - pred$se.fit*0.5)


    # Backtransformation of SQRT
    # same approximation with Var(XBeta) from GAM object
  } else if( back.transform == "sqrt" ) {

    pred <- pred$fit^2 - pred$se.fit


    # for ordinal / binomial regression
    # create "category numbers" for integer raster
  } else if( object$parameters$family[[2]] == "Binomial" & type == "class" ) {


    # optimal threshold for binarisation
    # see line 220 in 4_geoam_presentations_diverse_zhg.R
    pred <- as.integer( ifelse( pred > threshold, 1, 0) )


  } else if( object$parameters$family[[2]] == "PropOdds" & type == "class" ){

    pred <- as.integer( apply( pred,1, function(x){ min( which( cumsum(x) >= 0.5 ) ) } ) )

    if( is.factor( object$data[, object$parameters$response ] ) ){

      pred <- factor( pred, levels = sort(unique(pred)), labels = levels( object$data[, object$parameters$response ]))
    }

  }

  return(pred)
}

