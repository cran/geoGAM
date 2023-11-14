# Automatize model selectin for geoAM with boosting as presteps


### OVERVIEW ## #
# (behind the arrow = output of this step)

# 1. Search for relevant Factors -> List of rel. factors
# 2. Create a linear model from these -> fitted values
# 3. Compute Boosting with pred.values (2) as offset -> fitted values
# 4. Compute Boosting with first order interactions -> list of relevant baselearners (how?)
# 5. Fit full GAM, do backward selection by cross-validation ->
#GAM with non relevant baselearners excluded
# 6. Aggregate categories of factors (code in testing.R) -> gam.object
# 7. Compute prediction intervals with resampling original data and go
#    trough steps 1 to 6 again.

# require(grpreg)



## 1. Transform original data, center and scale -----------

f.transform.data <- function( resp, data, val.data, coords ){

  data.original <- data

  # Scale numerical data to [0;1], then center (- mean), but not the response
  l.sub.num <- unlist( lapply(data, is.factor) )
  l.sub.num[ names(l.sub.num) %in% c( resp, coords) ] <- TRUE

  # data.frame( no.transform = l.sub.num)

  fun.center <- function(x) { x.c <- (x - min(x) ) / ( max(x) - min(x))
  x.sc <- scale(x.c, center = TRUE, scale = FALSE)[, 1]
  return(x.sc)
  }
  data[, !l.sub.num ] <- lapply( data[, !l.sub.num, drop = FALSE], FUN = fun.center)

  # Center coordinates
  if( !is.null(coords)){
    data[, coords] <- data.frame(
      scale(data[, coords], center = TRUE, scale = FALSE))
  }

  # Add a separate intercept for boosting
  data$int <- rep(1, dim(data)[1] )


  ## Transform the data like the original data set
  # Scale numerical data to [0;1], then center (- mean), but not the response
  # Use means / max-min as in calibration data set!

  fun.center.v <- function(col, vali = val.data, kali = data.original ){

    x.c <- ( vali[, col] - min(kali[, col]) ) /
      ( max(kali[, col]) - min(kali[, col]) )

    # Mean from [0,1]-column from calibration-data set
    x.kali.m <- ( kali[, col] - min(kali[, col]) ) /
      ( max(kali[, col]) - min(kali[, col]))

    # center with this mean
    x.sc <- scale(x.c, center = mean( x.kali.m ), scale = FALSE)[, 1]
    return(x.sc)
  }

  val.data[, !l.sub.num]  <- data.frame(
    lapply( names(val.data[, !l.sub.num, drop = FALSE]), FUN = fun.center.v) )

  # Center coordinates with mean from calibration data
  if( !is.null(coords)){
    coord.means <-  colMeans( data.original[ , coords])
    val.data[, coords] <- data.frame(
      scale( val.data[, coords], center = coord.means, scale = FALSE))
  }

  # Add intercept
  val.data$int <- rep(1, nrow(val.data) )


  return( list( data, val.data ))
}


### 2. Search for relevant Factors -> List of rel. factors --------
#

# choose offset by group lasso
f.choose.offset <- function(resp, dat, t.sets, w.weights, v.family, verbose ){

  l.fact <- names(dat)[ unlist( lapply( dat, is.factor) ) == TRUE ]
  l.fact <- grep( paste0("^", resp, "$"), l.fact, invert = TRUE, value = TRUE)

  if( verbose > 1){ cat("\nFull list of factors: \n"); print(l.fact) }

  t.ids <- c()
  for(ii in 1:dim(t.sets)[1]){
    t.id <- which(t.sets[ii, 1:ncol(t.sets)] == 0)
    t.ids <- c(t.ids, t.id )
  }

  if( v.family[[2]] != "PropOdds" ){

    # desig matrix without intercept
    XX <- model.matrix( ~., dat[, l.fact, FALSE])[,-1]
    yy <- dat[, resp]

    # response needs to be integer with 0/1, use first level as = 0
    if( v.family[[2]] ==  "Binomial" ){ yy <- ifelse(yy == levels(yy)[1], 0, 1) }

    g.groups <- unlist( sapply(1:length(l.fact), function(n){
      rep( n, nlevels(dat[, l.fact[n]]) -1 )
    }) )

    g.cvfit <- cv.grpreg(XX, yy, group = g.groups, penalty = "grLasso",
                         family = paste(v.family[[3]][1]), cv.ind = t.ids)

    # choose lambda with min-CV + 1 SE, as in Glmnet optimum option
    l.se <- g.cvfit$cvse[ g.cvfit$min ] + g.cvfit$cve[ g.cvfit$min ]

    # in case SE > delta CV error, choose no offset
    if( l.se > g.cvfit$cve[1] ){
      l.off.fact <- f.fact.off <- ""
      p.offset.lm <- rep( 0, nrow(dat) )

    } else {

      idx.se <- min( which( g.cvfit$cve < l.se ) ) - 1

      # extract group numbers and list of chosen factors
      pp <- predict(g.cvfit, X = XX, which = idx.se, type = "groups",
                    lambda =  g.cvfit$lambda[ idx.se ])
      l.off.fact <- l.fact[as.numeric(pp)]
      f.fact.off <- paste(resp, " ~ ",paste(l.off.fact,collapse="+") )

      # create predictions for offset
      p.offset.lm <- predict(g.cvfit, X = XX, which = idx.se, type = "link",
                             lambda =  g.cvfit$lambda[ idx.se ])

      # use 1/2 * XBeta for Binomial (see Manual of mboost, page Family)
      if( v.family[[2]] == "Binomial" ){ p.offset.lm <- p.offset.lm*0.5}

    }

  } else {
    # for proportional odds model no Lasso available, use simple stepwise with BIC
    # for speed reasons and because stepwise CV of LM is much better, propably.

    f.po <- as.formula( paste( resp, "~",  paste( l.fact, collapse =  " + ")) )
    fit <- MASS::polr( f.po, data = dat, method = "logistic", model = TRUE)
    s.po <- step(fit, trace = 0, direction = "both", k = log(nrow(dat)))

    l.off.fact <- attr( s.po$terms, "term.labels")

    f.fact.off <- as.formula( paste( resp, "~",  paste( l.off.fact, collapse =  " + ")) )

    # compute glmboost, coefficents differ, but there is no other way to get
    # correct link XBeta for offset.
    p.offset.lm <- fitted( glmboost(f.fact.off, data = dat, family = PropOdds(),
                                    control = boost_control(mstop = 1000)) )
    g.cvfit <- c()
  }

  return(  list( l.off.fact, f.fact.off, p.offset.lm, g.cvfit )  )
}



# compute linear model, depending on response / family
f.lm <- function(formula, dat, v.family ){
  if(v.family[[2]] == "Binomial"){
    glm( as.formula(formula), data = dat, family = binomial(link= "logit"))

  } else if(v.family[[2]] == "PropOdds"){
    polr(as.formula(formula), data = dat, method = "logistic", model = TRUE)

  } else if(v.family[[2]] == "Multinomial"){

  } else { # Gaussian
    lm( formula, data = dat)
  }
}


# predict linear model according to response type
f.predict <- function(model, newdat, v.family){
  if(v.family[[2]] == "Binomial"){
    predict(model, newdat, type = "response")

  } else if(v.family[[2]] == "Multinomial" | v.family[[2]] == "PropOdds"){
    predict(model, newdat, type = "probs")

  } else { # Gaussian
    predict(model, newdat)
  }
}


f.lm.cv <- function( lm.model, t.sets, resp, v.family, n.cores){
  d.dat <- lm.model$model
  form.lm <- formula(lm.model)

  t.ids <- c()
  for(ii in 1:dim(t.sets)[1]){
    t.id <- which(t.sets[ii, 1:ncol(t.sets)] == 0)
    t.ids <- c(t.ids, t.id )
  }

  my.cv.lm.pred <- function(ii, t.t = t.ids, v.family) {
    d.dat <- lm.model$model
    # Calc model on training set
    m.lm.cv <- f.lm( form.lm, dat = d.dat[ which(t.t != ii),  ] , v.family = v.family)

    # Predict on test set
    return( f.predict( m.lm.cv, d.dat[ which(t.t == ii),  ], v.family = v.family ) )
  }

  l.pred <- mclapply(1:ncol(t.sets), my.cv.lm.pred, v.family = v.family,
                     mc.preschedule=FALSE, mc.set.seed=FALSE, mc.cores = n.cores)
  # l.pred <- lapply(1:10, my.cv.lm.pred)

  t.val <- f.prepare.validation.cv(dat = d.dat, pred = l.pred, resp = resp,
                                   v.family = v.family, t.sets = t.sets)

  stats <- round( f.validate.predictions( t.val[[1]], t.val[[2]], t.val[[3]]) )

  return( stats  )
}


# 2 stupid functions -> TODO replace by something more efficient!
f.prepare.validation.cv <- function( dat, pred, resp, t.sets, gam.fit = 0, v.family, se.fit = FALSE){
  # just take as t.fitted the one that is original response,
  # like that RMSE = BS

  t.ids <- c()
  for(ii in 1:dim(t.sets)[1]){
    t.id <- which(t.sets[ii, 1:ncol(t.sets)] == 0)
    t.ids <- c(t.ids, t.id )
  }

  t.fitted <- t.se <- rep(NA, nrow(dat))
  if( (v.family[[2]] == "Multinomial" & gam.fit == 0) | v.family[[2]] == "PropOdds"){
    t.fitted <- matrix( rep(t.fitted, max(as.numeric(dat[, resp])) ), ncol = max(as.numeric(dat[, resp]))  )
    for(ii in 1:ncol(t.sets)){

      # for Ignorace score!
      # match the prediction probability of the actual observed class
      #       l.pred.pr <- l.pred[[ii]]
      #       l.meas.cl <- as.numeric( d.dat[ which(t.ids == ii), resp ] )
      #       t.fitted[which(t.ids == ii)] <- sapply(seq(1,nrow(l.pred.pr)),
      #                                              function(x){ l.pred.pr[x, l.meas.cl[x]] })

      t.fitted[ which(t.ids == ii), ] <- pred[[ii]]
    }
    t.observed <- dat[, resp ]
    t.measure <- ifelse(v.family[[2]] == "Multinomial", "bss", "rps")
  } else {

    for(ii in 1:ncol(t.sets)){
      if(se.fit){
        t.fitted[which(t.ids == ii)] <- pred[[ii]]$pred
        t.se[ which(t.ids == ii)] <- pred[[ii]]$pred.se
      } else {
        t.fitted[which(t.ids == ii)] <- pred[[ii]]
      }
    }
    t.measure <- "st"

    if( v.family[[2]] == "Binomial" | gam.fit == 1 ){
      t.observed <- as.numeric(dat[, resp ])-1 # reference category = 2. level, change to 1
      t.measure <- "bs"
    } else {
      t.observed <- dat[, resp ]
    }
  }
  return( if(se.fit){ list(t.observed,t.fitted,t.measure, t.se) } else{ list(t.observed,t.fitted,t.measure) } )
}

f.prepare.validation.ext <- function( dat, pred, gam.fit = 0, v.family ){

  if( (v.family[[2]] == "Multinomial" & gam.fit == 0) | v.family[[2]] == "PropOdds"){
    t.observed <- dat
    t.measure <- ifelse(v.family[[2]] == "Multinomial", "bss", "rps")

  } else{
    t.measure <- "st"
    if( v.family[[2]] == "Binomial" | gam.fit == 1 ){
      t.observed <- as.numeric(dat)-1 # reference category = 2. level, change to 1
      t.measure <- "bs"
    } else {
      t.observed <- dat
    }
  }
  return( list(t.observed,pred,t.measure) )
}

## 3. Compute Boosting  --------------------------------------------------------
# with pred.values (2) as offset -> fitted values

f.boost <- function( resp, dat, p.offset.lm, val.data, f.fact.off, v.family,
                     gam.selection = 1, t.sets, non.stat, max.stop, coords){

  # Cutoff for baselearner selection for GAM
  sel.haupt <- 1
  sel.inter <- 1.2 # choose a bit less interactions

  l.fact <- names(dat)[ unlist( lapply( dat, is.factor) ) == TRUE ]
  l.fact <- grep( paste0("^", resp, "$"), l.fact, invert = TRUE, value = TRUE)

  l.num <- names(dat)[ unlist( lapply( dat, is.factor) ) == FALSE ]

  # Create Boosting-Baselearner-Formula

  # for Multnom expand baselearners (p. 33 of mboost manual)
  if( v.family[[2]] == "Multinomial"){
    X0 <- K0 <- diag(nlevels(dat[, resp]) - 1)
    colnames(X0) <- levels(dat[, resp])[-nlevels(dat[, resp])]
    t.exp <- "%O% buser(X0, K0, df = 2)"
  } else {
    t.exp <- ""
  }

  # Verschieden Formel-Bestandteile erstellen
  # Response
  f.resp <- paste(resp, " ~ ")

  # Create Intercept separate
  f.int <- paste( "bols(int, intercept = FALSE, df = 1)", t.exp)

  # Spatial Term
  if( !is.null(coords)) {
    f.spat <- paste( "bspatial(", paste(coords, collapse = ","), ", df = 5, knots = 12)", t.exp)
    } else{ f.spat <- c() }


  ## Factors, group if less than 6 levels (dummy coding for df=5)
  f.fact <- f.short <- c()
  for(fact in l.fact){
    if( length( levels(dat[, fact]) ) > 5) {
      if(length( levels(dat[, fact]) ) > 5 & v.family[[2]] == "Multinomial"){
        # factors cannot be fitted without intercept if nlevels > df*2 (?)
        f.fact <- c(f.fact, paste(" bols(", fact, ", df = 5, intercept = TRUE)", t.exp ))
      } else {
        f.fact <- c(f.fact, paste(" bols(", fact, ", df = 5, intercept = FALSE)", t.exp ))
      }
    } else {
      f.short <- c(f.short, fact)
    }
  }
  f.fact <- paste( f.fact, collapse = "+")

  # Assign SHORT level factors to groups
  # assume neighborhood in data frame, means similar topic to form a group
  f.sh <- c()
  if( length(f.short) >= 2 ){
    while( length(f.short) > 0 ){
      n.levels <- unlist( lapply( dat[, f.short], nlevels) )
      ii <- 1
      t.df <- 0
      # add up degrees of freedom
      while( t.df < 6 ){
        t.df <- sum( t.df, n.levels[ii]-1)
        ii <- ii+1
      }
      # in case factors with not enough levels left over (sum<6), add to last baselearner
      if( ii-1 < length(n.levels) ){
        t.leftover <- sum(n.levels[ii:length(n.levels)]) - length( ii:length(n.levels) )
        if( t.leftover < 6 ){ii <- length(n.levels)+1}
      }

      f.sh <- c( f.sh,
                 paste( "bols(", paste(f.short[1:(ii-1)], collapse = ", "),
                        ", df = 5, intercept = FALSE)", t.exp) )
      # remove used factors
      f.short <- f.short[ -c(1:(ii-1)) ]
    }
    f.factor <- paste(  paste( f.sh, collapse = "+") ) # short level factors

  } else if( length(f.short) == 1){
    f.factor <- paste(" bols(", f.short, ", df = 5, intercept = FALSE)", t.exp )
  } else {
    f.factor <- c()
  }

  if( nchar(f.fact) > 0 ){ f.factor <- paste( f.factor, f.fact, sep = "+")}  # long level factors
  f.factor <- sub( "^\\+|^ \\+", "", f.factor)


  ## Splines
  # Numerische bbs()
  f.num.bbs <- paste( "bbs(",
                      paste(
                        l.num[ -which( l.num %in% c(resp, coords, "int"))],
                        collapse = paste(", center = FALSE, df = 5)", t.exp, " + bbs("
                        )),
                      ", center = FALSE, df = 5)", t.exp,
                      sep = ""
  )

  # create Apex interaction for images 2014 and 2013

  # Numerische by Indikator bbs(); # bbs(x, by = z, df = 5)
  l.apx <- grep( "apxb", l.num, value = TRUE)
  int.apx <- grep( "apxi", l.fact, value = TRUE)
  if( length( l.apx) > 0 ){
    f.num.bbs.apx <- paste( "bbs(",
                            paste(
                              l.apx,
                              collapse = paste(", by=", int.apx, ", center = FALSE, df = 5)", t.exp, " + bbs("
                              )),
                            ", by=", int.apx, ", center = FALSE, df = 5)", t.exp,
                            sep = ""
    )
    f.num.bbs <- paste(f.num.bbs, "+", f.num.bbs.apx)
  }


  # Numerische vs. bspatial
  if( !is.null(coords) ){
  l.sub.num <- l.num[ -which( l.num %in% c(resp, coords, "int") )]
  f.spat.int <- paste( "bspatial(", paste(coords, collapse = ","), ", by=",
                       paste(
                         l.sub.num,
                         collapse = paste(
                           ", df=5, knots = 12) ", t.exp , " + bspatial(", paste(coords, collapse = ","), ", by=" )
                       ),
                       ", df=5, knots = 12)", t.exp,
                       sep = ""
  ) } else{ f.spat.int <- c() }

  # Create formula
  if( non.stat == 0 ){
    a.form <- as.formula( paste( f.resp, paste( c(f.int, f.num.bbs, f.spat,
                                                  f.factor),
                                                collapse = "+")) )
  } else {
    a.form <- as.formula( paste( f.resp, paste( c(f.int, f.num.bbs, f.spat,
                                                  f.factor , f.spat.int), collapse = "+" ) ) )
  }


  m.boost.1 <- gamboost( a.form, data = dat,
                         control = boost_control(mstop = max.stop),
                         offset = p.offset.lm,
                         weights = dat$w.weights,
                         family = v.family[[1]])

  # Cross-Validation to find mstop
  t.sets.m <- t.sets
  # cross validation sets need to be expanded to the Multinomial model frame
  if( v.family[[2]] == "Multinomial"){
    t.sets.m <- t.sets[ rep(1:nrow(t.sets), nlevels(dat[, resp])-1 ),  ]
  }
  m.boost.1.cv <- cvrisk(m.boost.1, folds = t.sets.m, papply = lapply)

  # # Use Mstop at optimum, or if small at minimum

  if( mstop( m.boost.1.cv) > 20 ){
    # Use -% of delta CV error as mstop
    redu.mstop <- ifelse( sum(p.offset.lm) == 0, 0.005, 0.01 )

    dlta10 <- mean( m.boost.1.cv[, mstop( m.boost.1.cv  ) ] ) +
      ( mean( m.boost.1.cv[, 1] ) - mean( m.boost.1.cv[, mstop( m.boost.1.cv  ) ] ) ) * redu.mstop

    l.se.stop <- max( which( (colMeans(m.boost.1.cv ) > dlta10 )[1:mstop( m.boost.1.cv) ] ) )
  } else{
    l.se.stop <- mstop( m.boost.1.cv)
  }

  bla <- capture.output(
    l.baselearn.haupt <- names(coef(m.boost.1[l.se.stop] )), type = "message" )



  l.baselearn.interact <- m.boost.2 <- m.boost.2.cv <- l.se.stop.2 <- c()


  l.valid <- list(
    RMSE_Boosting_cross_valid_step1 = round( c( sqrt( min( colMeans( m.boost.1.cv))) ), 5)#,
    # MSE_Boosting_cross_valid_step2 = round( c( min( colMeans( m.boost.2.cv)) ), 5)
  )


  # create emtpy prediction vector for no validation set available
  t.val <- list(); t.val[[2]] <- 0 #}

  return( list( l.baselearn.haupt, l.baselearn.interact, c(), l.valid,
                m.boost.1, m.boost.1.cv, m.boost.2, m.boost.2.cv, val.data,
                mstop_1_2 = c( l.se.stop, l.se.stop.2 ), pred = t.val[[2]]))

}



## ---------- Function to build tuning parameter matrix ---------------
#  Step from Boosting to GAM: for cross-validation by GAM

# !! ohne interaktionen, so far


# a) compute magnitude of residual plot
f.comp.magn <- function( w, boost.obj ){ # w = baselearner-string

  data <- model.frame(boost.obj, which = w)[[1]]
  pr <- predict(boost.obj, newdata = data, which = w)

  # Create residual plot data points and spline fun
  aa <- sort(data[[1]])
  bb <- pr[order(data[[1]], na.last = NA)]

  # remove outliers to compute magnitude of residual plot
  # use definition of boxplot
  outl <- boxplot(aa, plot = FALSE)$out
  bb.nout <- bb[ !aa %in% outl ]

  magn <- max(bb.nout)+abs(min(bb.nout) )
  # use selection probablity as weights on magnitude
  selpr <- as.numeric( summary(boost.obj)$selprob[
    match( w, attr( summary(boost.obj)$selprob, "names") ) ] )

  fx.spline <- splinefun(aa, bb)

  l.2 <- splinefun( aa, fx.spline(aa,deriv=2))

  l.intsum <- 0

  return( c( l.intsum, magn) )
}


#
# ll <- coef(m.boost.1 )[[ii]]
# df$diff[ii] <- max(ll) + abs( min( ll ))


# b) add matrix together
f.tuning.param <- function( boost.object, dat, resp, non.stat, t.sets){

  # List of splines baselearner
  bla <- capture.output(
    l.bbs <- names(coef(boost.object))[grep("bbs", names(coef(boost.object)))], type = "message" )

  if( length(l.bbs) > 0 ){
    # degree of linearity and magnitude for each
    l.lin <- matrix(
      unlist(
        lapply(l.bbs, f.comp.magn, boost.obj = boost.object) ),
      ncol = 2, byrow = TRUE )

    # selection probabilty for each
    l.selpr <- as.numeric( summary(boost.object)$selprob[
      match( l.bbs, attr( summary(boost.object)$selprob, "names") ) ] )

    df <- data.frame( name = l.bbs, selprob = l.selpr,
                      linear = l.lin[,1], magnitude = l.lin[,2], stringsAsFactors = FALSE)
  } else {

    df <- data.frame( name = character(), selprob = numeric(), linear = numeric(),
                      magnitude = numeric(), stringsAsFactors = FALSE)
  }

  # bols effektstaerke
  bla <- capture.output(
    l.bols <- names(coef(boost.object))[grep("bols", names(coef(boost.object)))], type = "message" )
  if( length(l.bols) == 0 ) {
    bols.magn <- c()
  } else {

    bols.magn <- unlist( lapply( 1:length(l.bols), function(ii){

      # split in single factors (that were grouped by bols to reach 5 df)
      l.f <- strsplit( gsub("bols\\(|intercept = FALSE, df = 5\\)", "", l.bols[ii]), ", " )[[1]]

      # get coefficient values
      bla <- capture.output(
        ll <- coef(boost.object)[
          names(coef(boost.object) ) == l.bols[ii] ][[1]], type = "message" )

      l.magn <- unlist( lapply( 1:length(l.f), FUN = function(ff){
        l.c <- ll[ grep( l.f[ff], names(ll)) ]

        # add reference level for first factor
        if( length(l.c) == 1){ l.c <- c(0, l.c)}

        magn <- max(l.c) + abs( min( l.c))
        names(magn) <- l.f[ff]

        return(magn)
      } ) )

      # apply selprop weight

      selpr <- as.numeric( summary(boost.object)$selprob[
        match( l.bols[ii], attr( summary(boost.object)$selprob, "names") ) ] )

      #l.magn <- l.magn*selpr

      return( l.magn )
    }) )

    df.bols <- data.frame( name = names(bols.magn), selprob = rep(1, length(bols.magn)),
                           linear = rep(1, length(bols.magn)), magnitude = bols.magn,
                           stringsAsFactors = FALSE )
    rownames(df.bols) <- NULL
    df <- rbind( df, df.bols)
  }


  # bspatial effektstaerke
  l.bspat <- names(coef(boost.object))[grep("bspat", names(coef(boost.object)))]
  if( non.stat == 1 | length(l.bspat) > 0 ){

    bspat.magn <- unlist( lapply( 1:length(l.bspat), function(ii){

      bla <- capture.output(
        ll <- unlist( coef(boost.object)[
          names(coef(boost.object) ) == l.bspat[ii] ] ), type = "message")

      magn <- max(ll) + abs( min( ll ))

      selpr <- as.numeric( summary(boost.object)$selprob[
        match( l.bspat[ii], attr( summary(boost.object)$selprob, "names") ) ] )

      #magn <- magn*selpr
      return( magn )
    }) )

    df.spat <- data.frame( name = l.bspat, magn.spat = bspat.magn, stringsAsFactors = FALSE )
    dim.1 <- unique( c( 0, as.numeric( round(quantile(df.spat$magn.spat, probs = seq(0,1,0.1)), 4) ) ) )
    dim.1[length(dim.1)] <- dim.1[length(dim.1)]*2

    # limit the number of magnitudes to cross validate
    if( length(dim.1) > 10 ){
      dim.1 <- dim.1[ c(TRUE,TRUE,TRUE,1:c(length(dim.1)-4) %% 2 == 0, TRUE) ]
    }

  } else{
    dim.1 <- 0
  }

  # set up tuning parameter matrix in three dimensions
  #   dim.1 <- c( 0, as.numeric( round(quantile(df.spat$magn.spat, probs = seq(0,1,0.1)), 4) ) )
  #   dim.1 <- c( 0, sort( unique( round( df$selprob, 3  ) ) ) )
  # dim.2 <- c( 0, as.numeric( round( quantile(df$linear, probs = seq(0,1,0.25)), 4 )  ) )
  # dim.3 <- c( 0, as.numeric( round(quantile(df$magnitude, probs = seq(0,1,0.1)), 4) ) )
  # take all values as thresholds
  dim.3 <-  c( 0, unique(sort(as.numeric( round(df$magnitude+0.00001, 6) ) )))
  dim.3[length(dim.3)] <- dim.3[length(dim.3)]*2

  if( length(dim.3) > 21 ){
    # do only limited number of cross validation models
    # first 15 steps of magnitude
    dim.red <- sort(dim.3, decreasing = TRUE)[c(1:10, 12, 14, 16, 18, 20)]

    # then only another 20 models.
    d <- ceiling( (length(dim.3) - 15) / 20 )
    idx <- cumsum(rep(d, 20))[ cumsum(rep(d, 10)) < (length(dim.3) - 15 )]

    dim.red <- c( dim.red,
                  sort(dim.3, decreasing = TRUE)[11:length(dim.3)][idx]
    )
    dim.3 <- c(0, sort( dim.red ))
  }

  # do not CV linearity threshold
  # dim.2 <- max(df$linear) + 1
  #   dim.1 <- 0
  if( non.stat == 1 | length(l.bspat) > 0 ){
    l.names.coef <- c( df$name, df.spat$name)
  } else {
    l.names.coef <- df$name
  }
  m.res <- matrix(0, nrow = length(dim.1)*length(dim.3), ncol = length(l.names.coef))

  #l.selpr.bspat <- summary(boost.object)$selprob[ grep("bspat|bols",
  #                                                     names(summary(boost.object)$selprob))]

  l.modelle <- c()
  l.modnr <- 0
  for( sp.mag in dim.1 ){
    for( mag in dim.3 ){

      l.modnr <- l.modnr + 1

      # collect modell parameter in a vector
      l.modelle <- c(l.modelle, paste( "spat", sp.mag, "-bbs", mag, sep = "--"))

      # choose relevant bbs and bols-baselearner
      #l.bbs.rel <- df$name[ df$selprob >= sprob & df$magnitude >= mag]
      l.bbs.rel <- df$name[ df$magnitude >= mag]

      # choose relvant bols/bspatial baselearner
      #l.bsp.rel <- names( l.selpr.bspat )[ l.selpr.bspat >= sprob]
      if( non.stat == 1 | length(l.bspat) > 0 ){
        l.bsp.rel <- df.spat$name[ df.spat$magn.spat >= sp.mag ]
      } else {
        l.bsp.rel <- c()
      }

      # set relevant ones to 1
      m.res[ l.modnr, match( c(l.bbs.rel, l.bsp.rel), l.names.coef)] <- 1
    }
  }

  # go through linearity thresholds
  #   l.modelle.lin <- c()
  #   m.res.fin <- matrix(, nrow = 0, ncol = ncol(m.res))
  #   for( lin in dim.2 ){
  #     m.res.lin <- m.res
  #     l.bbs.linear <- df$name[ df$linear >= lin ]
  #
  #     for( ll in l.bbs.linear){
  #       n.col <- match( ll, l.names.coef)
  #       m.res.lin[ m.res.lin[, n.col] == 1, n.col ] <- 2
  #     }
  #
  #     l.modelle.lin <- c( l.modelle.lin,
  #                         paste( l.modelle , "-lin.", round(lin,3), sep = ""))
  #
  #     m.res.fin <- rbind( m.res.fin, m.res.lin)
  #   }
  #
  #   d.res.fin <- data.frame( modell = l.modelle.lin, m.res.fin, stringsAsFactors =FALSE )
  d.res.fin <- data.frame( modell = l.modelle, m.res, stringsAsFactors = FALSE)
  names(d.res.fin)[-1] <- l.names.coef

  d.cv.list <- d.res.fin[ !duplicated(d.res.fin[,-1]), ]


  # iv) set null, where p > n, that cannot be fitted anymore
  # for larger data sets, put limit to max 800 (otherwise model selection to slow)
  # p only 80% of nrow(data)
  n.limit <- ifelse( 800 > min(colSums(t.sets))*0.8, min(colSums(t.sets))*0.8, 800 )
  l.remove <- c()
  for( l.row in 1:nrow(d.cv.list)){

    d.mod <- d.cv.list[l.row, ]
    l.baselearn  <- names(d.mod[,-1])[ d.mod[,-1] == 1 ]
    l.baselearn.haupt <- l.baselearn[ grep( "by", l.baselearn, invert = TRUE ) ]
    l.baselearn.interact <- l.baselearn[ grep( "by", l.baselearn ) ]

    l.fact <- names(dat)[ unlist( lapply( dat, is.factor) ) == TRUE ]
    l.fact <- grep( paste0("^", resp, "$"), l.fact, invert = TRUE, value = TRUE)


    ### Check for degrees of freedom
    # GAM cannot be fitted with coefficients > data (for 10fold CV)
    # reduce baselearner list until criteria met

    # Number of splines (incl interactions)
    n.bbs.h <- sum( grepl( "bbs", c( l.baselearn.haupt ) ) )
    # Number of spatial terms
    n.spat <- sum( grepl( "bspat", c( l.baselearn.haupt, l.baselearn.interact) ) )
    # Number of all factor levels
    n.categ <- sum( unlist( lapply( dat[, l.fact, FALSE], nlevels) ) )
    # Number of factor interaction
    n.inter <- 0
    if( length(l.fact) > 0 ){
      for( ff in 1:length(l.fact) ){
        l.bbs <- sum( grepl( "bbs",
                             l.baselearn.interact[ grepl( l.fact[ff], c( l.baselearn.interact)  ) ] )
        ) * 16 * nlevels( dat[ , l.fact[ff]])
        l.bols <- sum( grepl( "bols",
                              l.baselearn.interact[ grepl( l.fact[ff], c( l.baselearn.interact)  ) ] )
        ) * nlevels( dat[ , l.fact[ff]])

        n.inter <- n.inter + l.bols + l.bbs
      }
    }
    # bbs interactions with numeric
    l.bbs.num <- 0
    for( bb in l.baselearn.interact[ grepl( "bbs", l.baselearn.interact) ] ) {
      st <- strsplit( bb, split ="\\(|\\,|\\=", perl = TRUE)[[1]][4]
      st <- gsub(" ", "", st, fixed = TRUE)
      l.bbs.num <- l.bbs.num +  ( match( st, l.fact , nomatch = 0 ) == 0 ) *16
    }

    # Approximate number of coefficients (+ Intercept etc.)
    # + penalty for interactions (150 free rows)
    n.tot <- n.spat*255 + n.bbs.h*15 + n.categ + n.inter + l.bbs.num + 20

    if( n.tot > n.limit ){ l.remove <- c( l.remove, l.row) }

  }

  if( length(l.remove) > 0){ d.cv.list <- d.cv.list[ -l.remove, ] }

  # check for NA in rows
  if( ncol(d.cv.list) > 2 ){
    if( sum( apply( d.cv.list[, -1], 1, sum) == 0 ) >=1 ){
      l.rm <- which( apply( d.cv.list[, -1], 1, sum) == 0 )
      d.cv.list <- d.cv.list[ -l.rm, ]
    }
  }
  return( d.cv.list )

}




##  ---------- Function Rebuild formula from mboost to mgcv -----------

f.transform.formula <- function(
  l.baselearn.haupt,    # character vector of beaslearners of main effects
  l.baselearn.interact, # character vector of baselearners of interactions
  l.offset,             # char vector of factors in offset
  daten,                  # transformed original data.frame
  resp,
  l.linear = c(),
  w.weights,
  coords
  #
  # returns: list of 2: formula as input for gam() from mgcv, vector with lambdas
  #
){
  n.knots <- 12
  # kk <- n.knots + 4
  #wei <- w.weights
  ### Transform main effects

  # Prepare empty lists and boosting-lambda-weights
  l.lin.haupt <- l.spli.haupt <- l.lambd.haupt <- c()

  # Add offset factors
  l.lin.haupt <- l.offset

  # Add bbs that should be linearized
  for( baselearner in l.linear ){
    l.lin.haupt <- c(l.lin.haupt,
                     strsplit( baselearner, split ="\\(|\\,", perl = TRUE)[[1]][2])
  }


  # Add main effects
  for( baselearner in l.baselearn.haupt ){

    # a) linear main effects (bols)
    if( grepl("bols", baselearner) == TRUE){
      l.lin.haupt <- c(l.lin.haupt,
                       strsplit( baselearner, split ="\\(|\\,", perl = TRUE)[[1]][2])

      # If several factors in one baselearner, take them all
      len.fact <- length( strsplit( baselearner, split ="\\(|\\,", perl = TRUE)[[1]] )
      if(  len.fact > 4 ){
        for( ii in  (len.fact-4) ){
          l.lin.haupt <- c(l.lin.haupt,
                           strsplit( baselearner,
                                     split ="\\(|\\,", perl = TRUE)[[1]][c(2+ii)])
        }
      }
    }

    # b) splines main effects (bbs)
    if( grepl("bbs", baselearner) == TRUE){
      b.name <- strsplit( baselearner, split = "\\(|\\,", perl = TRUE)[[1]][2]
      name.new.h <- paste("s(", b.name, ", bs = 'ps', k = 16, m = c(3,2) )" )
      l.spli.haupt <- c( l.spli.haupt, name.new.h )

      # Extract lambdas from boosting
      # Example: <- extract( bbs(daten$promitt, df = 5, center = FALSE)$dpp(w.weights), what = "lambda")
      lambd <- extract( eval( parse( text = paste(
        sub( b.name, paste( "daten$", b.name, ", knots = 12", sep = ""), baselearner )
        , "$dpp(w.weights)", sep = ""))), what = "lambda" )
      names(lambd) <- c( name.new.h )
      l.lambd.haupt <- c( l.lambd.haupt, lambd )
    }

  } # end loop through baselearner main effect

  # c) Add bspatial part, if chosen by boosting selection
  if( sum( ( grepl("bspat", l.baselearn.haupt) +
             grepl("by", l.baselearn.haupt) ) == 1 ) == 1 ){
    f.spat.haupt <- paste0("te(", paste(coords, collapse = ","), ", bs = 'ps', m = c(3,2), k = 16)")

    # Add lambdas for spatial term
    l.sp <- extract( bspatial(daten[, coords[1]], daten[, coords[2]],
                              df = 5, center = FALSE, knots = 12)$dpp(w.weights), what = "lambda")
    names( l.sp ) <- f.spat.haupt
    l.lambd.haupt <- c( l.lambd.haupt, rep(l.sp, 2))


  } else { f.spat.haupt <- c() }

  # Add bspatial interactions, if chosen by either 1. or 2. boosting
  l.spat.inter <- l.lambd.spat.inter <- c()
  for( baselearner in c( l.baselearn.haupt, l.baselearn.interact )){

    if( sum(grepl("bspat",baselearner)+grepl("by",baselearner)) == 2  ){

      l.split <- strsplit( baselearner, split = "\\(|\\,|\\=", perl = TRUE)[[1]]
      name.new <- paste("te(", paste(coords, collapse = ","), ", by = ", l.split[5],
                        ", bs = 'ps', k = 6, m = c(3,2))" )
      l.spat.inter <- c( l.spat.inter, name.new)

      # Extract lambdas from boosting
      # Example: extract( bbs(daten$promitt, by = dat$cp_ph,
      #                   df = 5, center = FALSE)$dpp(wei), what = "lambda")
      # workaround, because of missinterpret of "eval", save covariates in new object
      by.covar <- daten[ , gsub(" ", "", l.split[5], fixed = TRUE)  ]
      i.string <- paste0("bspatial(daten$", coords[1] ,
                         ",daten$",coords[2],",by=by.covar,df=5,center=FALSE,knots=6)$dpp(w.weights)")

      lambd <- extract( eval( parse( text = i.string ) ) , what = "lambda" )
      names(lambd) <- c( name.new )
      l.lambd.spat.inter <- c( l.lambd.spat.inter, rep(lambd, 2) )
    }
  }



  # Concatenate main
  # Remove duplicates
  l.lin.haupt <- c( l.lin.haupt, l.spat.inter)
  l.lin.haupt <- gsub( " ", "", l.lin.haupt)
  l.lin.haupt <- l.lin.haupt[ !duplicated(l.lin.haupt) ]

  f.lin.haupt <- paste( l.lin.haupt, collapse = " + ")
  f.spli.haupt <- paste( l.spli.haupt, collapse = " + ")


  ### Transform interactions

  # Prepare empty lists and boosting-lambda-weights
  l.lin.inter <- l.spli.inter <- l.lambd.inter <- c()

  for( baselearner in l.baselearn.interact ){

    # a) linear effects (bols-interactions)
    if( grepl("bols", baselearner) == TRUE){
      l.split <- strsplit( baselearner, split ="\\(|\\,|\\=", perl = TRUE)[[1]]

      # if one partner is a factor, rename factor for later separate aggregation
      if(  is.factor(daten[, gsub(" ", "", l.split[4])] ) ){
        n.original <- gsub(" ", "", l.split[4])
        l.split[4] <- paste( n.original, "_",
                             substr(l.split[2], 1, 4 ), sep = "" )
        daten[, l.split[4] ] <- daten[ , n.original]
      }

      l.lin.inter <- c(l.lin.inter, paste(l.split[2], ":", l.split[4], sep = "") )
    }

    # b) splines effects (bbs-interactions)
    if( grepl("bbs", baselearner) == TRUE){
      l.split <- strsplit( baselearner, split = "\\(|\\,|\\=", perl = TRUE)[[1]]
      name.new <- paste("s(", l.split[2], ", by = ", l.split[4],
                        ", bs = 'ps', k = 16, m = c(3,2))" )
      l.spli.inter <- c( l.spli.inter, name.new)

      # Extract lambdas from boosting
      # Example: extract( bbs(daten$promitt, by = dat$cp_ph,
      #                   df = 5, center = FALSE)$dpp(w.weights), what = "lambda")
      # workaround, because of missinterpret of "eval", save covariates in new object
      i.covar <- daten[ , l.split[2] ]
      by.covar <- daten[ , gsub(" ", "", l.split[4], fixed = TRUE)  ]
      i.string <- "bbs(i.covar,by=by.covar,df=5,center=FALSE,knots=12)$dpp(w.weights)"

      lambd <- extract( eval( parse( text = i.string ) ) , what = "lambda" )
      names(lambd) <- c( name.new )

      # if interaction is a factor, add as many duplicates of lambda as there are levels
      n.levels <- 1
      d.values <- eval( parse(text= paste( "daten$", l.split[4], sep = "") ) )
      if( is.factor(d.values) ) {
        n.levels <- length( levels( d.values )) }
      l.lambd.inter <- c( l.lambd.inter, rep(lambd, n.levels) )
    }

  } # end loop through baselearner interactions


  # Concatenate interactions
  f.lin.inter <- paste( l.lin.inter, collapse = " + ")
  f.spli.inter <- paste( l.spli.inter, collapse = " + ")


  # Put everything to one formula
  # if( v.family[[2]] == "PropOdds" ){ r.response <- paste("as.numeric(", r.response, ")", sep = "")}
  l.final <- paste( paste( resp, " ~ 1 +"), paste( f.lin.haupt,
                                                   f.lin.inter, f.spli.haupt, f.spat.haupt, f.spli.inter, sep = "+") )

  l.list <- c( l.lin.haupt, l.lin.inter, l.spli.haupt, f.spat.haupt, l.spli.inter)

  # Remove NA at end
  l.final <- gsub( "NA", "", l.final)
  # Remove "+" at the end
  while( substr(l.final, nchar(l.final), nchar(l.final)) %in% c( " ", "+")){
    l.final <- substr(l.final, 1, nchar(l.final)-1) }

  f.final <- as.formula( l.final )

  l.lambdas.final <- c(l.lambd.haupt, l.lambd.spat.inter, l.lambd.inter)


  return( list( f.final, l.lambdas.final, daten) )

}


## ------- -Function to choose term to drop, returns new formula and lambdas----

f.gam.find.drop <- function(
  m.gam.full, v.family, coords ){

  # 1) p-values of linear terms
  # choose the one to drop

  # List of Elements in Formula
  l.terms <- strsplit( paste( formula(m.gam.full)[3]), "\\+")[[1]]
  # remove spaces and newLines
  l.terms <- gsub(" ", "", l.terms, fixed = TRUE)
  l.terms <- gsub( "\n", "", l.terms, fixed = TRUE)
  # remove intercept
  l.terms <- l.terms[ !l.terms %in% c("1", "")]

  # read original data from object
  dat <- m.gam.full$model

  names.lin.coef <- names(summary( m.gam.full)$p.pv)
  names.spline.coef <- names( summary( m.gam.full)$chi.sq )

  p.sum <- data.frame( name = l.terms,
                       name.output = rep("", length(l.terms)),
                       p = rep(1, length(l.terms)),
                       stringsAsFactors = FALSE)

  # Iterate through every term
  for( term.ii in 1:length(p.sum$name) ){
    term <- p.sum$name[ term.ii ]

    # F-Test linear (p-values)
    p.terms <- anova(m.gam.full)$pTerms.pv

    if( term %in% names(p.terms) ){
      # read p value from Anova-object
      l.index <- which( names(p.terms) == term )
      p.sum$name.output[ term.ii] <- names(p.terms)[l.index]
      p.sum$p[ term.ii ] <- as.numeric( p.terms[l.index] )
    }

    # F-Test linear interactions (formula swaps terms sometimes!!)
    if( grepl(":", term) & p.sum$name.output[ p.sum$name == term ] == "" ){
      # split and swap name
      swap.name <- paste( strsplit(term, ":")[[1]][2], ":",
                          strsplit(term, ":")[[1]][1], sep = "")

      # read p value from Anova-object
      l.index <- which( names(p.terms) == swap.name )
      if( is.na(l.index) ) l.index <- which( names(p.terms) == term )

      p.sum$name.output[ term.ii ] <- swap.name
      p.sum$p[ term.ii ] <- as.numeric( p.terms[l.index] )
    }

    # Smooth terms
    if( grepl( "s\\(", term) ){

      # Smooth interaction terms
      if( grepl( "by", term) ){
        name.by <- strsplit( term, "\\,|\\=")[[1]][3]

        # Create name as on output and save in object (for lambdas extraction)
        p.sum$name.output[ term.ii ] <-
          paste( strsplit( term, "\\,")[[1]][1], "):", name.by, sep = "")


        # in case interaction is factor
        if( is.factor( dat[, name.by] ) ){

          # Compute model without group of smooth with this factor (Nullmodel)

          # Remove from formula
          term.drop <- gsub("\\\n", "", term)
          f.terms <- terms(formula(m.gam.full))
          n.drop <- which(gsub(" ", "", attr(f.terms, "term.labels"),
                               fixed = TRUE) == term.drop)
          f.form.0 <- formula( drop.terms( terms(formula(m.gam.full)),
                                           dropx = n.drop, keep.response = TRUE) )
          # Drop lambdas
          l.lambd.0 <- m.gam.full$full.sp[ grep( p.sum$name.output[ term.ii ],
                                                 names( m.gam.full$full.sp ), fixed = TRUE, invert = TRUE) ]
          # Compute Nullmodel
          weights <- m.gam.full$model$`(weights)`
          dat <- cbind(dat, weights)
          m.gam.0 <- gam(f.form.0, data = dat, sp = l.lambd.0, weights = weights, family = v.family[[3]])

          p.sum$p[ term.ii ] <- anova(m.gam.0, m.gam.full, test = "F")[[ "Pr(>F)"]][2]

          # No test possible (?)
          # use mean of drop.1-anova
          if( is.na( p.sum$p[ term.ii ]) ){
            s.index <- grepl( p.sum$name.output[ term.ii ], attr( anova(m.gam.full)$chi.sq, "name" ), fixed = TRUE )
            p.sum$p[ term.ii ] <- mean( anova(m.gam.full)$s.pv[ s.index ] )
          }

        } else {

          # If interaction is not a factor
          # Save mean p-value from Anova object
          s.index <- which( p.sum$name.output[ term.ii ] == names.spline.coef )
          p.sum$p[ term.ii ] <- anova(m.gam.full)$s.pv[ s.index ]

        }

      } else {

        # Smooth terms (no interaction at all)
        p.sum$name.output[ term.ii ] <-
          paste( strsplit( term, "\\,")[[1]][1], ")", sep = "")
        s.index <- which( p.sum$name.output[ term.ii ] == names.spline.coef )
        p.sum$p[ term.ii ] <- anova(m.gam.full)$s.pv[ s.index ]

      }
    }

    t.spat.label <- paste0("te(", paste(coords,collapse=",") ,")")
    # if spatial
    if( (grepl( "te\\(", term)+!grepl("by", term)) == 2 ) {
      p.sum$p[ term.ii ] <-
        summary( m.gam.full)$s.pv[ names.spline.coef == t.spat.label ]
      # Save output name for lambda deletion
      p.sum$name.output[ term.ii ] <- t.spat.label
    }

    # if spatial interaction
    if( (grepl( "te\\(", term)+grepl("by", term)) == 2 ) {

      name.by <- strsplit( term, "\\,|\\=")[[1]][4]
      # Create name as on output and save in object (for lambdas extraction)
      p.sum$name.output[ term.ii ] <- nam <- paste( t.spat.label, ":", name.by, sep = "")

      p.sum$p[ term.ii ] <-
        summary( m.gam.full)$s.pv[ names.spline.coef == nam ]
    }
  }

  # Remove term with highest p-Value from formula
  p.sum$name <- as.character( p.sum$name )
  l.red.original <- p.sum[ p.sum$p == max(p.sum$p), "name"]
  if( length(l.red.original) > 1){ l.red.original <- l.red.original[1] }
  l.red <- gsub("\\\n", "", l.red.original)

  f.terms <- terms(formula(m.gam.full))
  n.drop <- which(gsub(" ", "", attr(f.terms, "term.labels"), fixed = TRUE) == l.red)
  f.drop <- formula( drop.terms( terms(formula(m.gam.full)),
                                 dropx = n.drop, keep.response = TRUE) )


  # Remove lambda from list if smooth was removed
  l.lambd.drop <- m.gam.full$full.sp
  if( sum( grepl("s\\(", l.red) ) >=1 ){

    # Select output name
    out.name <- p.sum$name.output[ p.sum$name == l.red.original ]
    # if interaction make grep (levels are appended to the name)
    if( grepl( ":", out.name ) ) {
      l.lambd.drop <- m.gam.full$full.sp[ grep( out.name,
                                                names( m.gam.full$full.sp ), fixed = TRUE, invert = TRUE ) ]
    } else {
      l.lambd.drop <- m.gam.full$full.sp[ -match( out.name,
                                                  names( m.gam.full$full.sp ) ) ]

    }
  }

  # if spatial is removed, remove te()-lambda
  if( sum( grepl("te\\(", l.red) + !grepl( "by", l.red ) ) == 2  ){
    l.grep <- grep( gsub( "\\)", "\\\\)", gsub( "\\(", "\\\\(", paste0("te(", paste(coords,collapse=",") ,")") )),
                    names(m.gam.full$full.sp))
    l.grep2 <- grep( "\\):", names(m.gam.full$full.sp)[l.grep])
    if(length(l.grep2)>0) l.grep <- l.grep[-l.grep2]
    l.lambd.drop <- m.gam.full$full.sp[ -l.grep ]

  }
  # if spatial interaction, remove
  if( sum( grepl("te\\(", l.red) + grepl( "by", l.red ) ) == 2 ){

    out.name <- paste0(p.sum$name.output[ p.sum$name == l.red.original ], c(1,2))

    l.lambd.drop <- m.gam.full$full.sp[ -match( out.name,
                                                names( m.gam.full$full.sp )) ]
  }

  return( list( f.drop, if( length( l.lambd.drop ) == 0 ){ NULL } else {l.lambd.drop}, l.red.original ) )
}




## -------- Function to cross validate GAM --------


# main CV function
f.gam.cv <- function(gam.model, resp, t.sets, mult.gam = 0, v.family, n.cores,
                     return.dat = FALSE, se.fit = FALSE){

  f.gam <- formula( gam.model )
  lambdas <- gam.model$full.sp
  d.dat <- gam.model$model[, -match("(weights)", names(gam.model$model) ), drop = FALSE]
  weights <- gam.model$model$`(weights)`
  d.dat <- cbind(d.dat, weights)

  t.ids <- c()
  for(ii in 1:dim(t.sets)[1]){
    t.id <- which(t.sets[ii, 1:ncol(t.sets)] == 0)
    t.ids <- c(t.ids, t.id )
  }

  my.cv.pred <- function(ii) {
    d.dat <- gam.model$model[, -match("(weights)", names(gam.model$model) ), drop = FALSE]
    weights <- gam.model$model$`(weights)`
    d.dat <- cbind(d.dat, weights)
    d.dat.cv <- d.dat[ which(t.ids != ii), ]

    # Calc model on training set
    m.gam.cv <- gam( f.gam, data = d.dat.cv, sp = lambdas, weights = weights, family = v.family[[3]] )

    # Predict on test set
    t.p <- predict( m.gam.cv, d.dat[ which(t.ids == ii),  ], type = "response")
    if( se.fit ){
      t.se <- predict( m.gam.cv, d.dat[ which(t.ids == ii),  ],
                       type = "response", se.fit = TRUE)$se.fit
      t.p <- data.frame( pred = t.p,  pred.se = t.se)}
    return(t.p)
  }

  l.pred <- mclapply(1:ncol(t.sets), my.cv.pred, mc.preschedule=FALSE, mc.set.seed=FALSE,
                     mc.cores = n.cores)

  t.val <- f.prepare.validation.cv(dat = d.dat, pred = l.pred, resp = resp, t.sets = t.sets,
                                   gam.fit = mult.gam, v.family = v.family, se.fit = se.fit)

  stats <- round( f.validate.predictions( t.val[[1]], t.val[[2]], t.val[[3]]), 5)

  if( se.fit ){
    return( list( stats, data.frame( data = t.val[[1]], pred = t.val[[2]],
                                     pred.se = t.val[[4]])) )
  } else if( return.dat ){
    return( list( stats, data.frame( data = t.val[[1]], pred = t.val[[2]]) ) )
  } else {

    return( stats  )
  }

}




## ------ Function to do backward selection ---------

f.gam.backward <- function( gam.model.full, t.sets, resp, v.family, n.cores, coords){

  # extract from gam.model
  dat <- gam.model.full$model


  # Count number of terms (minus intercept)
  # List of Elements in Formula
  l.terms <- strsplit( paste( formula(gam.model.full)[3]), "\\+")[[1]]
  l.terms <- l.terms[ !grepl( "\\\n ", l.terms) ]

  # remove intercept
  n.terms <- length(  l.terms[ !l.terms %in% c(" ", "")]  )-1

  # Add first
  model.drop <- list( gam.model.full )
  cv.drop <- list( f.gam.cv(gam.model=gam.model.full, resp = resp, t.sets = t.sets,
                            v.family = v.family, n.cores = n.cores) )
  l.drop <- c(" ")
  l.rmse <- cv.drop[[1]][[ v.family[[4]] ]]

  for( ii in 2:n.terms ){

    # find term to drop
    form.drop <- f.gam.find.drop(m.gam.full=model.drop[[ii-1]], v.family = v.family, coords = coords)
    # save name of droped for results
    l.drop <- c( l.drop, form.drop[[3]])

    # calc gam model
    weights <- gam.model.full$model$`(weights)`
    dat <- cbind(dat, weights)
    if( length(form.drop[[2]]) == 0){ snu.sp <- NULL}else{ snu.sp <- form.drop[[2]]}
    model.drop[[ii]] <- gam(form.drop[[1]], sp = snu.sp, data = dat,
                            weights = weights, family = v.family[[3]])

    # CV
    cv.drop[[ii]] <-  f.gam.cv(model.drop[[ii]],  resp = resp, t.sets = t.sets,
                               v.family = v.family, n.cores = n.cores)
    l.rmse <- c(l.rmse, cv.drop[[ii]][[ v.family[[4]] ]])
  }

  path.backward <- data.frame( name.drop = l.drop, RMSE = l.rmse)

  # Chose not minimum, but a bit more reduced model
  # If RMSE rounded to 3 digits is not different, take this model
  #n.minimum <- which( path.backward$RMSE == min(path.backward$RMSE) )
  n.optimum <- max( which( round( path.backward$RMSE, 3) ==
                             round( min(path.backward$RMSE), 3 ) ) )

  return( list( path.backward, model.drop[[n.optimum]], n.optimum)  )
}




## ----- Function to aggregate one single factor -----

f.gam.aggregate.one <- function( gam.model, term, inter = NULL, t.sets, verbose = 2 , resp,
                                 n.cores, v.family){

  # Rename term to be aggregated in gam.model with suffix
  # - "ag" (categ. main effect)
  # - "ag:-last letters of interaction" (categ. interaction terms)

  if( is.null( inter ) ) {
    term.suffix <- paste( term, "ag", sep = "")

    gam.model$model[, term.suffix ] <- gam.model$model[, term]

    form.suffix <- update( formula( gam.model), paste( "~. + ", term.suffix, "- ", term )  )
    weights <- gam.model$model$`(weights)`
    dat <- cbind( gam.model$model, weights)
    gam.model <- gam(form.suffix,
                     sp = gam.model$full.sp, data = dat, weights = weights, family = v.family[[3]] )
    term <- term.suffix
  }

  # List counter and results objects
  i.list <- 2
  model.aggreg <- list( gam.model )
  cv.aggreg <- list( f.gam.cv(model.aggreg[[1]],  resp = resp, t.sets = t.sets,
                              v.family = v.family, n.cores = n.cores) )
  l.aggreg <- c( " ")
  l.rmse <- cv.aggreg[[1]][[  v.family[[4]] ]]

  repeat{ # until all levels are aggregated

    dat <-  model.aggreg[[i.list-1]]$model

    sep.i <- ifelse( is.null(inter), "", ":" )
    term.inter <- paste(term, inter, sep = sep.i)

    # Save partial residuals
    t.part <- predict( model.aggreg[[i.list-1]], type="terms" ) ## get term estimates
    n.term <- which( dimnames( t.part)[[2]] == term.inter )
    # swap names if empty
    if( is.na(n.term[1]) ) {
      n.term <- which( dimnames( t.part)[[2]] == paste(inter, term, sep = ":") )}

    t.resid <- residuals(model.aggreg[[i.list-1]], type = "working") +
      t.part[, n.term ]

    # Put levels in list
    t.fact <- model.aggreg[[i.list-1]]$model[, term]
    l.lev <- levels(t.fact)

    # Build up matrix for distance measure pairwise t-Tests
    m.t <- matrix( data = NA, nrow =length(l.lev),
                   ncol = length(l.lev), dimnames = list(l.lev, l.lev))

    if( is.null(inter) ){

      # if not an interaction: pairwise t-test
      for( ii in 1:(length(l.lev)-1) ) {
        for( ss in (ii+1):length(l.lev) ) {

          xx <- t.resid[ t.fact == l.lev[ii] ]
          yy <- t.resid[ t.fact == l.lev[ss] ]

          m.t[ ii, ss] <- round( t.test(xx, yy)$p.value, 10 )
        }
      }
    } else {

      # if interaction with :factor, then pairwise linear model with
      # t-test of interaction
      for( ii in 1:(length(l.lev)-1) ) {
        for( ss in (ii+1):length(l.lev) ) {

          # save partials in new data frame, reduce data.frame to 2 levels
          d.sub <- dat
          d.sub$t.resid <- t.resid
          d.sub <- d.sub[ d.sub[, term ] %in% c(l.lev[ii], l.lev[ss]),  ]

          ff.lm <- paste( "t.resid ~ ", term, "+", inter, "+", term, ":", inter, sep = "")

          # if linear dependend
          if( is.na(lm( as.formula( ff.lm), data = d.sub,
                        weights = d.sub[ , "(weights)"] )$coefficients[ 4] ) ){
            m.t[ ii, ss] <- 1
          } else{

            m.t[ ii, ss] <- round( summary(
              lm( as.formula( ff.lm), data = d.sub, weights = d.sub[ , "(weights)"] ) )$coefficients[ 4, 4], 10)
          }
        }
      }
    }

    # read levels with largest t-Test
    l.ag1 <- l.lev[ which( m.t == max( m.t, na.rm = TRUE), arr.ind = TRUE )[1, "row"] ]
    l.ag2 <- l.lev[ which( m.t == max( m.t, na.rm = TRUE), arr.ind = TRUE )[1, "col"] ]

    # Aggregate levels in data frame
    name.aggreg <- paste( l.ag1, l.ag2, sep = "--")
    levels(t.fact)[ match(c(l.ag1, l.ag2),levels(t.fact))] <- name.aggreg
    dat[, term ] <- t.fact

    # Compute and CV of aggregation step
    weights <- gam.model$model$`(weights)`
    dat <- cbind(dat, weights)
    model.aggreg[[i.list]] <- gam(formula(gam.model),
                                  sp = gam.model$full.sp, data = dat,
                                  weights = weights,  family = v.family[[3]] )
    cv.aggreg[[i.list]] <- f.gam.cv( model.aggreg[[i.list]],  resp = resp, t.sets,
                                     v.family = v.family, n.cores = n.cores)

    # Collect results
    l.rmse <- c( l.rmse, cv.aggreg[[i.list]][[ v.family[[4]] ]])
    l.aggreg <- c( l.aggreg, name.aggreg)
    if( verbose > 2) { print( data.frame( l.aggreg, l.rmse )[ i.list,  ] ) }

    i.list <- i.list + 1
    if( length( levels(t.fact)) == 2 ) { break } # stop if only two levels left
  }

  # Compute minimum RMSE
  path.aggreg <- data.frame( name.aggr = l.aggreg, RMSE = l.rmse )
  if( verbose > 2){ print(path.aggreg) }

  # Compute optimum, tendency towards aggregation, here 4 digits only
  #n.min <- which( path.aggreg$RMSE == min(path.aggreg$RMSE) )
  n.opt <- max( which( round( path.aggreg$RMSE, 4) ==
                         round( min(path.aggreg$RMSE), 4 ) ) )

  return( list( path.aggreg, levels( model.aggreg[[n.opt]]$model[, term]),
                model.aggreg[[n.opt]])  )

}



## ---- Function to do factor aggregation for all factors -------

f.gam.aggregation.all <- function( gam.model, t.sets, verbose, resp, n.cores, v.family ){

  # get plain factors

  # List of Elements in Formula
  l.terms <- strsplit( paste( formula(gam.model)[3]), "\\+")[[1]]
  # remove spaces
  l.terms <- gsub(" ", "", l.terms, fixed = TRUE)
  # remove intercept
  l.terms <- l.terms[ l.terms != "1"]

  # read original data from object$
  m.agg <- gam.model
  d.agg <- gam.model$model
  n.fact <- 1; path <- list()

  for( term.ii in 1:length(l.terms) ){
    term <- l.terms[ term.ii ]


    #if( verbose == 2){ cat( "\n compute ", term, "\n \n")}

    # Linear main terms (are in the dataframe)
    if( term %in% names(d.agg)){
      # if factor
      if( is.factor( d.agg[, term]) && length(levels( d.agg[, term])) > 2 ){
        # compute best aggregation, give back new data.frame and new model
        aggr.obj <- f.gam.aggregate.one( m.agg, term, t.sets = t.sets, verbose = verbose,
                                         resp = resp, n.cores = n.cores, v.family = v.family)
        # Collect path and final levels
        path[[n.fact]] <- list( aggr.obj[[1]], aggr.obj[[2]])
        names( path )[n.fact] <- term
        m.agg <- aggr.obj[[3]]
        d.agg <- m.agg$model

        n.fact <- n.fact + 1
      }
    }

    # if linear interaction
    if( grepl(":", term )) {

      # Check if one is a factor  (otherwise do nothing)
      elem.1 <- is.factor( d.agg[ , str1 <- strsplit(term, split = ":")[[1]][1] ] )
      elem.2 <- is.factor( d.agg[ , str2 <- strsplit(term, split = ":")[[1]][2] ] )
      elem.1.l <- length( table( d.agg[, str1] ) ) > 2
      elem.2.l <- length( table( d.agg[, str2] ) ) > 2

      if( (elem.1 && elem.1.l) | ( elem.2 && elem.2.l ) ){
        # compute best aggregation, give back new data.frame and new model
        l.str <- c( str1, str2)
        aggr.obj <- f.gam.aggregate.one( m.agg, l.str[c(elem.1,elem.2)],
                                         l.str[!c(elem.1,elem.2)],
                                         t.sets = t.sets, verbose = verbose,
                                         resp = resp, n.cores = n.cores, v.family = v.family)
        # Collect path and final levels
        path[[n.fact]] <- list( aggr.obj[[1]], aggr.obj[[2]])
        names( path )[n.fact] <- term
        m.agg <- aggr.obj[[3]]
        d.agg <- m.agg$model

        n.fact <- n.fact + 1
      }
    }
  }

  return( list( path, m.agg) )

}



## --- Externe Validierung -----

f.gam.val.extern <- function( gam.model, resp, val.data, cal.data, mult.gam = 0, v.family, se.fit = FALSE, coords) {

  val.data.original <- val.data

  # Transformed calibration data
  daten <- gam.model$model

  ## Transform the data like the original data set
  # Scale numerical data to [0;1], then center (- mean), but not the response
  # Use means / max-min as in calibration data set!

  l.sub.fact <- l.sub.fact.1 <- unlist( lapply(val.data, is.factor) )
  l.sub.fact[ names(l.sub.fact) %in% c( resp, coords) ] <- TRUE


  fun.center.v <- function(col, vali = val.data, kali = cal.data, daten = gam.model$model ){

    x.c <- ( vali[, col] - min(kali[, col]) ) /
      ( max(kali[, col]) - min(kali[, col]) )

    # Mean from [0,1]-column from calibration-data set
    x.kali.m <- ( kali[, col] - min(kali[, col]) ) /
      ( max(kali[, col]) - min(kali[, col]))

    # center with this mean
    x.sc <- scale(x.c, center = mean( x.kali.m ), scale = FALSE)[, 1]
    return(x.sc)
  }

  val.data[, !l.sub.fact]  <- data.frame(
    lapply( names(val.data[, !l.sub.fact, drop = FALSE]), FUN = fun.center.v) )


  # Center coordinates with mean from calibration data
  if( !is.null(coords)){
    coord.means <-  colMeans( cal.data[ , coords])
    val.data[, coords] <- data.frame(
      scale( val.data[, coords], center = coord.means, scale = FALSE))
  }

  # Add intercept
  val.data$int <- rep(1, nrow(val.data) )

  # Add renamed factors for main effect
  val.data[ , paste( names( l.sub.fact.1 )[ l.sub.fact.1 ], "ag", sep = "") ] <-
    val.data[ , paste( names( l.sub.fact.1 )[ l.sub.fact.1 ] ) ]


  ##  Aggregate factors accordingly
  l.factors <- names( daten )[ unlist( lapply( daten, is.factor ) ) ]

  for( fact in l.factors ){

    if( fact %in% names(val.data) ){
      ii.val.fact <- which( names(val.data) == fact )
    } else {
      ii.val.fact <- which( names(val.data) == substr( fact, 0, nchar(fact)-5) )
      # Add, if renamed factor
      val.data[, fact ] <- val.data.original[ , ii.val.fact ]
    }

    n.lev.cal <- length( lev.cal <- levels( daten[, fact] ) )
    n.lev.val <- length( lev.val <- levels( val.data[, fact ]) )

    # if less levels in calibration set do the aggregation
    if( n.lev.cal < n.lev.val ){
      # Aggregate each level stepwise
      for( ii in 1:n.lev.val ){
        levels(val.data[, fact ])[ match( c( strsplit( lev.cal[ii], "--" )[[1]] ),
                                          levels( val.data[, fact ]) )]  <- lev.cal[ii]
      }
    }
  }

  resp <- val.data[ , resp ]
  val.dat.p <- val.data

  pred <- predict( gam.model, newdata = val.dat.p, type = "response")
  if( se.fit ){
    se.pred <- as.numeric( predict(gam.model, newdata = val.dat.p, type = "response", se.fit = TRUE)$se.fit )
  }

  t.val <- f.prepare.validation.ext(dat = resp, pred = pred, gam.fit = mult.gam, v.family = v.family)
  stats <- round( f.validate.predictions(t.val[[1]], t.val[[2]], t.val[[3]]), 5)

  d.pr <- data.frame( dat = resp, pred)
  names(d.pr) <- if( length(dim(pred)) == 1){ c( "dat", "pred")}else{ c("dat", paste0("pred.", 1:ncol(pred))) }

  if( se.fit ){ d.pr <- cbind( d.pr, se.pred)}

  return( list( Externe_Validierung = stats,
                Predictions = d.pr
  ) )
}



# ---- FINAL Selection Function -------
# to concatenate full model selection process

geoGAM <- function(
  # input
  response,
  covariates = names(data)[!(names(data) %in% c(response,coords))],
  data,
  coords = NULL, #c("x", "y"),
  weights = rep(1, nrow(data)),
  offset = TRUE, # 0 = no offset, 1 = group lasso offset
  max.stop = 300,
  non.stationary = FALSE,
  sets = NULL,
  seed = NULL,
  validation.data = NULL,
  verbose = 0, # print stuff, if = 1, a lot if = 2
  cores = min(detectCores(),10)
){

  # rename arguments troughout function
  # use real logical for offset and non.stat

  if( !is.null(coords) ){
    data$x <- data[, coords[1]]
    data$y <- data[, coords[2]]
  }

  resp <- response
  l.covar <- covars <- covariates
  calib.data <- data
  w.weights <- weights
  off.choice <- ifelse( offset, 1, 0)
  val.data <- validation.data
  non.stat <- ifelse( non.stationary, 1, 0)
  t.sets <- sets
  n.cores <- cores

  verbose <- ifelse( verbose > 0, 2, verbose)

  if( is.null(t.sets) ) t.sets <- mboost::cv(rep(1, nrow(data)), type = "kfold")
  if( is.null(dim(t.sets)) ){ # change to matrix format
    mx <- matrix(1, ncol = max(t.sets), nrow = length(t.sets) )
    for(xx in 1:length(t.sets)){
      mx[xx, t.sets[xx]] <- 0 }
    t.sets <- mx
  }

  # if windows set cores to 1
  if ( get_os() == "windows") n.cores <- 1


  dat <- calib.data
  # -------


  # tests on calibraion data
  if( !is.null(covariates)){
    if( !sum(c(response,covariates) %in% names(data)) == length(c(response,covariates)) ){
      stop("Not all covariates in data frame. Please check.")
    }
    if( !is.null(validation.data) & !sum(c(response,covariates) %in% names(validation.data)) == length(c(response,covariates)) ){
      warning("Not all covariates in validation data. External validation ignored.")
      validation.data <- NULL
    }
  } else {
    if( !(response %in% names(data)) ){
      stop("Response not in data frame.")
    }
    if( !(response %in% names(validation.data)) ){
      warning("Response not in validation data. External validation ignored.")
      validation.data <- NULL
    }
    if( !( names(data)[!(names(data) %in% response)] %in% names(validation.data)) ){
      warning("Covariates not complete in validation data. External validation ignored.")
      validation.data <- NULL
    }
  }


  # Check for response type, choose family
  if( is.ordered(calib.data[, resp])){
    v.family <- list( PropOdds(), "PropOdds", ocat(R=nlevels(dat[, resp])), "rps")
  } else if( nlevels(calib.data[, resp]) == 2 | length(table(calib.data[, resp]) ) == 2 )  {
    v.family <- list( Binomial(), "Binomial", binomial(), "bs")
  } else if( nlevels(calib.data[, resp]) > 2){
    v.family <- list(Multinomial(), "Multinomial", binomial(), "bss")

    cat("Your response is nominal.")
    stop("Model selection for multinomial models not yet implemented.")

  } else {
    v.family <- list(Gaussian(), "Gaussian", gaussian(),"rmse")
  }

  if( v.family[[2]] == "Binomial" & !is.factor(calib.data[, response]) ){
    calib.data[, response] <- factor( calib.data[, response],
                                      levels = sort(unique(calib.data[, response])),
                                      labels = c("absent", "present") ) }



  if( verbose > 0){ cat("Do model selection with family:  ", v.family[[2]])
    cat("\nNumber of observations for fitting: ", nrow(calib.data))
    cat("\nNumber of covariates: ", ncol(calib.data)-1) }

  if( verbose >= 2){ cat( "\n\n1. Transform data ... \n") }
  data.sel <- calib.data[ , c( resp, coords, l.covar) ]

  if( verbose > 2){
    cat("\n\nSummary of response: \n ")
    if( v.family[[2]] == "Gaussian" ){
      print( summary(data.sel[ , resp]) )
    } else{
      print( table(data.sel[ , resp] ))
    }
  }

  val.data.full <- val.data
  if( !is.null(val.data) ){
    val.data.sel <- val.data[ , c( resp, coords, l.covar) ]

  } else {
    val.data.sel <- data.sel
  }
  if( verbose > 2){ cat( "untransformed data: \n \n"); print(str(data.sel)) }

  data.trans <- f.transform.data( resp, data.sel, val.data.sel, coords = coords )
  if( !is.null(val.data)){ val.data <- data.trans[[2]]}
  data.trans <- data.trans[[1]]
  if( verbose > 2){ cat( "\n transformed data: \n \n"); print(str(data.trans)) }

  l.fact <- names(data.trans)[ unlist( lapply( data.trans, is.factor) ) == TRUE ]
  l.fact <- grep( paste0("^", resp, "$"), l.fact, invert = TRUE, value = TRUE)

  if(length(l.fact) == 0){ off.choice <- 3 }

  if( off.choice == 1 ){
    if( verbose >= 2){ cat( "\n2. Find factors as offset ... \n ") }

    # choose factors with cross validation
    tt <- f.choose.offset(resp = resp, data.trans, t.sets = t.sets,
                          w.weights = w.weights, v.family = v.family, verbose = verbose)

    if( verbose >= 1){ cat("\n\nChosen factors: \n" ); print(tt[[1]]) }

  } else {
    # Use no offset
    tt <- list()
    tt[[1]] <- tt[[2]] <- tt[[4]] <- ""
    tt[[3]] <- rep( 0, nrow(data.trans) )
  }


  if( verbose >=1){
    cat( "\nOffset used for boosting: \n ")
    print( summary( tt[[3]] ) ) }


  # Create homogenous cross validation sets
  if( !is.null(seed)){
    set.seed(seed)
    t.sets <- mboost::cv(rep(1, nrow(data.sel)), type = "kfold")
  }

  if( verbose >= 1){ cat( "\n3. Do boosting (time consuming) ... \n ") }

  t.boost <- f.boost( resp = resp, data.trans, p.offset.lm = tt[[3]], v.family = v.family,
                      val.data = val.data, f.fact.off = tt[[2]], t.sets = t.sets,
                      gam.selection = 1, non.stat = non.stat, max.stop = max.stop, coords = coords )


  if( verbose >= 1){ cat(" \n chosen baselearners main effects: \n ");
    print( t.boost[[1]]);
    #                      cat(" \n chosen baselearners interactions: \n");
    #                      print( t.boost[[2]] )
    if( !is.null(val.data)){
      cat( "\n cross validation RMSE: \n")
      print( t.boost[[4]][1])
    }
  }

  if( length(t.boost[[1]] > 3)){
    #cat( "\n Weights used in boosting : \n")
    #print( summary( extract( t.boost[[5]], what = "weights")) )

    if( verbose > 0) cat( "\n4. Cross validate tuning parameters for baselearner selection: \n")
    # find good baselearner set to continue

    # prop odds gam needs numeric factor
    if( v.family[[2]] == "PropOdds" ){ data.trans[ , resp] <- as.numeric( data.trans[, resp] ) }

    silent.Note <- capture.output(
      m.models.cv <- f.tuning.param(t.boost[[5]][mstop(t.boost[[6]])], data.trans,
                                    non.stat = non.stat, resp = resp, t.sets = t.sets ) )


    # create a function from the folowing.....
    # print(summary(w.weights))

    l.cv.res <- list()
    # loop through possible models
    for( r.row in 1:nrow(m.models.cv) ){
      m.one <- m.models.cv[ r.row, ]

      # factors in model
      l.f <- grepl( "\\(", names( m.one ) )
      l.fn <- names( m.one )[ m.one == 1 & !l.f ]

      l.baselearner <- names( m.one )[ m.one == 1 & l.f ]
      l.linear <- names( m.one )[ m.one == 2 ]




      l.baselearn.int <- l.baselearner[ grep( "by", l.baselearner, invert = FALSE) ]
      l.baselearn.haupt <- l.baselearner[ grep( "by", l.baselearner, invert = TRUE) ]
      weights <- w.weights

      silent.Note <- capture.output(
        m.gam.formula <- f.transform.formula(l.baselearn.haupt= l.baselearn.haupt,
                                             l.baselearn.interact= l.baselearn.int,
                                             l.offset = c( tt[[1]], l.fn ),
                                             daten = data.trans,
                                             resp = resp,
                                             l.linear = l.linear,
                                             w.weights = w.weights, coords = coords), type = "message"
      )

      dat <- cbind( m.gam.formula[[3]], weights)
      m.gam.tune <- gam( m.gam.formula[[1]], data = dat,
                         sp = m.gam.formula[[2]], weights = weights, family = v.family[[3]])

      cv.stats <- f.gam.cv(m.gam.tune, t.sets = t.sets, resp = resp, v.family = v.family,
                           n.cores = n.cores)

      #     ext.stats <- f.gam.val.extern( m.gam.tune,
      #                                    resp = resp, val.data = val.data.sel,
      #                                    cal.data = data.sel, v.family = v.family )[[1]]
      l.cv.res[[r.row]] <- list( cv.stat = cv.stats )#, ext.stat = c(0) )
    }

    # ....... return list of baselearner

    if( verbose > 0) cat("\nNumber of cross validated GAM: ", nrow(m.models.cv) )

    l.cv <- l.ext <- c()
    for( ii  in 1:length( l.cv.res ) ){
      l.cv <- c(l.cv, l.cv.res[[ii]][[1]][[ v.family[[4]] ]])
      # l.ext <- c( l.ext, l.cv.res[[ii]][[2]][[ v.family[[4]] ]])
    }

    # search for best cv model
    # take minimum RMSE and least complex model, if they are similar
    if( ncol(m.models.cv) > 2){
      m.models.cv$numb <- apply( m.models.cv[, -match("modell", names(m.models.cv))], 1,
                                 function(x){ sum( x == 1 )} )
    } else {
      m.models.cv$numb <- m.models.cv[, -match("modell", names(m.models.cv))]
    }

    l.min <- which( round(l.cv,2) == min( round(l.cv, 2) ) )
    l.complex <- which( m.models.cv$numb == min( m.models.cv$numb[ l.min ] ) )

    l.choose <- l.complex[ l.complex %in% l.min ]

    if(length(l.choose) > 1){ l.choose <- l.choose[ l.cv[ l.choose ] == min( l.cv[ l.choose ] ) ] }

    if( verbose > 0) cat( "\nCV GAM Models, line chosen (see csv):  ", l.choose, "\n")

    dpr <- cbind( m.models.cv[ c("modell", "numb")], l.cv )
    rownames(dpr) <- NULL
    if( verbose >= 2) print(dpr)

    m.one.fin <- m.models.cv[ l.choose, ]
    l.baselearner.fin <- names( m.one.fin[-ncol(m.one.fin)] )[ m.one.fin[-ncol(m.one.fin)] == 1 ]
    l.linear.fin <- names( m.one.fin )[ m.one.fin == 2 ]

    if( verbose > 0) {cat(" \n\nBaselearner after Tuning : \n")
      print( l.baselearner.fin) }

    l.basel.h <- l.baselearner.fin[ grep( "by", l.baselearner.fin, invert = TRUE) ]
    l.basel.h <- l.basel.h[ grepl("\\(", l.basel.h) ]
    l.basel.i <- l.baselearner.fin[ grep( "by", l.baselearner.fin)]

    l.basel.fact <- l.baselearner.fin[ !grepl("\\(", l.baselearner.fin) ]
    if( length(l.basel.fact) == 0){ l.fa <- tt[[1]] } else { l.fa <- c(tt[[1]], l.basel.fact)}


  } else{
    l.baselearner.fin <- t.boost[[1]]

    l.basel.h <- l.baselearner.fin[ grep( "by", l.baselearner.fin, invert = TRUE) ]
    l.basel.h <- l.basel.h[ grepl("\\(", l.basel.h) ]
    l.basel.i <- l.baselearner.fin[ grep( "by", l.baselearner.fin)]

    l.basel.fact <- l.baselearner.fin[ !grepl("\\(", l.baselearner.fin) ]
    if( length(l.basel.fact) == 0){ l.fa <- tt[[1]] } else { l.fa <- c(tt[[1]], l.basel.fact)}

  }


  if( verbose >= 2){ cat("\n5. Transform formula to mgcv style... \n ") }
  gam.formula <- f.transform.formula(l.baselearn.haupt = l.basel.h,
                                     l.baselearn.interact= l.basel.i,
                                     l.offset = l.fa,
                                     daten = data.trans,
                                     resp = resp,
                                     l.linear = l.linear.fin,
                                     w.weights = w.weights,
                                     coords = coords
  )
  if( verbose > 2){ print( gam.formula[[1]] ) }


  if( verbose >= 1){ cat("\n6. Compute full GAM ... \n ") }



  weights <- w.weights
  dat.gam <- cbind(gam.formula[[3]], weights)
  gam.full <- gam( gam.formula[[1]], data = dat.gam,
                   sp = gam.formula[[2]], weights = weights, family = v.family[[3]])
  if( verbose >= 1){ print( formula(gam.full) ) }
  #if( verbose >= 1){ print(summary(gam.full)) }

  #   if( verbose > 0) cat( "\n \n--- Cross validation of full GAM model --- \n ")
  #
  #   cv.gam <- f.gam.cv( gam.full, t.sets = t.sets,  resp = resp,
  #                       v.family = v.family, n.cores = n.cores)
  #   if( verbose > 0) print( round(cv.gam, 5) )


  # if only one term left in formula, don't do backward selection
  l.terms <- strsplit( paste( formula(gam.full)[3]), "\\+")[[1]]
  l.terms <- gsub(" ", "", l.terms, fixed = TRUE); l.terms <- gsub( "\n", "", l.terms, fixed = TRUE)
  l.terms <- l.terms[ !l.terms %in% c("1", "")]


  if( verbose >= 1){ cat("\n \n7. Compute GAM backward selection ... \n ") }
  if( length(l.terms) > 1){
    gam.back <- f.gam.backward(gam.model.full=gam.full, t.sets = t.sets, resp = resp,
                               v.family = v.family, n.cores = n.cores, coords = coords )
    if( verbose >= 1) { print( gam.back[[1]] );
      cat( "\n Optimal model: line ", gam.back[[3]]) }
  } else {
    gam.back <- list()
    gam.back[[3]] <- 1
    gam.back[[2]] <- gam.full
  }


  # List of Elements in Formula
  l.terms <- strsplit( paste( formula(gam.back[[2]])[3]), "\\+")[[1]]
  l.terms <- gsub(" ", "", l.terms, fixed = TRUE)
  l.terms <- l.terms[ l.terms != "" ]
  l.terms <- l.terms[ l.terms != "1"]
  l.terms <- l.terms[ !grepl("^s\\(|^te\\(", l.terms) ]

  if( length(l.terms) > 0 ){
    if( verbose >= 1){ cat("\n \n8. Compute aggregation of factor levels  ... \n ") }
    gam.aggr <- f.gam.aggregation.all( gam.back[[2]], t.sets = t.sets,
                                       verbose = verbose, resp = resp, n.cores = n.cores,
                                       v.family = v.family)
    #if( verbose >= 1) { cat( "\n path of aggregation \n"); print( gam.aggr[[1]]) }

  } else{

    gam.aggr <- list()
    gam.aggr[[1]] <- c()
    gam.aggr[[2]] <- gam.back[[2]]
  }

  if( grepl("^log\\.|^sqrt.|logit\\.", resp) ){
    fin.cv.f <- f.gam.cv(  gam.aggr[[2]],  resp = resp, t.sets = t.sets,
                           v.family = v.family, n.cores = n.cores, return.dat = TRUE,
                           se.fit = grepl("^log.|^sqrt.", resp) )
  }

  fin.cv <- f.gam.cv(  gam.aggr[[2]],  resp = resp, t.sets = t.sets,
                       v.family = v.family, n.cores = n.cores, return.dat = TRUE)

  if( verbose >= 1){ cat("\n \n9. Final Model ... \n ")
    print( summary( gam.aggr[[2]]))
    cat( "\n \n--- Cross validation of final reduced GAM model --- \n ")
    print( round(  fin.cv[[1]], 5)  )

  }


#   # Backtransformed CV
  if( grepl("^log\\.", resp) ){

    fin.cv.f[[2]]$data.transf <- fin.cv.f[[2]]$data
    fin.cv.f[[2]]$pred.transf <- fin.cv.f[[2]]$pred

    fin.cv.f[[2]]$pred <- as.numeric( exp(fin.cv.f[[2]]$pred.transf -
                                                      fin.cv.f[[2]]$pred.se*0.5) )
    fin.cv.f[[2]]$data <- calib.data[ , gsub("^log\\.", "", resp)]

    fin.cv <- fin.cv.f
  }
  if( grepl("^sqrt.", resp) ){

    fin.cv.f[[2]]$data.transf <- fin.cv.f[[2]]$data
    fin.cv.f[[2]]$pred.transf <- fin.cv.f[[2]]$pred

    fin.cv.f[[2]]$pred <- as.numeric( fin.cv.f[[2]]$pred.transf^2 - fin.cv.f[[2]]$pred.se )
    fin.cv.f[[2]]$data <- calib.data[ , gsub("^sqrt.", "", resp)]

    fin.cv <- fin.cv.f

  }
#   if( grepl( "^logit\\.", resp)) {
#     fin.cv.f[[2]]$pred.backtrnsf <- as.numeric(( exp(fin.cv.f[[2]]$pred) /
#                                                    (1 + exp(fin.cv.f[[2]]$pred)) )*100)
#     fin.cv.f[[2]]$dat.backtrnsf <- calib.data[ , gsub("^logit\\.", "", resp)]
#   }


  #   # WARNING: this code is not generic for any dataset!!!
  #   l.c <- ifelse( sum(grepl("timeset", names(calib.data))), c("site_id_unique", "timeset"), "site_id_unique")
  #   fin.cv.f[[2]] <- cbind(calib.data[ ,l.c, drop = FALSE], fin.cv.f[[2]])
#
#   if( grepl("^log\\.|^sqrt.|logit\\.", resp) ){
#     cat("\n--- computed on backtransformed response --- \n")
#     print( round( c( f.validate.predictions(d <- fin.cv.f[[2]]$dat.backtrnsf,
#                                             p <- fin.cv.f[[2]]$pred.backtrnsf, "st")), 5 ))
#   }


  if( !is.null(val.data)){
    cat( "\n \n--- External validation of final reduced GAM model --- \n ")
    l.val <- f.gam.val.extern( gam.aggr[[2]],
                               resp = resp, val.data = val.data.sel,
                               cal.data = data.sel, v.family = v.family,
                               se.fit = grepl("^log.|^sqrt.", resp), coords = coords)

    if( grepl("^log\\.", resp) ){

      l.val[[2]]$dat.transf <-  l.val[[2]]$dat
      l.val[[2]]$pred.transf <-  l.val[[2]]$pred

      l.val[[2]]$pred <- as.numeric( exp(l.val[[2]]$pred.transf - l.val[[2]]$se.pred*0.5) )
      l.val[[2]]$dat <- val.data.full[ , gsub("^log\\.", "", resp)]
    }
    if( grepl("^sqrt\\.", resp) ){

      l.val[[2]]$dat.transf <-  l.val[[2]]$dat
      l.val[[2]]$pred.transf <-  l.val[[2]]$pred

      l.val[[2]]$pred <- as.numeric( l.val[[2]]$pred.transf^2 - l.val[[2]]$se.pred )
      l.val[[2]]$dat <- val.data.full[ , gsub("^sqrt\\.", "", resp)]
    }
#     if( grepl( "^logit\\.", resp)) {
#       l.val[[2]]$pred.backtrnsf <- as.numeric(( exp(l.val[[2]]$pred) / (1 + exp(l.val[[2]]$pred)) )*100)
#       l.val[[2]]$dat.backtrnsf <- val.data.full[ , gsub("^logit\\.", "", resp)]
#     }


    #     # WARNING: this code is not generic for any dataset!!!
    #     l.c <- ifelse( sum(grepl("timeset", names(val.data.sel))), c("site_id_unique", "timeset"), "site_id_unique")
    #     l.val[[2]] <- cbind(val.data.full[ ,l.c, drop = FALSE], l.val[[2]])

    print( l.val[[1]] )

    # if( grepl("^log\\.|^sqrt.|logit\\.", resp) ){
    #   cat("\n--- computed on backtransformed response --- \n")
    #   print( round( c( f.validate.predictions(d <- l.val[[2]]$dat.backtrnsf,
    #                                           p <- l.val[[2]]$pred.backtrnsf, "st")), 5 ))
    # }


  }

  if( verbose > 0) cat( " \n \n ---- computation finished ... ! \n")


  result <- list(
    offset.grplasso = tt[[4]],
    offset.factors = tt[[1]],
    gamboost = t.boost[[5]],
    gamboost.cv = t.boost[[6]],
    gamboost.mstop = t.boost[[10]][1],
    gamback.cv = list(rmse.list = l.cv.res, optimal = l.choose, gamback.cv.model = gam.full),
    gamback.backward = list( rmse.list = gam.back[[1]], optimal =  gam.back[[3]],
                             gamback.backward.model =  gam.back[[2]]),
    gamback.aggregation = gam.aggr[[1]],
    gam.final = gam.aggr[[2]],
    gam.final.cv = fin.cv[[2]],
    gam.final.extern = if( !is.null(val.data) ){ l.val[[2]] } else { NULL },
    data = data,
    validation.data = validation.data,
    parameters = list( response = response, covariates = covariates, family = v.family, weights = weights, off.choice = off.choice,
                       coords = coords, non.stationary = non.stationary, sets = sets,
                       seed = seed, max.stop = max.stop)
  )

  class(result) <- "geoGAM"
  return(result)
}



# Compute prediction statistics - Function from GEOROB by A. Papritz
# Angepasst, dass Gewichte in Berechnung RMSE reinkommen
f.validate.predictions <- function (data, pred, statistic = c("st", "bs", "bss", "rps"),
                                    cv.weights = rep(1, length(data)))
{
  statistic = match.arg(statistic)
  if( !statistic %in% c("cbs", "bss", "rps") ){
    t.sel <- complete.cases(data, pred)
    if (sum(t.sel) < length(t.sel))
      warnings("missing values encountered when validating predictions")
    data <- data[t.sel]
    pred <- pred[t.sel]
  }

  result <- switch(statistic, st = {
    error <- (data - pred)*cv.weights
    statistics <- c(me = mean(error),
                    mede = median(error),
                    rmse = sqrt(mean(error^2)),
                    made = mad(error, center = 0),
                    R2 = 1 - ( sum(error^2) / sum((data - mean(data))^2) ) )
  }, bs = {
    # brier score
    if( is.factor(data) ){ data <- as.numeric(data)-1 }

    error2 <- ((data - pred)^2)*cv.weights
    # ignorance score (Wilks-2011, S 341.)
    pred.ig <- pred
    pred.ig[ data == 0 ] <- 1 - pred[ data == 0 ]
    igscore <- -log( pred.ig )*cv.weights

    statistics <- c(bs = mean(error2), is = mean(igscore) )

  }, bss = {
    # sum of pairwise brier scores
    data.m <- t(sapply( as.numeric( data ),
                        FUN = function(x){ c(rep(0, x-1), 1,
                                             rep(0, max(as.numeric(data))-x)) } ))

    statistics <- c( bss = mean( sapply( 1:length(data),
                                         function(x) ( (pred[x,] - data.m[x,])^2 ) ) ) )


  }, rps =  {
    # ranked probabilty score
    # cumulative sums over predicted and observed
    ym  <- t(apply( pred, 1, cumsum))
    om <- t(sapply( max(as.numeric(data)) - as.numeric( data ) + 1,
                    FUN = function(x){ c(rep(0, max(as.numeric(data))-x), rep( 1, x)  ) } ))

    # ignorance score (Wilks-2011, S. 351)
    igscore <- -log( sapply(1:nrow(pred), FUN = function(x){ pred[x, as.numeric( data[x] ) ] } )
                     + 0.00000000000001)

    statistics <- c( rps = mean( sapply( 1:length(data),
                                         function(x) sum( (ym[x,] - om[x,])^2 ) ) ),
                     is = mean(igscore) )
  }

  )
  return(result)
}





# get windows back
# > quantregForest:::get_os
get_os <- function ()
{
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf["sysname"]
    if (os == "Darwin")
      os <- "osx"
  }
  else {
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
