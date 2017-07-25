## File: simulate_multivariate.r 
## Description: Functions to conduct multivariate simulation using 
##              different approaches (indep, indep on ICA, copula, BN).
## Date: April 2016
## By: Cugliari, Rolland


simulatemulti <- function(data, method, n, ...){
  if(missing(n)) n <- nrow(data)
  switch(method,
    indep    = simudataindep(data, n),
    indepPCA = simudataindepPCA(data, n, ...),
    copula   = simudatacopula(data, n, ...)
    #,    BN       = simudataBN(data, n)
  )
}

## Simulation from independent data 
simudataindep <- function(data, n_simu) {
  if(missing(n_simu)) n_simu <- nrow(data)
  
  data_sim <- matrix(runif(n_simu * ncol(data)), nrow = n_simu)

  for(j in 1:ncol(data)) {
    data_sim[, j] <- icdf(u = data_sim[, j], 
                          x = data[, j], 
                          n = n_simu)
  }
  
  colnames(data_sim) <- colnames(data)
  return(data_sim)
}

## Simulation from independent data 
simudataindepPCA <- function(data, n_simu) {
  if(missing(n_simu)) n_simu <- nrow(data)
  prfit <- prcomp(data, scale. = TRUE, retx = TRUE)

  x_sim <- simudataindep(data = prfit$x, n_simu)
  
  data_sim_sc <- x_sim %*% solve(prfit$rotation)
  data_sim <- sweep(data_sim_sc %*% diag(prfit$scale), 2, 
                    FUN = "+", prfit$center)
  
  colnames(data_sim) <- colnames(data)
  return(data_sim)
}

## Simulation from C- and D-vine copula models
simudatacopula <-function(data, n_simu, 
                    # Options for CopSelect
                    familyset     = NA,      # test all
                    type          = "DVine", # or CVine
                    selectioncrit = "AIC",   # or BIC
                    indeptest     = TRUE,    # test independence
                    plot          = FALSE    # should plot results
                    ){
  require(CDVine)
  
  if(missing(n_simu)) n_simu <- nrow(data)
  
  cdata <- apply(data, 2, fpx) # transform data into copula data
  
  ## Select a D-vine structure (each of the bivariate copulas 
  #  when the independence test is rejected.
  Vine_structure <-  CDVineCopSelect(cdata, 
                         familyset     = familyset, 
                         type          = type, 
                         selectioncrit = selectioncrit,
                         indeptest     = indeptest)

  ## Estimate a D-vine structure
  fit <-  CDVineMLE(cdata, 
                    family = Vine_structure$family, 
                    start  = Vine_structure$par, 
                    start2 = Vine_structure$par2, 
                    type   = type)
  
  cdata_sim <- CDVineSim(n_simu, 
                         family = Vine_structure$family, 
                         par   = fit$par, 
                         par2  = fit$par2, 
                         type  = type)
  
  if(plot) CDVineTreePlot(data   = cdata,
                          family = Vine_structure$family,
                          par    = fit$par,
                          par2   = fit$par2,
                          type   = "DVine",
                          edge.labels = c('family', "emptau"),
                          main = "")

  data_sim <- matrix(0, nrow = n_simu, ncol = ncol(data))
  for(j in 1:ncol(data)) {
    data_sim[, j] <- icdf(u = cdata_sim[, j], 
                          x = data[, j], 
                          n = n_simu)
  }
 
  colnames(data_sim) <- colnames(data)
  return(data_sim)
}


## Function: fpx
## Description : Obtain the ranks of a vector (A factor of
##               0.5 is used to avoid extremes values of ranks).
fpx <- function(x) (rank(x, ties.method = "max") - 0.5) / length(x)


## Function: icdf [1st paper version]
## Description : Generalized inverse of the empirical cumulative
##               function. (Should be checked!!!!)
# icdf <- function(u, x, n) {
#   xstar <- numeric(n)
#   for(i in seq_along(xstar)){
#     cand <- ((rank(x) - 0.5) / length(x) - u[i]) <= 0
#     xstar[i] <- ifelse(!any(cand), min(x), max(x[cand]))
#   }
#   return(xstar)  
# }

## Function: icdf
## Description : Generalized inverse of the empirical cumulative
##               function.
icdf <- function(u, x, n) {
  
  freq <- fpx(x)
  Fn   <- splinefun(x, freq, method = "monoH.FC")
  
  xstar <- numeric(n)
  
  for(i in seq_along(xstar)){
    xstar[i] <- uniroot(function(x) Fn(x) - u[i], 
                        range(x), extendInt = "upX", 
                        f.lower = - u[i], f.upper = 1 - u[i])$root
  }
  return(xstar)  
}

## Function: panel.contour
## Description : Nice representation of dependent data using pair plots
##               and contours of bivariate copula
panel.contour <- function(x, y, bw=2, size=100){
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(-3,3,-3,3), new=TRUE)
  BiCopMetaContour(x, y, bw, size, axes=FALSE)
}
