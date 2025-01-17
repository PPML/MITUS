#' Basic Functions Used In Package
#' This script includes several functions that are useful for data manipulation
#' throughout the package.

#'The first function takes the odds of a risk, then calculates the join
#'probability of that odds and and an inputted Odds Ratio that the user inputs, then
#'generates the corresponding risk.
#'@name ORAdd
#'@param val known risk for population of interest
#'@param OR odds ratio by which the known risk will be modified
#'@return new risk for population
ORAdd <- function(val,OR) {
  x <- (val/(1-val))*OR; x/(1+x)
}

#'This function calculates the logit of the value or parameter inputted.
#'@name lgt
#'@param x value of interest, often an odds
#'@return logit of input x
lgt <-  function(x) log(x/(1-x));

#'This function calculates the inverse logit of the value or parameter inputted.
#'@name invlgt
#'@param x value of interest, often an odds
#'@return inverse logit of input x
invlgt <- function(x) 1/(1+exp(-x))

#'This function is used to simulate a logistic increase between two years of interest
#'to a user-specified maximum. Useful for simulating slow increase/decrease of a
#'parameter over time.
#'@name LgtCurve
#'@param StYr Year to start the logit curve
#'@param Endyr Year to end the logit curve
#'@param EndVal maximum value that the logistic curve can reach
#'@return vector of values demonstrating a logistic increase from zero to EndVal

LgtCurve <- function(StYr,Endyr,EndVal) {
  z <- log(1/0.005-1)
  zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),
             z*(1+2*(2100-Endyr)/(Endyr-StYr)),
             by=(2*z)/(Endyr-StYr)/12)
  zz  <- as.numeric(EndVal)/(1+exp(-zz))
  if(StYr>1950) {
    zz[1:((StYr-1950)*12)] <- 0
  }
  zz  }

#'This function creates a smooth curve from a vector of values.
#'@name SmoCurve
#'@param vec vector of values that the user would like to have smoothed into a curve
#'@return expanded vector of values that determine the shape of the smoothing spline
SmoCurve <- function(vec) {
  jj <- predict(smooth.spline(x=1:length(vec),y=vec,spar=0.2),
                x=seq(1,length(vec),1/12))$y;
  jj[jj<0] <- 0 ; jj
}

#'This function creates a smooth curve from a vector of values.
#'@name SmoCurve_knots
#'@param vec vector of values that the user would like to have smoothed into a curve
#'@param nk number of knots
#'@return expanded vector of values that determine the shape of the smoothing spline
SmoCurve_knots <- function(vec,nk) {
  jj <- predict(smooth.spline(x=1:length(vec),y=vec,spar=0.2, nknots = nk),
                x=seq(1,length(vec),1/12))$y;
  jj[jj<0] <- 0 ; jj
}

#'This function creates a smooth curve from a vector of values.
#'@name SmoCurve_decade_month
#'@param vec vector of values that the user would like to have smoothed into a curve
#'@return expanded vector of values that determine the shape of the smoothing spline
SmoCurve_decade_month <- function(vec) {
  jj <- predict(smooth.spline(x=1:length(vec),y=vec,spar=0.2),
                x=seq(1,length(vec),1/120))$y;
  jj[jj<0] <- 0 ; jj
}
#'This function creates a smooth curve from a vector of values.
#'@name SmoCurve_decade_year
#'@param vec vector of values that the user would like to have smoothed into a curve
#'@return expanded vector of values that determine the shape of the smoothing spline
SmoCurve_decade_year <- function(vec) {
  jj <- predict(smooth.spline(x=1:length(vec),y=vec,spar=0.2),
                x=seq(1,length(vec),1/10))$y;
  jj[jj<0] <- 0 ; jj
}
#'This function is for the expit
#'@name expit
#'@param x number
#'@return expit of that number x
expit <- function(x) {
  exp(x)/(1+exp(x))
}


#'This function is used to return the basis of a spline function.
#'@name bspline
#'@param x param 1
#'@param k knots
#'@param i param 2
#'@param m param 3
#'@return basis of a spline
bspline <- function(x,k,i,m) {
  if (m==-1) {
    res <- as.numeric(x<k[i+1] & x>=k[i])
  } else {
    z0  <- (x-k[i]) / (k[i+m+1]-k[i]);
    z1  <- (k[i+m+2]-x) / (k[i+m+2]-k[i+1])
    z0*bspline(x,k,i,m-1) + z1*bspline(x,k,i+1,m-1)
  }  }

#'Dirichlet multinomial density function
#'@name dDirMult
#'@param M params of the dirichlet, i.e. model strain distribution
#'@param rho correlation parameter = 1/sample size
#'@param n category counts from the survey
#'@return M density
#'@export
dDirMult <- function(M,n,Rho) {
  if(dim(as.matrix(M))[2]==1) {
    M <- M/sum(M)
  } else {
    M <- M/rowSums(M)
  }
  M <- as.matrix(M)
  ### this condition is necessary to check for the vector inputs from
  ### likelihoods which have a variable number of bins for the distribution
  ### this is due to missing/censored data in TB reporting
  if (ncol(M)==1){
    sum(lgamma(n+M/Rho))-sum(lgamma(M/Rho))
  } else{
    rowSums(lgamma(n+M/Rho))-rowSums(lgamma(M/Rho))
  }
}
#'Simple function to return the mode of a vector
#'@name getmode
#'@param x some vector (character or numeric)
#'@return mode
#'@export
getmode <- function(v) {
  uniqv <- unique(v)
  max_match <- max(tabulate(match(v, uniqv)))
  if (max_match == 1){
    print("There is no mode. Inspect manually.")
  }
  else {
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
}

trunc_number_n_decimals <- function(numberToTrunc, nDecimals){
  numberToTrunc <- numberToTrunc + (10^-(nDecimals+5))
  splitNumber <- strsplit(x=format(numberToTrunc, digits=20, format=f), split="\\.")[[1]]
  decimalPartTrunc <- substr(x=splitNumber[2], start=1, stop=nDecimals)
  truncatedNumber <- as.numeric(paste0(splitNumber[1], ".", decimalPartTrunc))
  return(truncatedNumber)
}

### LIST OF LOCATION VECTORS IN TERMS OF ABSOLUTE INCIDENCE OTIS
ordered_locs<-c(
  "CA", "TX", "NY", "FL", "IL",
  "NJ", "GA", "PA", "MD", "VA",
  "MA", "NC", "WA", "AZ", "OH",
  "MN", "TN", "HI", "IN", "MI",
  "LA", "AL", "SC", "MO", "OR", #25
  "MS", "AR", "OK", "NV", "KY",
  "CO", "AK", "CT", "IA", "WI",
  "NM", "DC", "KS", "NE", "DE",
  "RI", "UT", "ID", "ME", "ND",
  "NH", "SD", "WV", "MT", "VT",
  "WY"
)
