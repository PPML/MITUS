#' Basic Functions Used In Package
#' This script includes several functions that are useful for data manipulation
#' throughout the package.

#'The first function takes the odds of a risk, then calculates the join
#'probability of that odds and and an inputted Odds Ratio that the user inputs, then
#'generates the corresponding risk.
#'@param val known risk for population of interest
#'@param OR odds ratio by which the known risk will be modified
#'@return new risk for population


ORAdd <- function(val,OR) {
  x <- (val/(1-val))*OR; x/(1+x)
}

#'This function calculates the logit of the value or parameter inputted.
#'@title lgt
#'@name logit
#'@param x value of interest, often an odds
#'@return logit of input x

lgt <-  function(x) log(x/(1-x));

#'This function calculates the inverse logit of the value or parameter inputted.
#'@param x value of interest, often an odds
#'@return inverse logit of input x

invlgt <- function(x) 1/(1+exp(-x))

#'This function is used to simulate a logistic increase between two years of interest
#'to a user-specified maximum. Useful for simulating slow increase/decrease of a
#'parameter over time.
#'@param StYr Year to start the logit curve
#'@param Endyr Year to end the logit curve
#'@param EndVal maximum value that the logistic curve can reach
#'@return vector of values demonstrating a logistic increase from zero to EndVal

LgtCurve <- function(StYr,Endyr,EndVal) {
  z <- log(1/0.005-1)
  zz  <- seq(-z*(1+2*(StYr-1950)/(Endyr-StYr)),
              z*(1+2*(2050-Endyr)/(Endyr-StYr)),
             by=(2*z)/(Endyr-StYr)/12)
  zz  <- as.numeric(EndVal)/(1+exp(-zz))
  if(StYr>1950) {
    zz[1:((StYr-1950)*12)] <- 0
  }
  zz  }

#'This function creates a smooth curve from a vector of values.
#'@param vec vector of values that the user would like to have smoothed into a curve
#'@return expanded vector of values that determine the shape of the smoothing spline

SmoCurve <- function(vec) {
  jj <- predict(smooth.spline(x=1:length(vec),y=vec,spar=0.2),
                x=seq(1,length(vec),1/12))$y;
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
#'@param x
#'@param k
#'@param i
#'@param m
#'@return basis of a spline

bspline <- function(x,k,i,m) {
  if (m==-1) {
    res <- as.numeric(x<k[i+1] & x>=k[i])
} else {
    z0  <- (x-k[i]) / (k[i+m+1]-k[i]);
    z1  <- (k[i+m+2]-x) / (k[i+m+2]-k[i+1])
    z0*bspline(x,k,i,m-1) + z1*bspline(x,k,i+1,m-1)
}  }
