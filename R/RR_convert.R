################################################################################
########### THIS CODE IS TO CONVERT THE USER INPUTTED RISK RATIOS    ###########
########### INTO THE LEVELS NECESSARY FOR THE MODEL PARAMETERS FOR
########### THE GENERIC RISK GROUPS.
########### INPUTS: SINGLE RISK RATIO BY USER
########### OUTPUTS: RATE OF SLOW PROGRESSION OF TB                  ###########
###########          PROBABILITY OF FAST PROGRESSION OF TB           ###########
###########          PROBABILITY OF FAST PROGRESSION OF TB WITH PI   ###########
###########
################################################################################

################################################################################
########### SET THE MAXIMUM ACCEPTABLE RISK RATIO THAT THE MODEL     ###########
########### WILL ACCEPT AND RUN PROPERLY.
max_risk_factor <- 40
min_risk_factor <- 1
################################################################################
###########   ASK THE USER TO INPUT A RISK RATIO FOR THE MODEL RUN   ###########
start<-function(risk_factor=0) {
  assign_RR <- function(risk_factor){
    if (risk_factor < max_risk_factor)
      risk_factor <- risk_factor
    else risk_factor <- readline("Enter an integer for the increased risk of TB reactivation risk for your risk group of interest:")
  }
  assign_RR (risk_factor)
    while (risk_factor < min_risk_factor | risk_factor > max_risk_factor ){
      if (risk_factor !=0)
        cat("The ratio must be above", min_risk_factor, "and below", max_risk_factor, "! Please enter a new integer value:")
      assign_RR(risk_factor)
    }
}
########### AS IT IS CURRENTLY CODED THE MODEL CANNOT HANDLE A       ###########
########### NEGATIVE OR PROTECTIVE RR. ALSO THE ERROR MESSAGE
########### IS ALSO CURRENTLY BROKEN FOR
risk_factor <- risk_factor

################################################################################
################################################################################
########### CONVERT THE USER INPUTTED RATIO INTO USABLE RATIOS FOR   ###########
########### THE MODEL; THESE ARE BASED ON LINEAR RELATIONSHIPS       ###########
########### WHICH ARE EXTRAPOLATED FROM THE RELATIONSHIP OF HIV      ###########
########### PARAMETERS FROM THE NATIONAL MODEL. VALUES ARE CALCULATED###########
########### FOR THE HIGHEST LEVEL OF REACTIVATION RISK IN THIS WAY.  ###########
################################################################################
################################################################################

################################################################################
###########   ODDS RATIO FOR PROBABILITY OF FAST PROGRESSION FOR RF  ###########
################################################################################
ORpfastRF   <-risk_factor
################################################################################
###########    ODDS RATIO FOR PROB OF FAST PROGRESSION W/PI FOR RF   ###########
################################################################################
ORpfastPIRF <-risk_factor
################################################################################
###########      RATE RATIO FOR RATE OF SLOW PROGRESSION FOR RF      ###########
################################################################################
RRslowRF    <- risk_factor

################################################################################
###########  CREATE A VECTOR FOR ALL FOUR LEVELS OF THE RISK FACTOR  ###########
################################################################################
###########  DEPENDENT ON THE USER INPUTTED RATIO THE MODEL WILL     ###########
###########  DETERMINE THE NUMBER OF TB REACTIVATION LEVELS TO USE   ###########

cut4 <- 20
cut3 <- 10
cut2 <- 5


  if (risk_factor > cut4){
    m <- log(risk_factor)/3
} else if (risk_factor < cut4 & risk_factor > cut3){
    m <- log(risk_factor)/2
} else if (risk_factor < cut3 & risk_factor > cut2){
    m <- log(risk_factor)
} else {
    m <-1
}

vORpfastRF   <-c(1, )

vORpfastPIRF <-c(1,)

vRRslowRF    <-c(1, )


