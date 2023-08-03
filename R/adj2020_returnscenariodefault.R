def_returnScenario <- function(){

# default is from DHS OIS LPR Annual flow found here:
# https://www.dhs.gov/sites/default/files/2023-02/2022_0405_plcy_lawful_permanent_residents_fy2021v2.pdf
Immig <- list(
  "return_months" = 845:862,
  "multiplier" = 1.33
)
# default is informed by the household pulse survey found here:
# https://www.cdc.gov/nchs/covid19/rands/reduced-access-to-care.htm
# Mean of the prior is set to 0.40 reduction.
# Extrapolating the rate of decline in the survey, we estimate 18 months
# to return to normal.
rDxt <- list(
  "return_months" = 844:861,
  "multiplier" = 1
)

# default value based on data from COVID-19 US State Policy (CUSP) database
# data found here: https://statepolicies.com
# 12 months of strong impact,
# 12 months return to normal
Trans <- list(
  "return_months" = 855:866,
  "multiplier" = 1
)

# default value based on CDC Wonder Multiple Cause of Death data
# found here: https://wonder.cdc.gov/mcd-icd10-expanded.html
CaseFat <- list(
  "return_months" = 844:890,
  "multiplier" = 1
)

ReturnParams <- list("Immig" = Immig,
                     "rDxt" = rDxt,
                     "Trans" = Trans,
                     "CaseFat" = CaseFat)

return(ReturnParams)
}
