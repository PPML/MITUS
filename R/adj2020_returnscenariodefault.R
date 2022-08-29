def_returnScenario <- function(){

Immig <- list(
  "return_months" = 865:888,
  "multiplier" = 1
)

rDxt <- list(
  "return_months" = 865:888,
  "multiplier" = 1
)

Trans <- list(
  "return_months" = 865:888,
  "multiplier" = 1
)

CaseFat <- list(
  "return_months" = 865:888,
  "multiplier" = 1
)

ReturnParams <- list("Immig" = Immig,
                     "rDxt" = rDxt,
                     "Trans" = Trans,
                     "CaseFat" = CaseFat)

return(ReturnParams)
}
