% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_simulations.R
\name{OutputsInt}
\alias{OutputsInt}
\title{wrapper function for the above function}
\usage{
OutputsInt(
  loc,
  ParMatrix,
  n_cores = 1,
  output_month = 11,
  startyr = 1950,
  endyr = 2050,
  Int1 = 0,
  Int2 = 0,
  Int3 = 0,
  Int4 = 0,
  Int5 = 0,
  Scen1 = 0,
  Scen2 = 0,
  Scen3 = 0,
  prg_chng,
  ttt_list,
  par2020 = c(0.4232265, 0.3707595, 0.1984619, 1.1158255),
  care_cascade = def_care_cascade()
)
}
\arguments{
\item{loc}{two-digit abbreviation for location}

\item{ParMatrix}{parameters to use in the simulation}

\item{n_cores}{how many cores to use}

\item{output_month}{when to output results}

\item{startyr}{year to start the simulation}

\item{endyr}{year to end the simulation}

\item{Int1}{boolean for intervention 1}

\item{Int2}{boolean for intervention 2}

\item{Int3}{boolean for intervention 3}

\item{Int4}{boolean for intervention 4}

\item{Int5}{boolean for intervention 5}

\item{Scen1}{boolean for scenario 1}

\item{Scen2}{boolean for scenario 2}

\item{Scen3}{boolean for scenario 3}

\item{prg_chng}{vector of program change values}

\item{ttt_list}{list of targeted testing and treatment values}

\item{par2020}{vector of 2020 adjustment parameters}
}
\value{
out outputs
}
\description{
wrapper function for the above function
}
