% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FLselex.R
\name{selex.fwd}
\alias{selex.fwd}
\title{function to forcast nyears under different selex patterns}
\usage{
selex.fwd(
  sel,
  stock,
  sr = NULL,
  fyears = 50,
  Fref = NULL,
  nyears = 3,
  plim = 0.975,
  fbar = FALSE,
  vbgf = NULL
)
}
\arguments{
\item{sel}{list of selex parameters of class of FLPars()}

\item{stock}{stock object of class FLStock}

\item{sr}{optional spawner-recruitment function FLSR}

\item{fyears}{number of forecase years}

\item{Fref}{option to input current sel refpts = c("F0","Fmsy","F0.1","Fspr30","Fsq")}

\item{nyears}{number of years for reference conditions}

\item{plim}{set fbar for ages with Selectivy >= plim (default 0.975)}

\item{fbar}{option to not correct for Fbar}
}
\value{
FLStocks object
}
\description{
function to forcast nyears under different selex patterns
}
