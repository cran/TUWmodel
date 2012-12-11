\name{example_TUWmodel}
\docType{data}
\alias{example_TUWmodel}
\alias{Q_Vils}
\alias{P_Vils}
\alias{T_Vils}
\alias{PET_Vils}
\alias{SWE_Vils}
\alias{areas_Vils}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Data-sample}
\description{
Data for examples in the package
}
\usage{
 data(example_TUWmodel)
}
\format{The data are time series for the river Vils at Vils (Laende):
\enumerate{
\item \code{Q_Vils} vector of observed daily discharges (mm/day);
\item \code{P_Vils} matrix of observed daily precipitation (mm/day) for 6 zones;
\item \code{T_Vils} matrix of observed daily temperatures (degC) for 6 zones;
\item \code{PET_Vils} matrix of potential evapotranspiration (mm/day) for 6 zones;
\item \code{SWE_Vils} matrix of observed snow water equivalent (mm) for 6 zones;
\item \code{areas_Vils} areas of the 6 zones.
}
}
%\source{}
%\references{}
%\author{Alberto Viglione}
%\note{}

\examples{
data(example_TUWmodel)

ls()

t <- as.Date(rownames(P_Vils))
plot(t, apply(P_Vils, 1, weighted.mean, w=areas_Vils), 
     type="S", xlab="", ylab="Precipitation [mm/day]")

plot(t, apply(PET_Vils, 1, weighted.mean, w=areas_Vils), 
     type="l", xlab="", ylab="Potential evapotranspiration [mm/day]")

plot(t, T_Vils[,1], type="l", xlab="", ylab="Temperature [degC]")
 lines(t, T_Vils[,6], col=3)

plot(t, SWE_Vils[,6], col=3, type="l", xlab="", ylab="Snow water equivalent [mm]")
 lines(t, SWE_Vils[,1], col=1)

plot(as.Date(names(Q_Vils)), Q_Vils, type="l", xlab="", ylab="Runoff [mm/day]")
}
\keyword{datasets}
