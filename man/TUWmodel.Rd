\name{TUWmodel}
\alias{TUWmodel}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Lumped hydrological model developed at the Vienna University of Technology for education purposes}
\description{
 \code{TUWmodel} is a lumped conceptual rainfall-runoff model, following the structure of the HBV model. 
             The model runs on a daily time step and consists of a snow routine, a soil moisture routine and a flow routing routine. 
             See Parajka, J., R. Merz, G. Bloeschl (2007) Uncertainty and multiple objective calibration in regional water balance modelling: case study in 320 Austrian catchments, Hydrological Processes, 21, 435-446.
}
\usage{
 TUWmodel (prec, airt, ep, area, param=c(1.2,1.2,2,-2,0,0.9,100,3.3,0.5,9,105,50,2,10,26.5), incon=c(50,0,2.5,2.5), itsteps=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{prec}{vector/matrix of precipitation input (mm/day) (\code{ncol} = number of zones)}
 \item{airt}{vector/matrix of air temperatures (degC)}
 \item{ep}{vector/matrix of potential evapotranspiration (mm/day)}
 %\item{snowd}{vector of snow depths (same length as \code{prec}) (cm)}
 \item{area}{catchment area (km2) (vector if more zones)}
 \item{param}{vector of parameters:
\enumerate{
\item \code{SCF} snow correction factor (0.9-1.5);
\item \code{DDF} degree day factor (0.0-5.0 mm/degC/day);
\item \code{Tr} threshold temperature above which precipitation is rain (1.0-3.0 degC);
\item \code{Ts} threshold temperature below which precipitation is snow (-3.0-1.0 degC);
\item \code{Tm} threshold temperature above which melt starts (-2.0-2.0 degC);
\item \code{LPrat} parameter related to the limit for potential evaporation (0.0-1.0);
\item \code{FC} field capacity, i.e., max soil moisture storage (0-600 mm);
\item \code{BETA} the non linear parameter for runoff production (0.0-20.0);
\item \code{k0} storage coefficient for very fast response (0.0-2.0 days);
\item \code{k1} storage coefficient for fast response (2.0-30.0 days);
\item \code{k2} storage coefficient for slow response (30.0-250.0 days);
\item \code{lsuz} threshold storage state, i.e., the very fast response start if exceeded (1.0-100.0 mm);
\item \code{cperc} constant percolation rate (0.0-8.0 mm/day);
\item \code{bmax} maximum base at low flows (0.0-30.0 day);
\item \code{croute} free scaling parameter (0.0-50.0 day2/mm);
}
see Parajka, J., R. Merz, G. Bloeschl (2007) Uncertainty and multiple objective calibration in regional water balance modelling: case study in 320 Austrian catchments, Hydrological Processes, 21, 435-446, doi:10.1002/hyp.6253.}
 \item{incon}{initial conditions for the model:
\code{SSM0} soil moisture (mm);
\code{SWE0} snow water equivalent (mm);
\code{SUZ0} initial value for fast (upper zone) response storage (mm);
\code{SLZ0} initial value for slow (lower zone) response storage (mm)}
 \item{itsteps}{length of the output (if NULL all the time series are used)}
}
\details{
No details for now.
}
\value{
 \code{TUWmodel} gives a vector of simulated runoff as \code{q} (m3/s), and the following vector/matrices:
\enumerate{
\item \code{qzones} simulated runoff for each zone (m3/s);
\item \code{swe} snow water equivalent (mm);
\item \code{q0} surface runoff (m3/s);
\item \code{q1} subsurface runoff (m3/s);
\item \code{q2} baseflow (m3/s);
\item \code{rmoist} relative soil moisture (between 0 and 1);
\item \code{rain} liquid precipitation (mm/day)
\item \code{snow} solid precipitation (mm/day)
\item \code{eta} actual evapotranspiration (mm/day)
}
}
%\references{}
%\author{Alberto Viglione}
%\note{}

%\seealso{}

\examples{
## ------------------------------------------------------------------------------------------------------ ##
## Load the data
data(example_TUWmodel)


## Simulate runoff and plot observed vs simulated series
## Lumped case (weighted means of the inputs)
simLump <- TUWmodel(prec=apply(P_Vils, 1, weighted.mean, w=areas_Vils), 
                    airt=apply(T_Vils, 1, weighted.mean, w=areas_Vils),
                    ep=apply(PET_Vils, 1, weighted.mean, w=areas_Vils), 
                    area=sum(areas_Vils),
             param=c(1.02,1.70,2,0,-0.336,0.934,121,2.52,0.473,9.06,142,50.1,2.38,10,25))

plot(as.Date(names(Q_Vils)), Q_Vils, type="l", xlab="", ylab="Discharges [mm/day]")
 lines(as.Date(rownames(T_Vils)), simLump$q, col=2)
legend("topleft", legend=c("Observations","Simulations"), col=c(1,2), lty=1, bty="n")

plot(as.Date(rownames(SWE_Vils)), apply(SWE_Vils, 1, weighted.mean, w=areas_Vils), 
     type="l", xlab="", ylab="Snow Water Equivalent [mm]")
 lines(as.Date(rownames(T_Vils)), simLump$swe, col=2)

## Distribute case (6 zones)
simDist <- TUWmodel(prec=P_Vils, airt=T_Vils, ep=PET_Vils, area=areas_Vils,
             param=c(1.02,1.70,2,0,-0.336,0.934,121,2.52,0.473,9.06,142,50.1,2.38,10,25))

plot(as.Date(names(Q_Vils)), Q_Vils, type="l", xlab="", ylab="Discharges [mm/day]")
 lines(as.Date(rownames(T_Vils)), simDist$q, col=2)
legend("topleft", legend=c("Observations","Simulations"), col=c(1,2), lty=1, bty="n")

plot(as.Date(rownames(SWE_Vils)), apply(SWE_Vils, 1, weighted.mean, w=areas_Vils),
     type="l", xlab="", ylab="Snow Water Equivalent [mm]")
 lines(as.Date(rownames(T_Vils)), apply(simDist$swe, 1, weighted.mean, w=areas_Vils), col=2)
}
\keyword{univar}
