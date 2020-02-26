\name{TUWmodel_dual}
\alias{TUWmodel_dual}
\alias{hbvmodel_dual}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Lumped dual layer hydrological model developed at the Vienna University of Technology for education purposes}
\description{
 \code{TUWmodel_dual} is a lumped conceptual rainfall-runoff model with dual representation of soil layer developed at the TUW, 
             following the structure of the HBV model. The dual soil layer consists of a skin soil layer which represents the layer observed 
             by satellite soil moisture sensors and a root zone soil storage which is identical with the original TUWmodel concept.
             The model runs on a daily or shorter timestep and has similar structure as the original TUWmodel function.
             More information about the dual layer concept and its implementation can be found in Parajka, J., V. Naeimi, G. Bloeschl 
             and J. Komma (2009) Matching ERS scatterometer based soil moisture patterns with simulations of a conceptual dual layer 
             hydrologic model over Austria, Hydrol. Earth Syst. Sci., 13, 259-271.
}
\usage{
 TUWmodel_dual (prec, airt, ep, area=1, 
            param=c(1.2,1.2,2,-2,0,0.9,100,3.3,0.5,9,105,50,2,10,26.5,1,0.8,10), 
            incon=c(50,0,2.5,2.5,0), itsteps=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{prec}{vector/matrix of precipitation input [mm/timestep] (\code{ncol} = number of zones)}
 \item{airt}{vector/matrix of air temperatures [degC]}
 \item{ep}{vector/matrix of potential evapotranspiration [mm/timestep]}
 %\item{snowd}{vector of snow depths (same length as \code{prec}) (cm)}
 \item{area}{if more zones, vector of the percentage of area for each zone (or proportional to it, i.e., if the sum is different from 1, it will be rescaled to be 1)}
 \item{param}{vector/matrix of parameters (\code{ncol} = number of zones):
\enumerate{
\item \code{SCF} snow correction factor [-] (e.g., 0.9-1.5);
\item \code{DDF} degree day factor [mm/degC/timestep] (e.g., 0.0-5.0 mm/degC/day);
\item \code{Tr} threshold temperature above which precipitation is rain [degC] (e.g., 1.0-3.0 degC);
\item \code{Ts} threshold temperature below which precipitation is snow [degC] (e.g., -3.0-1.0 degC);
\item \code{Tm} threshold temperature above which melt starts [degC] (e.g., -2.0-2.0 degC);
\item \code{LPrat} parameter related to the limit for potential evaporation [-] (e.g., 0.0-1.0);
\item \code{FC} field capacity, i.e., max soil moisture storage [mm] (e.g., 0-600 mm);
\item \code{BETA} the non linear parameter for runoff production [-] (e.g., 0.0-20.0);
\item \code{k0} storage coefficient for very fast response [timestep] (e.g., 0.0-2.0 days);
\item \code{k1} storage coefficient for fast response [timestep] (e.g., 2.0-30.0 days);
\item \code{k2} storage coefficient for slow response [timestep] (e.g., 30.0-250.0 days);
\item \code{lsuz} threshold storage state, i.e., the very fast response start if exceeded [mm] (e.g., 1.0-100.0 mm);
\item \code{cperc} constant percolation rate [mm/timestep] (e.g., 0.0-8.0 mm/day);
\item \code{bmax} maximum base at low flows [timestep] (e.g., 0.0-30.0 days);
\item \code{croute} free scaling parameter [timestep^2/mm] (e.g., 0.0-50.0 days^2/mm);
\item \code{fc_skin} field capacity, i.e., max soil moisture storage, of the top soil skin layer [mm] (e.g., 0.1-10 mm);
\item \code{f_eta} ???? [??] (e.g., 0.7-0.95 ??);
\item \code{k_moist} ???? [??] (e.g., 5-15 ??);
}
}
 \item{incon}{vector/matrix of initial conditions for the model (\code{ncol} = number of zones):
\code{SSM0} initial value of soil moisture [mm];
\code{SWE0} initial value of snow water equivalent [mm];
\code{SUZ0} initial value for fast (upper zone) response storage [mm];
\code{SLZ0} initial value for slow (lower zone) response storage [mm];
\code{SSL0} initial value of top soil skin layer storage [mm]}
 \item{itsteps}{length of the output (if NULL all the time series are used)}
}
\details{
%No details for now.
More details about the model structure are given in 

 Parajka, J., V. Naeimi, G. Bloeschl and J. Komma (2009) Matching ERS scatterometer based soil moisture patterns with simulations of a conceptual dual layer hydrologic model over Austria, Hydrol. Earth Syst. Sci., 13, 259-271.

}
\value{
 \code{TUWmodel_dual} gives a vector of simulated runoff as \code{q} [mm/timestep], and the following vector/matrices:
\enumerate{
\item \code{qzones} simulated runoff for each zone [mm/timestep];
\item \code{q0} surface runoff [mm/timestep];
\item \code{q1} subsurface runoff [mm/timestep];
\item \code{q2} baseflow [mm/timestep];
\item \code{rain} liquid precipitation [mm/timestep];
\item \code{snow} solid precipitation [mm/timestep];
\item \code{melt} snowmelt [mm/timestep];
\item \code{moist} soil moisture [mm];
\item \code{swe} snow water equivalent [mm];
\item \code{eta} actual evapo-transpiration [mm/timestep];
\item \code{suz} upper storage zone [mm];
\item \code{slz} lower storage zone [mm];
\item \code{ssl} top soil skin layer storage [mm];
}
}
%\references{}
%\author{Alberto Viglione}
%\note{}

%\seealso{}

\examples{
## Load the data
data(example_TUWmodel)


## Simulate runoff and plot observed vs simulated series
## Lumped case (weighted means of the inputs)
simLump <- TUWmodel_dual(prec=apply(P_Vils, 1, weighted.mean, w=areas_Vils), 
                     airt=apply(T_Vils, 1, weighted.mean, w=areas_Vils),
                     ep=apply(PET_Vils, 1, weighted.mean, w=areas_Vils), 
                     area=sum(areas_Vils),
              param=c(1.02,1.70,2,0,-0.336,
                      0.934,121,2.52,
                      0.473,9.06,142,
                      50.1,2.38,10,25,1,0.8,10))

plot(as.Date(names(Q_Vils)), Q_Vils, type="l", xlab="", ylab="Discharges [mm/day]")
 lines(as.Date(rownames(T_Vils)), simLump$q, col=2)
legend("topleft", legend=c("Observations","Simulations"), col=c(1,2), lty=1, bty="n")

plot(as.Date(rownames(SWE_Vils)), apply(SWE_Vils, 1, weighted.mean, w=areas_Vils), 
     type="l", xlab="", ylab="Snow Water Equivalent [mm]")
 lines(as.Date(rownames(T_Vils)), simLump$swe, col=2)

## Distribute input case (6 zones)
simDist <- TUWmodel_dual(prec=P_Vils, airt=T_Vils, ep=PET_Vils, area=areas_Vils/sum(areas_Vils),
             param=c(1.02,1.70,2,0,-0.336,
                     0.934,121,2.52, 
                     0.473,9.06,142,
                     50.1,2.38,10,25,1,0.8,10))

plot(as.Date(names(Q_Vils)), Q_Vils, type="l", xlab="", ylab="Discharges [mm/day]")
 lines(as.Date(rownames(T_Vils)), simDist$q, col=2)
legend("topleft", legend=c("Observations","Simulations"), col=c(1,2), lty=1, bty="n")

plot(as.Date(rownames(SWE_Vils)), apply(SWE_Vils, 1, weighted.mean, w=areas_Vils),
     type="l", xlab="", ylab="Snow Water Equivalent [mm]")
 lines(as.Date(rownames(T_Vils)), apply(simDist$swe, 1, weighted.mean, w=areas_Vils), col=2)

## Distributed input and parameters case
parametri <- matrix(rep(c(1.02,1.70,2,0,-0.336,
                          0.934,121,2.52,
                          0.473,9.06,142,
                          50.1,2.38,10,25,1,0.8,10), 6), ncol=6)
parametri[2,] <- c(1.4, 1.7, 1.9, 2.2, 2.4, 3.0)
simDist2 <- TUWmodel_dual(prec=P_Vils,
                      airt=T_Vils, 
                      ep=PET_Vils, 
                      area=areas_Vils/sum(areas_Vils),
                param=parametri)

plot(as.Date(names(Q_Vils)), Q_Vils, type="l", xlab="", ylab="Discharges [mm/day]")
 lines(as.Date(rownames(T_Vils)), simDist2$q, col=2)
legend("topleft", legend=c("Observations","Simulations"), col=c(1,2), lty=1, bty="n")
}
\keyword{univar}
