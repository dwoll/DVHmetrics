\name{getNTCP}
\alias{getNTCP}
\alias{getNTCP.DVHs}
\alias{getNTCP.DVHLst}
\alias{getNTCP.DVHLstLst}
\title{Normal tissue complication probability (NTCP)}
\description{
Calculate normal tissue complication probability (NTCP) from Lyman's probit model, Niemierko's logit model, the Poisson model, or the Kaellman relative seriality model. May be based on EQD2.
}
\usage{
getNTCP(x,
        NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL, NTCPs=NULL,
        EUDa=NULL, EUDfn=NULL, EUDab=NULL,
        NTCPtype=c("probit", "logit", "poisson", "relative_seriality"), ...)

\method{getNTCP}{DVHs}(x,
        NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL, NTCPs=NULL,
        EUDa=NULL, EUDfn=NULL, EUDab=NULL,
        NTCPtype=c("probit", "logit", "poisson", "relative_seriality"), ...)

\method{getNTCP}{DVHLst}(x,
        NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL, NTCPs=NULL,
        EUDa=NULL, EUDfn=NULL, EUDab=NULL,
        NTCPtype=c("probit", "logit", "poisson", "relative_seriality"), ...)

\method{getNTCP}{DVHLstLst}(x,
        NTCPtd50=NULL, NTCPm=NULL, NTCPn=NULL, NTCPgamma50=NULL, NTCPs=NULL,
        EUDa=NULL, EUDfn=NULL, EUDab=NULL,
        NTCPtype=c("probit", "logit", "poisson", "relative_seriality"), ...)
}
\arguments{
  \item{x}{One cumulative DVH (object of class \code{DVHs}, multiple cumulative DVHs from one patient with multiple structures (object of class \code{DVHLst}), or multiple cumulative DVHs from many patients, each with multiple structures (object of class \code{DVHLstLst}). See \code{\link[DVHmetrics]{readDVH}}.}
  \item{NTCPtd50}{Tolerance dose with 50\% complication probability.}
  \item{NTCPm}{Probit/logit Parameter m. Equal to 1 / (NTCPgamma50*sqrt(2*pi)).}
  \item{NTCPn}{Parameter n. Equal to 1/a with exponential gEUD parameter a.}
  \item{NTCPgamma50}{Poisson parameter gamma50. Equal to 1 / (NTCPm*sqrt(2*pi))}
  \item{NTCPs}{Relative seriality parameter s.}
  \item{EUDa}{If gEUD should be based on EQD2: Exponential parameter a.}
  \item{EUDfn}{If gEUD should be based on EQD2: Number of fractions.}
  \item{EUDab}{If gEUD should be based on EQD2: alpha/beta ratio for the relevant tissue.}
  \item{NTCPtype}{\code{"probit"} - Lyman probit model, \code{"logit"} - Niemierko logit model, \code{"poisson"} - Poisson model, \code{"relative_seriality"} - Kaellmann relative seriality model.}
  \item{...}{Ignored. Used to catch additional arguments passed from \code{\link[DVHmetrics]{getMetric}}.}
}
\details{
For the logit, probit, and Poisson method, gEUD is used for DVH reduction. This is equivalent to the Kutcher-Burman DVH reduction scheme. The probit model is given in equation (1), the logit model in eqation (2), and the Poisson model in equation (3) in Kaellman (1992), with gEUD plugged in for D. The relative seriality model is given in equation (18).
}
\value{
A data frame with variables \code{NTCP}, \code{patID}, and \code{structure}.
}
\references{
Kaellman, P., Agren, A., & Brahme, A. (1992). Tumor and normal tissue responses to fractionated non-uniform dose delivery. International Journal of Radiation Biology, 62(2), 249-262.

Kutcher, G. J., Burman, C., Brewster L., Goitein, M., & Mohan, R. (1991). Histogram reduction method for calculating complication probabilities for threedimensional treatment planning evaluations. International Journal of Radiation Oncology Biology Physics, 21(1), 137-146.

Lyman, J. T. (1985). Complication probability as assessed from dose volume histograms. Radiation Research, 104(2), S13-19.

Niemierko, A. (1999). A generalized concept of equivalent uniform dose. Medical Physics, 26(6), 1100.

Rancati et al. (2004). Fitting late rectal bleeding data using different NTCP models: results from an Italian multi-centric study (AIROPROS0101). Radiotherapy Oncology, 73, 21-32.
}
\seealso{
\code{\link[DVHmetrics]{getTCP}},
\code{\link[DVHmetrics]{getEUD}},
\code{\link[DVHmetrics]{getMetric}}
}
\examples{
## treatment was in 2 Gy fractions
getNTCP(dataMZ[[1]][["HEART"]],
        NTCPtd50=48, NTCPm=0.6, NTCPn=0.5, NTCPtype="probit")

getNTCP(dataMZ[[1]][["HEART"]],
        NTCPtd50=52.3, NTCPgamma=1.28, NTCPs=1, NTCPtype="relative_seriality")
}
