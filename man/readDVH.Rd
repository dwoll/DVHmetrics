\name{readDVH}
\alias{readDVH}
\title{Read DVH text files}
\description{
Reads single or multiple DVH text files as exported from Cadplan or Varian Eclipse. Supports cumulative and differential DVHs.
}
\usage{
readDVH(x, type=c("Eclipse", "Cadplan"), planInfo=FALSE)
}
\arguments{
  \item{x}{\code{character} vector giving paths to DVH text files. May contain globbing symbols understood by \code{\link[base]{Sys.glob}}. If missing and in interactive mode, \code{readDVH} opens a file selector widget. Under Windows, this widget allows selecting multiple files simultaneously.}
  \item{type}{\code{character} string. Indicates which program the DVH text files were exported from. Supported: \code{"Cadplan"} or \code{"Eclipse"}.}
  \item{planInfo}{Experimental: Either \code{FALSE} or \code{character} string. In the latter case, \code{readDVH} tries to extract additional information from the \code{Plan} field in the DVH file, e.g., the prescription dose for a sum plan or the boost quadrant. Undocumented, see source.}
}
\details{
Absolute dose values need to be given in Gy or cGy, absolute volume values in cm\^3.

Differential DVHs are automatically converted to cumulative DVHs, but the differential DVH information is kept.

Sum plans are supported.

For Eclipse, dose and volume can be both absolute and relative.
For Cadplan, dose can be absolute and relative, volume must be absolute.

Tested with text files exported from Cadplan version 8 and Varian Eclipse version 10 and 11.
}
\value{
Returns an object of class \code{DVHLstLst}. This is a list (one component with class \code{DVHLst} for each original file from one patient) of lists (each component is an object of class \code{DVHs}). A \code{DVHs} object is a list with the following components:
\describe{
\item{\code{dvh}}{\code{matrix} - DVH values}
\item{\code{patID}}{\code{character} string - patient ID}
\item{\code{type}}{\code{character} string - cumulative or differential DVH}
\item{\code{plan}}{\code{character} string - plan name}
\item{\code{structure}}{\code{character} string - structure name}
\item{\code{structVol}}{\code{numeric} - structure volume}
\item{\code{doseUnit}}{\code{character} string - measurement unit dose}
\item{\code{volumeUnit}}{\code{character} string - measurement unit volume}
\item{\code{doseRx}}{\code{numeric} - prescription dose}
\item{\code{isoDoseRx}}{\code{numeric} - iso-dose percentage}
\item{\code{doseMin}}{\code{numeric} - minimum dose from DVH file}
\item{\code{doseMax}}{\code{numeric} - maximum dose from DVH file}
\item{\code{doseAvg}}{\code{numeric} - average dose from DVH file}
\item{\code{doseMed}}{\code{numeric} - median dose from DVH file}
\item{\code{doseSD}}{\code{numeric} - dose standard deviation from DVH file}
}
}
\seealso{
\code{\link[base]{Sys.glob}},
\code{\link[DVHmetrics]{print.DVHs}},
\code{\link[DVHmetrics]{showDVH}},
\code{\link[DVHmetrics]{getMetric}},
\code{\link[DVHmetrics]{checkConstraint}}
}
\examples{
\dontrun{
# not run
# pick DVH files interactively
res <- readDVH()
res

# read all txt files in subdirectory DVH
res <- readDVH("DVH/*.txt", type="Eclipse")
res
}
}
