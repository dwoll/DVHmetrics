\name{readDVH}
\alias{readDVH}
\title{Read DVH text files}
\description{
Reads single or multiple DVH text files as exported from Varian Eclipse(TM), CadPlan(TM), or MasterPlan(TM). Supports cumulative and differential DVHs.
}
\usage{
readDVH(x, type=c("Eclipse", "Cadplan", "Masterplan"), planInfo=FALSE, add)
}
\arguments{
  \item{x}{\code{character} vector giving paths to DVH text files. May contain globbing symbols understood by \code{\link[base]{Sys.glob}}. If missing and in interactive mode, \code{readDVH} opens a file selector widget. Under Windows, this widget allows selecting multiple files simultaneously.}
  \item{type}{\code{character} string. Indicates which program the DVH text files were exported from. Supported: \code{"Cadplan"} (tested with version 6.4.7), \code{"Eclipse"} (tested with versions 10-13), or \code{"Masterplan"} (tested with version 4.3).}
  \item{planInfo}{Experimental: Either \code{FALSE} or \code{character} string. In the latter case, \code{readDVH} tries to extract additional information from the \code{Plan} field in the DVH file, e.g., the prescription dose for a sum plan or the boost quadrant. Undocumented, see source.}
  \item{add}{\code{DVHLstLst} object. Existing object that should be merged with the new data from the files.}
}
\details{
Absolute dose values need to be given in Gy or cGy, absolute volume values in cm^3.

Differential DVHs are automatically converted to cumulative DVHs, but the differential DVH information is kept.

Sum plans are supported.
}
\value{
Returns an object of class \code{DVHLstLst}. This is a list (one component with class \code{DVHLst} for each original file from one patient) of lists (each component is an object of class \code{DVHs}). A \code{DVHs} object is a list with the following components:
\describe{
\item{\code{dvh}}{\code{matrix} - DVH values}
\item{\code{patID}}{\code{character} string - patient ID}
\item{\code{date}}{\code{character} string - date of DVH export}
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
# pick DVH files interactively
res <- readDVH()
res

# read all txt files in subdirectory DVH
res <- readDVH("DVH/*.txt", type="Eclipse")
res
}
}
