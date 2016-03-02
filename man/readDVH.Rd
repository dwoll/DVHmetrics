\name{readDVH}
\alias{readDVH}
\title{Read DVH text files}
\description{
Reads single or multiple DVH text files as exported from Varian Eclipse(TM), CadPlan(TM), OnCentra MasterPlan(TM), Philipps Pinnacle3 (TM), Elekta Monaco (TM), Tomo HiArt (TM), or RaySearch Labs RayStation (TM). Supports cumulative and differential DVHs.
}
\usage{
readDVH(x,
        type=c("Eclipse", "Cadplan", "Masterplan",
               "Pinnacle", "Monaco", "HiArt", "RayStation"),
        planInfo=FALSE, courseAsID=FALSE, add)
}
\arguments{
  \item{x}{\code{character} vector giving paths to DVH text files. May contain globbing symbols understood by \code{\link[base]{Sys.glob}}. If missing and in interactive mode, \code{readDVH} opens a file selector widget. Under Windows, this widget allows selecting multiple files simultaneously. For \code{type="Pinnacle"}, \code{x} should be one of the following: A directory with information for one patient, a directory with several sub-directories (one for each patient), or a zip file of such directories. Under Windows, if \code{x} is missing and \code{type="Pinnacle"}, \code{readDVH} opens a folder selector widget.}
  \item{type}{\code{character}. Indicates which program the DVH text files were exported from. Supported: \code{"Cadplan"} (tested with version 6.4.7), \code{"Eclipse"} (tested with Varian Eclipse version 8), \code{"Masterplan"} (tested with OnCentra MasterPlan version 4.3), \code{"Pinnacle"} (tested with Pinnacle3 version 9, see Details), \code{"Monaco"} (tested with Elekta Monaco version 5), \code{"HiArt"} (TomoTherapy HiArt), or \code{"RayStation"} (RaySearch Labs RayStation).}
  \item{planInfo}{Experimental: Either \code{FALSE} or \code{character} string. In the latter case, \code{readDVH} tries to extract additional information from the \code{Plan} field in the DVH file, e.g., the prescription dose for a sum plan or the boost quadrant. Undocumented, see source.}
  \item{courseAsID}{\code{logical}. If \code{TRUE}, the \code{Course} entry in the header section of a DVH file is appended to the regular patient ID. Currently supported only for \code{type="Eclipse"}.}
  \item{add}{\code{DVHLstLst} object. Existing object that should be merged with the new data from the files.}
}
\details{
Absolute dose values need to be given in Gy or cGy, absolute volume values in cm^3.

Differential DVHs are automatically converted to cumulative DVHs, but the differential DVH information is kept.

Sum plans are supported.

For RayStation, only cumulative DVHs with absolute volume are currently supported. Volume is assumed to be measured in cm^3.

To export data from Tomo HiArt, copy to clipboard and then save to file from a text editor. Support for Tomo HiArt files is currently limited to those with absolute dose and relative volume. Please send an anonymized sample file if you need this functionality. Since HiArt files do not contain info on patient ID, the current workaround is to generate a random ID.

Pinnacle3 files have to be exported using its own scripting facility such that information from one patient is contained in one directory. A suitable export script is available on request from the package authors. The directory layout for one patient has to be as follows (experimental, likely to change in future versions):
\itemize{
\item{Files (CSV format with column headers):
\itemize{
\item{\code{DoseInfo.csv} (variables \code{"PrescriptionDose cGy"}, \code{"NumberOfFractions"}, \code{"Dosis cGy"})}
\item{\code{PatInfo.csv} (variables \code{"LastName"}, \code{"FirstName"}, \code{"MedicalRecordNumber"})}
\item{\code{PlanInfo.csv} (variable \code{"PlanName"})}}}
\item{Directory: \code{Data}:
\itemize{
\item{\code{Info.csv} (variables \code{"Filename"}, \code{"RegionOfInterestName"},  \code{"DoseMin cGy"}, \code{"DoseMax cGy"}, \code{"DoseMean cGy"}, \code{"Volume ccm"})}
\item{DVH1.csv, DVH2.csv, ...  - the actual DVH data files with names defined in \code{Info.csv} variable \code{"Filename"}. They should look like
\verb{NumberOfDimensions = 2;
NumberOfPoints = 431;
Points[] =\{
  0,0
  10,0
  ...
  4000,100
\};
  }
}}}}
}
\value{
Returns an object of class \code{DVHLstLst}. This is a list (one component with class \code{DVHLst} for each original file from one patient) of lists (each component is an object of class \code{DVHs}). A \code{DVHs} object is a list with the following components:
\describe{
\item{\code{dvh}}{\code{matrix} - cumulative DVH values}
\item{\code{dvhDiff}}{\code{matrix} - differential DVH values, only created a) if original file contained a differential DVH or b) by \code{\link[DVHmetrics]{convertDVH}}}
\item{\code{patID}}{\code{character} string - patient ID}
\item{\code{date}}{\code{character} string - date of DVH export}
\item{\code{type}}{\code{character} string - cumulative or differential DVH}
\item{\code{plan}}{\code{character} string - plan name}
\item{\code{course}}{\code{character} string - course - currently Eclipse only}
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
\code{\link[DVHmetrics]{checkConstraint}},
\code{\link[DVHmetrics]{convertDVH}}
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
