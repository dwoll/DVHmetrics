\name{mergeDVH}
\alias{mergeDVH}
\title{Merge existing DVH objects}
\description{
Combine several existing DVH objects into one object.
}
\usage{
mergeDVH(...)
}
\arguments{
  \item{...}{\code{DVHLstLst} objects.}
}
\details{
The first object determines whether the resulting object is organized by patient or by structure. Objects need not originally come from the same treatment planning system.
}
\value{
Returns an object of class \code{DVHLstLst}.
}
\examples{
\dontrun{
# pick some DVH files interactively
a <- readDVH(type="Cadplan")

# pick other DVH files interactively
b <- readDVH(type="Eclipse")

# combine DVH data
res <- mergeDVH(a, b)
res
}
}
