\name{get_name_elem}
\alias{get_name_elem}
\title{Return element of a comparison string}
\description{Utility function to return the first or second element of a comparison string used to indicate which two meshes are being compared.}
\usage{
get_name_elem(x, which=1, sep=" <-> ")
}
\arguments{
  \item{x}{\code{character}. The string that indicates the comparison.}
  \item{which}{\code{numeric}. Indicates whether to return the first or second element.}
  \item{sep}{\code{character}. The string that indicates the separator between the first and second element.}
}
\details{Thin front-end to \code{\link[base]{gsub}}. Only useful internally or for Shiny GUI.}
\value{\code{character} element as the first or second element of the string.}
\seealso{
\code{\link[base]{gsub}},
\code{\link[MeshAgreement]{get_mesh_pairs}}
}
\examples{
get_name_elem("MESH_1 <-> MESH_2", which=1)
}
