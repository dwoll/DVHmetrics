\name{MeshAgreement-package}
\alias{MeshAgreement-package}
\alias{MeshAgreement}
\docType{package}
\title{Agreement measures for 3D meshes}
\description{
Calculates distance-based and volume-overlap-based
agreement measures for 3D meshes. Overall agreement for
a set of meshes is calculated as the aggregate agreement
for all pairwise comparisons.
Includes web-based graphical user interface.
}
\details{
\tabular{ll}{
Package: \tab MeshAgreement\cr
Type: \tab Package\cr
Version: \tab 0.1.0.999\cr
Date: \tab 2022-10-28\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Daniel Wollschlaeger and Heiko Karle

Maintainer: Daniel Wollschlaeger <wollschlaeger@uni-mainz.de>
}
\references{
Fotina et al. Critical discussion of evaluation parameters for inter-observer variability in target definition for radiation therapy. Strahlenther Onkol 2012; 188: 160-167.

Hanna  et al. Geometrical Analysis of Radiotherapy Target Volume Delineation: a Systematic Review of Reported Comparison Methods. Clin Oncol 2010; 22, 515-525.

Heimann et al. Comparison and Evaluation of Methods for Liver Segmentation From CT Datasets. IEEE Trans Med Imaging 2009; 28: 1251-1265.

Sherer et al. Metrics to evaluate the performance of auto-segmentation for radiation treatment planning: A critical review. Radiother Oncol 2021; 160: 185-191.

\url{https://www.cgal.org/}
}
\keyword{package}
\examples{
pairL <- get_mesh_pairs(data_heart_meshL)
get_mesh_agree_pair(pairL[[1]])
}
