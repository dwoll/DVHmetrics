\name{MeshAgreement-package}
\alias{MeshAgreement-package}
\alias{MeshAgreement}
\docType{package}
\title{Agreement measures for 3D meshes}
\description{
Calculates distance-based and volume-overlap-based
agreement measures for triangular 3D meshes. These include
the Hausdorff distance, the average surface distance, the
Dice similarity coefficient, and the Jaccard similarity
coefficient as documented in Stockinger et al. (2021)
<doi:10.1186/s13014-021-01965-5>. Overall agreement for
a set of meshes is calculated as the aggregate agreement
for all pairwise comparisons. Based on algorithms
provided by the 'VCGLIB' <http://vcg.isti.cnr.it/vcglib/>
and 'CGAL' <https://www.cgal.org/> libraries.
Includes web-based graphical user interface.
}
\details{
\tabular{ll}{
Package: \tab MeshAgreement\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2022-11-14\cr
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

\url{http://vcg.isti.cnr.it/vcglib/}
}
\keyword{package}
\examples{
heartL <- mesh3dL_to_cgalMeshL(data_heart_obsL)
pairL  <- get_mesh_pairs(heartL)
get_mesh_agree_pair(pairL[[1]], do_ui=FALSE)
}
