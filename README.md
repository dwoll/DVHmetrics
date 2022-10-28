# MeshAgreement

Agreement measures for 3D structures saved as mesh files. R package containing an interactive [R](https://www.r-project.org/) [Shiny](https://shiny.rstudio.com/) app. You can upload mesh files (STL, PLY, OBJ, OFF) to generate agreement measures for all pairwise comparisons, as well as the corresponding aggregated agreement. The intended application is to compare delineated structures for radiotherapy treatment planning.

**Currently unstable under heavy development.**

# Implemented agreement measures

 * Pairwise distance-based and volume-overlap-based metrics
     * DCOM: Distance between centers of mass
     * ASD: Average surface distance
     * RMSD: Root mean squared surface distance
     * HD_max: Hausdorff distance - max of both directed HDs
     * HD_avg: Hausdorff distance - average of both directed HDs
     * HD_95:  Hausdorff distance - average of both 95th percentiles of directed surface distances
 * JSC: Jaccard similarity coefficient
 * DSC: Dice similarity coefficient

# Required packages

`MeshAgreement` heavily relies on packages developed by Stéphane Laurent, currently only available from GitHub:

  * [PolygonSoup](https://github.com/stla/PolygonSoup)
  * [SurfaceReconstruction](https://github.com/stla/SurfaceReconstruction)
  * [MeshesTools](https://github.com/stla/MeshesTools)
  * [Boov](https://github.com/stla/Boov)

Required packages on CRAN:

  * [shiny](https://CRAN.R-project.org/package=shiny)
  * [bs4Dash](https://CRAN.R-project.org/package=bs4Dash)
  * [shinyWidgets](https://CRAN.R-project.org/package=shinyWidgets)
  * [dplyr](https://CRAN.R-project.org/package=dplyr)
  * [Rvcg](https://CRAN.R-project.org/package=Rvcg)
  * [RcppCGAL](https://CRAN.R-project.org/package=RcppCGAL)

# Literature (selection)

 * Fotina et al. Critical discussion of evaluation parameters for inter-observer variability in target definition for radiation therapy. Strahlenther Onkol 2012; 188: 160-167.
 * Hanna  et al. Geometrical Analysis of Radiotherapy Target Volume Delineation: a Systematic Review of Reported Comparison Methods. Clin Oncol 2010; 22, 515-525.
 * Heimann et al. Comparison and Evaluation of Methods for Liver Segmentation From CT Datasets. IEEE Trans Med Imaging 2009; 28: 1251-1265.
 * Sherer et al. Metrics to evaluate the performance of auto-segmentation for radiation treatment planning: A critical review. Radiother Oncol 2021; 160: 185-191."
