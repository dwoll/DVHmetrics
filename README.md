# DVHmetrics

## Analyze Dose-Volume Histograms and Check Constraints using R

`DVHmetrics` provides functionality for the analysis of dose-volume histograms in radiation oncology:

 * Read dose-volume-histogram (DVH) text files
 * Calculate arbitrary DVH metrics, show DVH diagrams
 * Calculate gEUD, BED, convert DVH doses to EQD2
 * Calculate TCP and NTCP according to LKB (Lyman probit), Niemierko logit, or Kaellman Poisson (relative seriality) model
 * Check and visualize quality assurance constraints for the DVH
 * Includes a [`shiny`](http://shiny.rstudio.com/)-based web application for most of the functionality - also available online at [shiny.imbei.uni-mainz.de:3838/DVHmetrics/](http://shiny.imbei.uni-mainz.de:3838/DVHmetrics/).
 * For further explanations and an example walkthrough, see the [package vignette](http://cran.rstudio.com/web/packages/DVHmetrics/vignettes/DVHmetrics.pdf)
