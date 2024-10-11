# TODO

 * Vignette: The radiation effect (E) or cell kill can be defined as follows: E = n*(alpha*d + beta*d^2), divide both sides by alpha
 * getMetric(dvh_m, metric="D5cc") works even though volumeUnit is NA
 * invalid option "error" -> use paste() first inside stop()
 * x <- readDVH("c:/Users/Daniel/Downloads/DVH poumon volume absolute.txt", type="HiArt", volume_from_dvh=TRUE, hiart=list(doseRx=50))
 Fehler in mapply(FUN = f, ..., SIMPLIFY = FALSE) : # list(list(doseRx=50)) ?
  Eingaben mit Länge 0 können nicht mit Eingaben anderer Länge gemischt werden
 * `readDVH()`: use `tryCatch()` to diagnose problems with input files containing umlauts / accents that are not in UTF-8 encoding / use con <- file(..., encoding="UTF8"), then readLines(con, encoding="")

## Medium term

 * check that expensive calculations (e.g., conversion to differential DVH) are done only once, and then stored
 * conformity index CI?
 * allow `checkConstraint(x, "DNTCP < 50%")` instead of `0.5Gy`
 * `getMetric()`, `checkConstraint()`, `readDVH()` -> multicore support
 * `readDVH()` -> if `x` is a directory, read all files in it
 * `readDVH()` -> if `x` is a zip file, read all files in it
 * read files exported from iPlan, BrainLab

## Long term

 * shiny app: better way to save old data instead of temp-file
 * `getDMEAN()` `smooth.spline()` -> `tol` must be strictly positive
 * `getDMEAN()` bandwidth choice for `locpoly()` when `dpill()` fails
 * `getDMEAN()` adaptive binning grid size for `locpoly()` with small bandwidth
 * merge DVH files from the same patient ID
 * add Eclipse fields: Comment, Exported by, Description, Plan Status, Approval Status, Dose Coverage, Sampling Coverage, Equivalent Sphere Diameter, Conformity Index, Gradient Measure
 * convert `RadOnc` objects to `DVHmetrics` objects and back
 * `showConstraint()` -> draw constraint arrows as custom geoms (`ggplot2` 2.1.0 https://github.com/hadley/ggplot2/blob/master/vignettes/extending-ggplot2.Rmd)
 * harmonize structures by reading equivalence file, possibly with regex

# Build

Sys.setenv(R_GSCMD="C:\\Program Files\\gs\\gs9.26\\bin\\gswin32c.exe")

 * "c:\program files\r\r-4.2.1\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.3.8.tar.gz
 * "c:\program files\r\r-4.2.1\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-4.2.1\bin\x64\Rcmd.exe" check DVHmetrics_0.3.9.tar.gz --as-cran
 * `install.packages("d:/daniel_work/workspace/DVHmetrics_0.3.9.tar.gz", repos=NULL, type="source")`

 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.3.9.tar.gz
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" check DVHmetrics_0.3.9.tar.gz --as-cran
 * `install.packages("h:/workspace/DVHmetrics_0.3.9.tar.gz", repos=NULL, type="source")`
