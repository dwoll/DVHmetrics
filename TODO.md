# TODO

 * `readDVH(..., type="HiArt")` -> remove all non-numbers as these may creep in during copying from the clipboard, add option to externally supply absolute volumes / prescription dose
 * `readDVH()`: use `tryCatch()` to diagnose problems with input files containing umlauts / accents that are not in UTF-8 encoding

## Medium term

 * `getNTCP()` relative seriality model (equation 18 in Kaellman et al.), from differential DVH with density / per unit dose
 * check that expensive calculations (e.g., conversion to differential DVH) are done only once, and then stored
 * conformity index CI?
 * shiny app: shinydashboard?
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

# build
 * "c:\program files\r\r-3.4.2patched\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.3.7.tar.gz
 * "c:\program files\r\r-3.4.2patched\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-3.4.2patched\bin\x64\Rcmd.exe" check DVHmetrics_0.3.7.tar.gz --as-cran
 * `install.packages("d:/daniel_work/workspace/DVHmetrics_0.3.7.tar.gz", repos=NULL, type="source")`

 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.3.7.tar.gz
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" check DVHmetrics_0.3.7.tar.gz --as-cran
 * `install.packages("h:/workspace/DVHmetrics_0.3.7.tar.gz", repos=NULL, type="source")`
