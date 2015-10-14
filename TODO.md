# TODO
## Medium term

 * showDVH(..., addMSD=TRUE) -> linearly interpolate dose on a grid first?
   MasterPlan does not repeat identical volume values with different doses
 * document showDVH(..., addMSD=TRUE) in vignettes
 * allow `checkConstraint(x, "DNTCP < 50%")` instead of `0.5Gy`
 * `getMetric()`, `checkConstraint()`, `readDVH()` -> multicore support
 * `readDVH()` -> if `x` is a directory, read all files in it
 * `readDVH()` -> if `x` is a zip file, read all files in it
 * read files exported from iPlan, BrainLab

## Long term

 * expand Jupyter notebook
 * Shiny app: better way to save old data instead of temp-file
 * `getDMEAN()` `smooth.spline()` -> `tol` must be strictly positive
 * `getDMEAN()` bandwidth choice for `locpoly()` when `dpill()` fails
 * `getDMEAN()` adaptive binning grid size for `locpoly()` with small bandwidth
 * merge DVH files from the same patient ID
 * add Eclipse fields: Comment, Exported by, Description, Plan Status, Approval Status, Dose Coverage, Sampling Coverage, Equivalent Sphere Diameter, Conformity Index, Gradient Measure
 * convert `RadOnc` objects to `DVHmetrics` objects and back
 * `showConstraint()` -> draw constraint arrows as custom geoms (coming in `ggplot2` 1.1.0)
 * harmonize structures by reading equivalence file, possibly with regex

# build
 * "c:\program files\r\r-3.2.1\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.3.3.tar.gz
 * "c:\program files\r\r-3.2.1\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-3.2.1\bin\x64\Rcmd.exe" check DVHmetrics_0.3.3.tar.gz --as-cran
 * `install.packages("d:/daniel_work/workspace/DVHmetrics_0.3.3.tar.gz", repos=NULL, type="source")`

 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.3.3.tar.gz
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" check DVHmetrics_0.3.3.tar.gz --as-cran
 * `install.packages("h:/workspace/DVHmetrics_0.3.3.tar.gz", repos=NULL, type="source")`
