# TODO
## Medium term

 * `getEQD2()`, `getIsoEffD()` methods for `DVHs`, `DVHLst`, `DVHLstLst`
 * DVHshiny vignette -> document BED/EQD2 pane
 * allow `checkConstraint(x, "DNTCP < 50%")` instead of `0.5Gy`
 * `getMetric()`, `checkConstraint()`, `readDVH()` -> multicore support
 * `readConstraint()` needs to be platform dependent (`file.choose()`)
 * read files exported from TomoTherapy Hi-Art, iPlan, BrainLab
 * `readDVH()` -> if `x` is a directory, read all files in it
 * `readDVH()` -> if `x` is a zip file, read all files in it

## Long term

 * Jupyter notebook
 * Shiny app: better way to save old data instead of temp-file
 * `getDMEAN()` `smooth.spline()` -> `tol` must be strictly positive
 * `getDMEAN()` bandwidth choice for `locpoly()` when `dpill()` fails
 * `getDMEAN()` adaptive binning grid size for `locpoly()` with small bandwidth
 * merge DVH files from the same patient ID
 * add Eclipse fields: Comment, Exported by, Description, Course, Plan Status, Approval Status, Dose Coverage, Sampling Coverage, Equivalent Sphere Diameter, Conformity Index, Gradient Measure
 * convert `RadOnc` objects to `DVHmetrics` objects and back
 * `showConstraint()` -> draw constraint arrows as custom geoms
 * harmonize structures by reading equivalence file, possibly with regex

# build
 * "c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.2.tar.gz
 * "c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" check DVHmetrics_0.2.tar.gz --as-cran
 * `install.packages("d:/daniel_work/workspace/DVHmetrics_0.2.tar.gz", repos=NULL, type="source")`

 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.2.tar.gz
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" check DVHmetrics_0.2.tar.gz --as-cran
 * `install.packages("h:/workspace/DVHmetrics_0.2.tar.gz", repos=NULL, type="source")`
