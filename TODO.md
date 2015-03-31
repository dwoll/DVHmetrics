# TODO
## Short term

 * Constraint "DNTCP < 0.5Gy" -> "DNTCP < 0.5P"
 * Vignette -> BED, EQD2, isoeffective dose, NTCP, TCP
 * getMetric(), checkConstraint(), readDVH() -> multicore support
 * getDVHmean() smooth.spline() -> tol must be strictly positive
 * getDVHmean() bandwidth choice for locpoly() when dpil() fails
 * getDVHmean() adaptive binning grid size for locpoly() with small bandwidth

## Long term

 * shiny app: better way to save old data instead of temp-file
 * merge DVH files from the same patient ID
 * add Eclipse fields: Comment, Exported by, Description, Course, Plan Status, Approval Status, Dose Coverage, Sampling Coverage, Equivalent Sphere Diameter, Conformity Index, Gradient Measure
 * `readDVH()` -> if `x` is a directory, read all files in it
 * make `getMetrics()`, `showDVH()`, `checkConstraints()`, `showConstraints()` work with `RadOnc` objects
 * read files exported from Philips Pinnacle3, TomoTherapy Hi-Art, iPlan, BrainLab
 * draw constraint arrows as custom geoms
 * harmonize structures by reading equivalence file, possibly with regex

# build
 * "c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.1.tar.gz
 * "c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" check DVHmetrics_0.1.0.900.tar.gz --as-cran
 * `install.packages("d:/daniel_work/workspace/DVHmetrics_0.1.0.900.tar.gz", repos=NULL, type="source")`

 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.1.tar.gz
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" check DVHmetrics_0.1.0.900.tar.gz --as-cran
 * `install.packages("h:/workspace/DVHmetrics_0.1.tar.gz", repos=NULL, type="source")`
