# TODO
 * DVHshiny -> switch for `getMetric(..., interp)` option
 * merge DVH files from the same patient ID
 * add Eclipse fields: Comment, Exported by, Description, Course, Plan Status, Approval Status, Dose Coverage, Sampling Coverage, Equivalent Sphere Diameter, Conformity Index, Gradient Measure
 * `readDVH()` -> if `x` is a directory, read all files in it
 * differential -> cumulative: doses are interval mid-points?
 * make `getMetrics()`, `showDVH()`, `checkConstraints()`, `showConstraints()` work with `RadOnc` objects
 * gEUD, BED, NTCP, TCP
 * read files exported from Philips Pinnacle3, Helax, TomoTherapy Hi-Art, iPlan
 * draw constraint arrows as custom geoms
 * harmonize structures by reading equivalence file, possibly with regex

# build
 * "c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.1.tar.gz
 * "c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-3.1.0\bin\x64\Rcmd.exe" check DVHmetrics_0.1.tar.gz --as-cran
 * `install.packages("d:/daniel_work/workspace/DVHmetrics_0.1.0.900.tar.gz", repos=NULL, type="source")`

 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" INSTALL --build DVHmetrics_0.1.tar.gz
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" build DVHmetrics --resave-data --compact-vignettes="both"
 * "c:\program files\r\r-devel\bin\x64\Rcmd.exe" check DVHmetrics_0.1.tar.gz --as-cran
 * `install.packages("h:/workspace/DVHmetrics_0.1.tar.gz", repos=NULL, type="source")`
