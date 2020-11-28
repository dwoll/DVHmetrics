### About DVHmetrics

This web application is a front-end for the [R](http://www.r-project.org/) package [`DVHmetrics`](https://github.com/dwoll/DVHmetrics) that provides functionality for radiation oncology: Use it to read dose-volume-histogram (DVH) text files, to calculate free DVH metrics, to show DVH diagrams, to check and visualize quality assurance constraints for the DVH.

#### Authors

[`DVHmetrics`](http://CRAN.R-project.org/package=DVHmetrics) and this web application are written by:

Daniel Wollschläger [<wollschlaeger@uni-mainz.de>](mailto:wollschlaeger@uni-mainz.de)  
Institute for Medical Statistics, Epidemiology und Informatics ([IMBEI](http://www.imbei.de/))  
University Medical Center Mainz  
Germany

Heiko Karle [<karle@uni-mainz.de>](mailto:karle@uni-mainz.de")  
Department of Radiation Oncology  
University Medical Center Mainz  
Germany

Data courtesy of Department of Radiation Oncology (Prof. Dr. Schmidberger), University Medical Center Mainz, Germany.

#### Documentation

For documentation on how to use the package and this web application, see:  
[DVHmetrics vignette](http://cran.rstudio.com/web/packages/DVHmetrics/vignettes/DVHmetrics.pdf)  
[Web application vignette](http://cran.rstudio.com/web/packages/DVHmetrics/vignettes/DVHshiny.pdf)

Source code at: [github.com/dwoll/DVHshiny/](https://github.com/dwoll/DVHmetrics/tree/master/inst/DVHshiny)

#### Acknowledgements

Thanks to Marcus Stockinger and Michael R. Young for testing and feedback, to Sandra Bührdel, Hannes Rennau, Ulrich Wolf, Bjorne Riis, Nico Banz, and Michael R. Young for sample DVH files.
Created with the [Shiny](http://shiny.rstudio.com/) web application framework powered by [RStudio](http://www.rstudio.com/).

Uses functionality provided by the R packages [ggplot2](http://cran.r-project.org/package=ggplot2), [reshape2](http://cran.r-project.org/package=reshape2), [grid](http://cran.r-project.org/package=grid), and [KernSmooth](http://cran.r-project.org/package=KernSmooth).

#### References

Chang W, Cheng J, Allaire JJ, Xie Y, McPherson J (2020). shiny: Web Application Framework for R.  
R package version 1.5.0.    
[http://CRAN.R-project.org/package=shiny](http://CRAN.R-project.org/package=shiny)

David Granjon (2019). bs4Dash: A 'Bootstrap 4' Version of 'shinydashboard'. R package version 0.5.0. [https://CRAN.R-project.org/package=bs4Dash](https://CRAN.R-project.org/package=bs4Dash)

R Core Team (2020). R: A language and environment for statistical computing.  
R Foundation for Statistical Computing, Vienna, Austria.  
[http://www.R-project.org/](http://www.R-project.org/)

Wand M (2019). KernSmooth: Functions for Kernel Smoothing Supporting Wand & Jones (1995).  
R package version 2.23-16.  
[http://CRAN.R-project.org/package=KernSmooth](http://CRAN.R-project.org/package=KernSmooth)

Wickham H (2007). Reshaping Data with the reshape Package.  
Journal of Statistical Software, 21(12), 1-20.  
[http://www.jstatsoft.org/v21/i12/](http://www.jstatsoft.org/v21/i12/)  
[http://cran.r-project.org/package=reshape2](http://cran.r-project.org/package=reshape2)

Wickham H (2016). ggplot2: elegant graphics for data analysis. 2nd ed.  
New York: Springer.  
[http://cran.r-project.org/package=ggplot2](http://cran.r-project.org/package=ggplot2)
