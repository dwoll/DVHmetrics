run_gui <-
function(...) {
    appDir <- system.file("MeshAgreement", package="MeshAgreement")
    if(!nzchar(appDir)) {
        stop("Could not find Shiny directory. Try re-installing 'MeshAgreement'.", call.=FALSE)
    }

    have_shiny   <- requireNamespace("shiny",   quietly=TRUE)
    have_bs4Dash <- requireNamespace("bs4Dash", quietly=TRUE)
    have_DT      <- requireNamespace("DT",      quietly=TRUE)
    have_rgl     <- requireNamespace("rgl",     quietly=TRUE)
    
    if(!have_shiny)   { warning("Package 'shiny' not found.") }
    if(!have_bs4Dash) { warning("Package 'bs4Dash' not found.") }
    if(!have_DT)      { warning("Package 'DT' not found.") }
    if(!have_rgl)     { warning("Package 'rgl' not found.") }
    
    if(have_shiny && have_bs4Dash && have_DT && have_rgl) {
        shiny::runApp(appDir, ...)
    } else {
        wants  <- c("shiny", "bs4Dash", "DT", "rgl")
        have   <- c(have_shiny, have_bs4Dash, have_DT, have_rgl)
        ip_str <- sprintf("install.packages('%s')", paste0(wants[!have], collapse="', '"))
        stop(sprintf("Run %s to install missing packages.", ip_str))
    }
}
