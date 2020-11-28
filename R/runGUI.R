runGUI <-
function(...) {
    appDir <- system.file("DVHshiny", package="DVHmetrics")
    if(!nzchar(appDir)) {
        stop("Could not find Shiny directory. Try re-installing 'DVHmetrics'.", call.=FALSE)
    }

    if(requireNamespace("DT", quietly=TRUE)) {
        ## check if we have bs4Dash for newer GUI
        if(requireNamespace("bs4Dash", quietly=TRUE)) {
            ## breaking changes introduced in bs4Dash 2.0.0
            ## check which version is available
            bs4Dash_version <- packageVersion("bs4Dash")
            if(compareVersion("2.0.0", as.character(bs4Dash_version)) == -1) {
                shiny::runApp(appDir, ...)
            } else {
                appDir_bs4Dash_old <- paste0(appDir, "_bs4Dash_05")
                shiny::runApp(appDir, ...)
            }
        } else {
            warning("Package 'bs4Dash' not found - running legacy version")
            appDir_legacy <- paste0(appDir, "_legacy")
            shiny::runApp(appDir_legacy, ...)
        }
    } else {
        stop("Could not find package 'DT'. Please install first.", call.=FALSE)
    }
}
