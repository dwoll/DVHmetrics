runGUI <-
function(...) {
    appDir <- system.file("DVHshiny", package="DVHmetrics")
    if(!nzchar(appDir)) {
        stop("Could not find Shiny directory. Try re-installing 'DVHmetrics'.", call.=FALSE)
    }

    if(requireNamespace("DT", quietly=TRUE)) {
        shiny::runApp(appDir, ...)
    } else {
        stop("Could not find package 'DT'. Please install first.", call.=FALSE)
    }
}
