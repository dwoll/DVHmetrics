runGUI <-
function(...) {
    shiny::runApp(appDir=system.file("DVHshiny", package="DVHmetrics"), ...)    
}
