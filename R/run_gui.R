run_gui <-
function(...) {
    appDir <- system.file("MeshAgreement", package="MeshAgreement")
    if(!nzchar(appDir)) {
        stop("Could not find Shiny directory. Try re-installing 'MeshAgreement'.", call.=FALSE)
    }

    shiny::runApp(appDir, ...)
}
