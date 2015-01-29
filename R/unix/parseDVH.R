## under Unix-like OS, we don't have choose.files() and Filters
## TODO: fix for fPath=NA -> no initial /
parseDVH <- function(x) {
    files <- Sys.glob(x)
    if(length(files) >= 1L) {
        ## read in files into a list of character vectors
        DVHraw <- lapply(files, readLines)
    
        ## name them using patient IDs
        getPatID <- function(txt) {
            IDline <- txt[grep("^(Patient ID|Case)[[:blank:]]*:", txt)]
            IDres  <- sub("^.+?:[[:blank:]]+([[:alnum:][:punct:][:blank:]]+$)", "\\1", IDline, perl=TRUE)
            collWS(trimWS(IDres, side="both"))
        }
    
        ## patient id's as names for list components
        names(DVHraw) <- lapply(DVHraw, getPatID)
        DVHraw
    } else {
        warning("No files were selected")
        NULL
    }
}
