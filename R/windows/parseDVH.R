parseDVH <- function(x) {
    files <- if(!missing(x)) {
        Sys.glob(x)
    } else if(interactive()) {  # && (.Platform$OS.type == "windows"))
        ## are we are in interactive mode AND under Windows?
    	## we are under Windows since this sits in a platform-specific directory
        ## choose files interactively
        myFilt <- rbind(Filters, txtCsvDat=c("Data files (*.txt, *.csv, *.dat)",
                                             "*.txt;*.csv;*.dat"))
        files <- choose.files(filters=myFilt[c("txtCsvDat", "All"), ], index=1)
    } else {
        character(0)
    }

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
