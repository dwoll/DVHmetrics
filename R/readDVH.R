#####---------------------------------------------------------------------------
## returns a list (1 component per DVH file) of lists (1 component = 1 list per structure)
readDVH <- function(x, type=c("Eclipse", "Cadplan", "Masterplan", "Pinnacle"),
                    planInfo=FALSE, add) {
    type <- match.arg(type)

    dvhRawL <- if(missing(x)) {
        parseDVH(type=type)
    } else {
        parseDVH(x, type=type)
    }
    
    parseFun <- switch(type,
                       Eclipse=parseEclipse,
                       Cadplan=parseCadplan,
                       Masterplan=parseMasterplan,
                       Pinnacle=mergePinnaclePat)
    
    dvhLL <- if(length(dvhRawL) >= 1L) {
        res <- Map(parseFun, dvhRawL, planInfo=planInfo)
        Filter(Negate(is.null), res)
    } else {
        warning("No files selected")
        NULL
    }

    ## check if result should be added to existing object
    dvhLL <- if(!missing(add)) {
        if(inherits(add, "DVHLstLst")) {
            ## check if add is also organized by patient
            if(attributes(add)$byPat) {
                c(dvhLL, add)
            } else {
                addByPat <- reorgByPat(add, byPat=TRUE)
                c(dvhLL, addByPat)
            }
        } else {
            warning("add is not a DVHLstLst object - not added")
            dvhLL
        }
    } else {
        dvhLL
    }

    if(length(unique(names(dvhLL))) < length(dvhLL)) {
        warning(c("Some DVHs are for the same patient ID -",
                  "this will lead to problems in constraint checking"))
    }

    ## organized by patient (top level)
    if(!is.null(dvhLL)) {
        attr(dvhLL, which="byPat") <- TRUE
        class(dvhLL) <- "DVHLstLst"
    }

    return(dvhLL)
}

#####---------------------------------------------------------------------------
## add new DVHs to existing ones
mergeDVH <- function(...)  {
    dots <- list(...)
    stopifnot(all(vapply(dots, function(x) inherits(x, "DVHLstLst"), logical(1))))
    
    ## make sure all DVHs are organized either by patient or by structure
    orgs <- vapply(dots, function(x) attributes(x)$byPat, logical(1))
    org1 <- orgs[1]    # organization of the first DVH

    dvhLL <- if(any(orgs != org1)) {
        res <- lapply(dots, function(x) {
            if(attributes(x)$byPat == org1) {
                x
            } else {
                reorgByPat(x, org1)
            }
        })

        unlist(res, recursive=FALSE)
    } else {
        unlist(dots, recursive=FALSE)
    }

    ## organized by patient (top level)
    attr(dvhLL, which="byPat") <- org1
    class(dvhLL) <- "DVHLstLst"

    return(dvhLL)
}

#####---------------------------------------------------------------------------
## trim whitespace on beginning/end of string
trimWS <- function(x, side="both")  {
    side <- match.arg(side, c("left", "right", "both"))
    pattern <- switch(side, left="^\\s+", right="\\s+$", both="^\\s+|\\s+$")
    gsub(pattern, "", x)
}

#####---------------------------------------------------------------------------
## collapse whitespace into one space
collWS <- function(x)  {
    gsub("[[:blank:]]+", " ", x)
}

#####---------------------------------------------------------------------------
## remove whitespace everywhere
removeWS <- function(x)  {
    gsub("[[:blank:]]+", "", x)
}
