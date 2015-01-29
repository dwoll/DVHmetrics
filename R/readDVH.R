#####---------------------------------------------------------------------------
## returns a list (1 component per DVH file) of lists (1 component = 1 list per structure)
readDVH <- function(x, type=c("Eclipse", "Cadplan", "Masterplan"), planInfo=FALSE) {
    type <- match.arg(type)

    dvhRawL <- if(missing(x)) {
        parseDVH()
    } else {
        parseDVH(x)
    }

    parseFun <- switch(type,
                       Eclipse=parseEclipse,
                       Cadplan=parseCadplan,
                       Masterplan=parseMasterplan)

    dvhLL <- lapply(dvhRawL, parseFun, planInfo=planInfo)

    if(length(unique(names(dvhLL))) < length(dvhLL)) {
        warning(c("Some DVHs are for the same patient ID -",
                  "this will lead to problems in constraint checking"))
    }

    ## organized by patient (top level)
    attr(dvhLL, which="byPat") <- TRUE
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
