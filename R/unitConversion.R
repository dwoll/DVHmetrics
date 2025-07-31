## conversion of dose and volume units
getConvFac <-
function(conversion="CGY2GY") {
    conversion <- toupper(removeWS(conversion))

    ## check how the conversion factor is indicated
    idxCGY2GY  <- conversion == "CGY2GY"
    idxCGY2EVG <- conversion == "CGY2EV/G"

    idxGY2CGY  <- conversion == "GY2CGY"
    idxGY2EVG  <- conversion == "GY2EV/G"

    idxEVG2CGY <- conversion == "EV/G2CGY"
    idxEVG2GY  <- conversion == "EV/G2GY"

    idxCGY2CGY <- conversion == "CGY2CGY"
    idxGY2GY   <- conversion == "GY2GY"
    idxEVG2EVG <- conversion == "EV/G2EV/G"

    idxCC2CC   <- conversion == "CC2CC"

    ## did we catch all requested conversion types?
    idxAll <- idxCGY2GY  | idxCGY2EVG |
              idxGY2CGY  | idxGY2EVG  |
              idxEVG2CGY | idxEVG2GY  |
              idxCGY2CGY | idxGY2GY   | idxEVG2EVG | idxCC2CC

    if(!all(idxAll)) {
        warning(c('Conversion type(s) "',
                  toString(conversion[!idxAll]),
                  '" not found - conversion factor set to NA'))
    }

    convFac <- rep(NA_real_, length(conversion))

    ## conversion factors for dose units
    convFac[idxCGY2GY]  <- 1/100
    convFac[idxCGY2EVG] <- NA_real_

    convFac[idxGY2CGY]  <- 100
    convFac[idxGY2EVG]  <- NA_real_

    convFac[idxEVG2CGY] <- NA_real_
    convFac[idxEVG2GY]  <- NA_real_

    convFac[idxCGY2CGY] <- 1
    convFac[idxGY2GY]   <- 1
    convFac[idxEVG2EVG] <- 1

    convFac[idxCC2CC]   <- 1

    return(convFac)
}

## determine unit from conversion string
getUnits <-
function(x="CGY2GY", first=TRUE) {
    if(!is.character(x)) {
        warning("Unit not recognized - input must have form like CGY2GY")
        return(NA_character_)
    }

    units    <- strsplit(x, "2", fixed=TRUE)  # first and second part of string
    unitLens <- lengths(units)                # count parts
    if(!all(unitLens == 2L)) {                # check that there are two parts
        warning("Unit not recognized - input must have form like CGY2GY")
        return(NA_character_)
    }

    knownUnits <- c("CGY", "GY", "EV/G", "CC")
    isKnown    <- vapply(units, function(x) { all(x %in% knownUnits) }, logical(1))
    if(!all(isKnown)) {
        warning(c("Unit not recognized - needs to be one of\n",
                  paste(knownUnits, collapse=" ")))
        return(NA_character_)
    }

    if(first) {
        sapply(units, head, n=1)
    } else {
        sapply(units, tail, n=1)
    }
}
