melt_dw <- function(x, valname="value", level=1) {
    if(!is.list(x)) {
        setNames(data.frame(x), valname)
    } else {
        xL      <- lapply(x, melt_dw, level=level + 1, valname=valname)
        lengths <- vapply(xL, nrow, integer(1))
        names   <- if(!is.null(names(x))) {
            names(x)
        } else { seq_along(x) }

        x_comb <- do.call("rbind", xL)
        x_comb[[paste0("L", level)]] <- rep(names, lengths)
        rownames(x_comb) <- NULL
        x_comb
    }
}

reshape_dfL <- function(y, idvar) {
    varyvar <- setdiff(names(y), idvar)
    yL <- reshape(y,
                  direction="long",
                  idvar=idvar,
                  varying=varyvar,
                  v.names="value",
                  timevar="variable")

    yL[["variable"]] <- varyvar[yL[["variable"]]]
    rownames(yL) <- NULL
    yL
}
