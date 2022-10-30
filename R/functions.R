get_name_pair <- function(x, y, sep=" <-> ") {
    paste0(x, sep, y) 
}

get_name_elem <- function(x, which=1L, sep=" <-> ") {
    stopifnot(which %in% c(1L, 2L))
    stopifnot(grepl(sep, x))
    if(which == 1) {
        gsub(paste0("^(.+)", sep, ".+$"), "\\1", x)
    } else {
        gsub(paste0("^.+", sep, "(.+$)"), "\\1", x)
    }
}

read_mesh_one <- function(x, name,
                          reconstruct=c("none", "AFS", "Poisson"),
                          spacing) {
    reconstruct <- match.arg(reconstruct)
    
    mesh_name <- if(missing(name)) {
        basename(file_path_sans_ext(x))
    } else {
        basename(file_path_sans_ext(name))
    }
    
    mesh_0 <- readMeshFile(x)
    mesh_r <- if(reconstruct == "none") {
        toRGL(Mesh(vertices=mesh_0$vertices,
                   faces=mesh_0$faces))
    } else if(reconstruct == "AFS") {
        AFSreconstruction(mesh_0$vertices)
    } else if(reconstruct == "Poisson") {
        PoissonReconstruction(mesh_0$vertices, spacing=spacing)
    }
    
    vol_0 <- try(meshVolume(mesh_r))
    ctr_0 <- try(meshCentroid(mesh_r))
    vol <- if(!inherits(vol_0, "try-error")) {
        vol_0
    } else {
        NA_real_
    }
    
    ctr <- if(!inherits(ctr_0, "try-error")) {
        ctr_0
    } else {
        rep(NA_real_, 3)
    }
    
    list(name    =mesh_name,
         mesh    =mesh_r,
         volume  =vol,
         centroid=ctr)
}

read_mesh <- function(x, name,
                      reconstruct=c("none", "AFS", "Poisson"),
                      spacing) {
    reconstruct <- match.arg(reconstruct)
    
    mesh_names <- if(missing(name)) {
        basename(file_path_sans_ext(x))
    } else {
        basename(file_path_sans_ext(name))
    }
    
    Map(read_mesh_one,
        setNames(x, mesh_names),
        mesh_names,
        reconstruct=reconstruct,
        spacing=spacing)
}

print_mesh_one <- function(x) {
    vol_fmt_str <- if(is.na(x$volume))        { "%s" }           else { "%.2f" }
    ctr_fmt_str <- if(any(is.na(x$centroid))) { "[%s, %s, %s]" } else { "[%.2f, %.2f, %.2f]" }
    cat(sprintf("Mesh: %s", x$name),
        "\n",
        capture.output(print(x$mesh)),
        "\n",
        sprintf(paste0("Volume: ", vol_fmt_str), x$volume),
        "\n",
        sprintf(paste0("Centroid: ", ctr_fmt_str),
                x$centroid[1],
                x$centroid[2],
                x$centroid[3]),
        "\n",
        "\n", sep="")
}

print_mesh <- function(x) {
    invisible(Map(print_mesh_one, x))
}

## all pairs from list of meshes
get_mesh_pairs <- function(x, sep=" <-> ") {
    stopifnot(length(x) >= 1L)
    l <- if(length(x) >= 2L) {
        pairs_idx <- t(combn(seq_along(x), 2))
        lapply(seq_len(nrow(pairs_idx)), function(i) {
            idx1  <- pairs_idx[i, 1]
            idx2  <- pairs_idx[i, 2]
            mesh1 <- x[[idx1]]
            mesh2 <- x[[idx2]]
            list(name=get_name_pair(mesh1$name, mesh2$name, sep=sep),
                 mesh1=mesh1,
                 mesh2=mesh2)
        })
    } else {
        mesh1 <- x[[1]]
        mesh2 <- x[[1]]
        list(list(name=get_name_pair(mesh1$name, mesh2$name, sep=sep),
                  mesh1=mesh1,
                  mesh2=mesh2))
    }
    
    pair_names <- lapply(l, function(x) { x$name })
    setNames(l, pair_names)
}

## union and intersection for list of two meshes x
get_mesh_ui_pair <- function(x) {
    union_0     <- try(MeshesUnion(       list(x$mesh1$mesh, x$mesh2$mesh), clean=TRUE))
    intersect_0 <- try(MeshesIntersection(list(x$mesh1$mesh, x$mesh2$mesh), clean=TRUE))
    
    if(!inherits(union_0, "try-error") && !inherits(intersect_0, "try-error")) {
        union     <- toRGL(union_0)
        intersect <- toRGL(intersect_0)
    } else {
        union     <- NULL
        intersect <- NULL
    }
    
    list(name=get_name_pair(x$mesh1$name, x$mesh2$name),
         union=union,
         intersection=intersect)    
}

get_mesh_ui <- function(x) {
    pairL <- get_mesh_pairs(x)
    Map(get_mesh_ui_pair, pairL)
}

get_mesh_metro_pair <- function(x, chop=TRUE, ...) {
    metro <- vcgMetro(x$mesh1$mesh,
                      x$mesh2$mesh,
                      ...)
    
    if(chop) {
        metro[["forward_hist"]]  <- NULL
        metro[["backward_hist"]] <- NULL
    }
    
    metro[["name"]] <- get_name_pair(x$mesh1$name, x$mesh2$name)
    metro
}

get_mesh_metro <- function(x, chop=TRUE, ...) {
    pairL <- get_mesh_pairs(x)
    Map(get_mesh_metro_pair, pairL, chop=chop, ...)
}

## distance measures, union, intersection for each mesh pair
get_mesh_agree_pair <- function(x, metro, ui, chop=TRUE, ...) {
    ## distance-based measures
    if(missing(metro)) {
        metro <- get_mesh_metro_pair(x, chop=chop, ...)
    }
    
    DCOM   <- sqrt(sum((x$mesh2$centroid - x$mesh1$centroid)^2))
    HD_max <- max(c(metro$ForwardSampling$maxdist, metro$BackwardSampling$maxdist))
    HD_avg <- (metro$ForwardSampling$maxdist + metro$BackwardSampling$maxdist) / 2
    if((sum(metro$distances1^2) > (.Machine$double.eps^0.5)) &&
                 (sum(metro$distances2^2) > (.Machine$double.eps^0.5))) {
        HD_95  <- (quantile(metro$distances1, probs=0.95) + quantile(metro$distances2, probs=0.95)) / 2
        RMSD   <- sqrt(mean(c(metro$distances1^2, metro$distances2^2)))
        ASD    <-      mean(c(metro$distances1,   metro$distances2))
    } else {
        HD_95  <- NA_real_
        RMSD   <- NA_real_
        ASD    <- NA_real_
    }
    # n1     <- length(metro$distances1) # metro$ForwardSampling$nsamples
    # n2     <- length(metro$distances2) # metro$BackwardSampling$nsamples
    # w1     <- n1 / (n1+n2)
    # w2     <- n2 / (n1+n2)
    # ASD    <- w1*metro$ForwardSampling$meandist + w2*metro$BackwardSampling$meandist
    # RMSD   <- sqrt(w1*metro$ForwardSampling$RMSdist^2 + w2*metro$BackwardSampling$RMSdist^2)
    
    ## volume-overlap-based measures
    ## check if union/intersection are supplied
    if(missing(ui)) {
        ui <- get_mesh_ui_pair(x)
    }
    
    if(!is.null(ui$union) && !is.null(ui$intersection)) {
        vol_mesh1     <- x$mesh1$volume
        vol_mesh2     <- x$mesh2$volume
        vol_union     <- meshVolume(ui$union)
        vol_intersect <- meshVolume(ui$intersection)
        
        JSC <- vol_intersect   / vol_union
        DSC <- 2*vol_intersect / (vol_mesh1 + vol_mesh2)
    } else {
        JSC <- NA_real_
        DSC <- NA_real_
    }
    
    data.frame(mesh1=x$mesh1$name,
               mesh2=x$mesh2$name,
               DCOM=DCOM,
               HD_max=HD_max,
               HD_avg=HD_avg,
               HD_95=unname(HD_95),
               ASD=ASD,
               RMSD=RMSD,
               JSC=JSC,
               DSC=DSC)
}

get_mesh_agree <- function(x, chop=TRUE, ...) {
    pairL  <- get_mesh_pairs(x)
    uiL    <- Map(get_mesh_ui_pair,    pairL)
    metroL <- Map(get_mesh_metro_pair, pairL, chop=chop, ...)
    agreeL <- Map(get_mesh_agree_pair,
                  pairL,
                  metroL,
                  chop=chop)

    d <- do.call("rbind", agreeL)
    rownames(d) <- NULL
    d
}

get_mesh_agree_long <- function(x) {
    vars_varying <- c("DCOM",
                      "HD_max", "HD_avg", "HD_95", "ASD", "RMSD",
                      "JSC", "DSC")
    
    vars_id <- names(x)[!(names(x) %in% vars_varying)]
    
    dL <- reshape(x,
                  direction="long",
                  idvar=vars_id,
                  varying=vars_varying,
                  v.names="observed",
                  timevar="metric")
    
    rownames(dL) <- NULL
    
    dL[["metric"]] <- factor(dL[["metric"]],
                             levels=seq_along(vars_varying),
                             labels=vars_varying)
    
    dL
}

get_mesh_agree_aggr <- function(x, na.rm=FALSE) {

    d_agreeL <- get_mesh_agree_long(x)
    d_agreeL[["observed_ln"]] <- log(d_agreeL[["observed"]])
    
    d_mean   <- aggregate(observed    ~ metric, FUN=mean,   data=d_agreeL, na.rm=na.rm)
    d_median <- aggregate(observed    ~ metric, FUN=median, data=d_agreeL, na.rm=na.rm)
    d_sd     <- aggregate(observed    ~ metric, FUN=sd,     data=d_agreeL, na.rm=na.rm)
    d_var    <- aggregate(observed    ~ metric, FUN=var,    data=d_agreeL, na.rm=na.rm)
    d_varlog <- aggregate(observed_ln ~ metric, FUN=var,    data=d_agreeL, na.rm=na.rm)

    d_aggr <- Reduce(function(x, y) { suppressWarnings(merge(x, y, by="metric")) },
                     list(d_mean, d_median, d_sd, d_var, d_varlog))

    names(d_aggr)       <- c("metric", "Mean", "Median", "SD", "VAR", "VAR_log")
    d_aggr[["CV"]]      <- d_aggr[["SD"]] / d_aggr[["Mean"]]
    d_aggr[["CV_ln"]]   <- sqrt(exp(d_aggr[["VAR_log"]]) - 1)
    d_aggr[["VAR"]]     <- NULL
    d_aggr[["VAR_log"]] <- NULL
    
    d_aggr
}

get_mesh_agree_aggr_long <- function(x) {
    vars_varying <- c("Mean", "Median", "SD", "CV", "CV_ln")
    vars_id      <- names(x)[!(names(x) %in% vars_varying)]
    
    dL <- reshape(x,
                  direction="long",
                  idvar=vars_id,
                  varying=vars_varying,
                  v.names="observed",
                  timevar="statistic")
    
    rownames(dL) <- NULL
    
    dL[["statistic"]] <- factor(dL[["statistic"]],
                                levels=seq_along(vars_varying),
                                labels=vars_varying)

    dL
}
