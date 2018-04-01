#------------------------------------------------------------------------------#
#' Spatial Analysis by Distance IndicEs (SADIE)
#'
#' \code{sadie} performs the SADIE procedure. It computes different indices and
#' probabilities based on the distance to regularity for the observed spatial
#' pattern and a specified number of random permutations of this pattern. Both
#' kind of clustering indices described by Perry et al. (1999) and Li et al.
#' (2012) can be computed.
#'
#' By convention in the SADIE procedure, clustering indices for a donor unit
#' (outflow) and a receiver unit (inflow) are positive and negative in sign,
#' respectively.
#'
#' @param data A data frame or a matrix with only three columns: the two first
#'     ones must be the x and y coordinates of the sampling units, and the last
#'     one, the corresponding disease intensity observations. It can also be a
#'     \code{\link{count}} or an \code{\link{incidence}} object.
#' @param index The index to be calculated: "Perry", "Li-Madden-Xu" or "all".
#'     By default, only Perry's index is computed for each sampling unit.
#' @param nperm Number of random permutations to assess probabilities.
#' @param rseed Is a random seed used for the permutations? (only useful for
#'     checking purposes).
#' @param seed Seed used when \code{rseed = FALSE} (only useful for checking
#'     purposes).
#' @param threads Number of threads to perform the computations.
#' @param ... Not yet implemented.
#' @param method Method for the transportation algorithm.
#'
#' @references
#'
#' Perry JN. 1995. Spatial analysis by distance indices. Journal of Animal
#' Ecology 64, 303–314. \href{http://dx.doi.org/10.2307/5892}{doi:10.2307/5892}
#'
#' Perry JN, Winder L, Holland JM, Alston RD. 1999. Red–blue plots for detecting
#' clusters in count data. Ecology Letters 2, 106–113.
#' \href{http://dx.doi.org/10.1046/j.1461-0248.1999.22057.x}{doi:10.1046/j.1461-0248.1999.22057.x}
#'
#' Li B, Madden LV, Xu X. 2012. Spatial analysis by distance indices: an
#' alternative local clustering index for studying spatial patterns. Methods in
#' Ecology and Evolution 3, 368–377.
#' \href{http://dx.doi.org/10.1111/j.2041-210X.2011.00165.x}{doi:10.1111/j.2041-210X.2011.00165.x}
#'
#' @examples
#' my_count <- count(aphids, mapping(x = xm, y = ym))
#' my_res <- sadie(my_count)
#' my_res
#' summary(my_res)
#' plot(my_res)
#'
#' my_df <- aphids[, c("xm", "ym", "i")]
#' sadie(my_df)
#'
#' my_incidence <- incidence(tomato_tswv$field_1929[tomato_tswv$field_1929$t == 1, ])
#' my_incidence <- clump(my_incidence, unit_size = c(x = 3, y = 3))
#' plot(my_incidence)
#' my_res <- sadie(my_incidence, index = "all", threads = 4)
#' my_res
#' summary(my_res)
#' plot(my_res) # Identical to: plot(my_res, index = "Perry")
#' plot(my_res, index = "Li-Madden-Xu")
#'
#' @name sadie
#' @export
#------------------------------------------------------------------------------#
sadie <- function(data, ...) {
    UseMethod("sadie")
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @method sadie data.frame
#' @export
#------------------------------------------------------------------------------#
sadie.data.frame <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                             nperm = 100, rseed = TRUE, seed = 12345,
                             threads = 1, ..., method = "shortsimplex") {

    index <- match.arg(index)
    n_col <- ncol(data)
    # data structure:
    # - 1st and 2nd columns: x and y coordinates, respectively.
    # - 3rd column: observed disease intensity data.
    stopifnot(n_col == 3)
    cost <- as.matrix(dist(data[1:2]))

    # More interesting
    if (!rseed) set.seed(seed)

    N <- nrow(data)
    data[[n_col]] <- as.numeric(data[[n_col]]) # Just in case it's integer, to have the same nature for data and fdata

    # Create a homogeneous (flatten) version of data
    fcount <- rep(mean(data[[n_col]]), N)

    # Use optimal transportation algorithm
    opt_transport <- wrap_transport(data[[n_col]], fcount, cost, method)
    opt_transport <- as.matrix(opt_transport, N)

    # Computation the costs
    cost_flows <- vapply(1:N, function(x) {
        costTotiCPP(x, opt_transport, cost, type = "both")
    }, numeric(1L))
    # Due to the convention only positive values are added
    # up. Da = total observed distance (or total cost).
    Da <- sum(cost_flows[cost_flows >= 0]) # Donors

    # Deal with indices
    if (index == "all") index <- c("Perry", "Li-Madden-Xu")
    info_P   <- list(clust = NA_real_, Ea = NA_real_)
    info_LMX <- list(clust = NA_real_, Dis_all = matrix(NA_real_))
    idx_P    <- NA_real_
    idx_LMX  <- NA_real_
    Ea       <- NA_real_ # Ea = randomized distances
    Ia       <- NA_real_ # Da = observed,
    Pa       <- NA_real_
    new_prob <- NA_real_
    Dis_all  <- NA_real_
    if (any(index == "Perry")) {
        cat("Computation of Perry's indices:\n")
        info_P <- clust_P(opt_transport, cost, N, nperm,
                          start = data[[n_col]], end = fcount,
                          method, cost_flows, threads)
        idx_P <- info_P[["clust"]]
        Ea    <- info_P[["Ea"]]
        Ia    <- Da / mean(Ea)
        Pa    <- sum(Ea > Da) / length(Ea)
    } else {
        warning("Ia and Pa cannot be computed if Perry's indices are not.")
    }
    if (any(index == "Li-Madden-Xu")) {
        cat("Computation of Li-Madden-Xu's indices:\n")
        info_LMX <- clust_LMX(opt_transport, cost, N, nperm,
                              start = data[[n_col]], end = fcount,
                              method, cost_flows, threads)
        idx_LMX  <- info_LMX[["clust"]]
        new_prob <- info_LMX[["prob"]]
        Dis_all  <- info_LMX[["Dis_all"]]
    }

    info_clust <- data.frame(data,
                             cost_flows = cost_flows,
                             idx_P      = idx_P,
                             idx_LMX    = idx_LMX,
                             prob       = new_prob)

    # Summary indices
    summary_idx <- data.frame(overall = c(mean(abs(idx_P)),
                                          mean(abs(idx_LMX))),
                              inflow  = c(mean(idx_P[idx_P < 0]),
                                          mean(idx_LMX[idx_LMX < 0])),
                              outflow = c(mean(idx_P[idx_P > 0]),
                                          mean(idx_LMX[idx_LMX > 0])))
    rownames(summary_idx) <- c("Perry's index", "Li-Madden-Xu's index")

    #    return(Sadie())
    res <- list(info_clust = info_clust,
                Da = Da,
                Ia = Ia,
                Pa = Pa,
                Ea_perry = Ea,
                Ea_li = Dis_all, # Needed ?
                summary_idx = summary_idx,
                nperm = nperm,
                rseed = rseed,
                seed = NA_real_) #TODO: to correct
    attr(res, "class") <- "sadie"
    attr(res, "call")  <- match.call()
    res
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @export
#------------------------------------------------------------------------------#
sadie.matrix <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                         nperm = 100, rseed = TRUE, seed = 12345,
                         threads = 1, ..., method = "shortsimplex") {
    sadie.data.frame(as.data.frame(data), index, nperm, rseed, seed,
                     threads, ..., method = method)
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @export
#------------------------------------------------------------------------------#
sadie.count <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                        nperm = 100, rseed = TRUE, seed = 12345,
                        threads = 1, ..., method = "shortsimplex") {
    mapped_data <- map_data(data)
    mapped_data <- mapped_data[, c("x", "y", "i")] # no t
    sadie.data.frame(mapped_data, index, nperm, rseed, seed, threads, ...,
                     method = method)
}

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @export
#------------------------------------------------------------------------------#
sadie.incidence <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                        nperm = 100, rseed = TRUE, seed = 12345,
                        threads = 1, ..., method = "shortsimplex") {
    mapped_data <- map_data(data)
    mapped_data <- mapped_data[, c("x", "y", "i")] # no t, no n
    #mapped_data[["n"]] <- NULL # n is not used in SADIE procedure.
    sadie.data.frame(mapped_data, index, nperm, rseed, seed, threads, ...,
                     method = method)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.sadie <- function(x, ...) {
    cat("\nCall:\n")
    print(attr(x, "call"))
    cat("\nIa: ", format(x$Ia, digits = 1, nsmall = 4),
        " (Pa = ", format.pval(x$Pa), ")\n\n", sep = "")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
summary.sadie <- function(object, ...) {
    res <- object
    class(res) <- "summary.sadie"
    res
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.summary.sadie <- function(x, ...) {
    cat("\nCall:\n")
    print(attr(x, "call"))
    n <- 6L
    cat("\nFirst ", n, " rows of clustering indices:\n", sep = "") # Cf Image pkg for a nice display
    print(head(x$info_clust, n = n))
    #printCoefmat(x$info_clust)
    #cat("number of permutations: ", object$nperm, "\n",
    #    "random seed: ", object$rseed, "\n", sep = "")
    cat("\nSummary indices:\n")
    print(x$summary_idx)
    cat("\nMain outputs:")
    cat("\nIa: ", format(x$Ia, digits = 1, nsmall = 4),
        " (Pa = ", format.pval(x$Pa), ")\n", sep = "")
    cat("\n'Total cost': ", x$Da, "\n",
        "Number of permutations: ", x$nperm, "\n",
        "Random Seed?: ", x$rseed, "(", x$seed, ")\n\n", sep = "")
}

#------------------------------------------------------------------------------#
# Utilities
#------------------------------------------------------------------------------#
wrap_transport <- function(start, end, cost, method = "shortsimplex",
                           control = list(), ...) {
    res <- transport::transport(start, end, costm = cost, method = method,
                                control = control, ...)
    class(res) <- c("transport", class(res))
    res
}

#------------------------------------------------------------------------------#
#' @method as.matrix transport
#' @export
#------------------------------------------------------------------------------#
as.matrix.transport <- function(x, dim_mat) {#, dimnames = NULL, ...) { # N: Number of sampling units... to CPP to optimize the code?
    as_matrix_transport(x, dim_mat)
}

#------------------------------------------------------------------------------#
# standardised and dimensionless clustering index (nui) for a donor
# Perry
#------------------------------------------------------------------------------#
clust_P <- function(flow, cost, dim_mat, nperm, start, end,
                    method = "shortsimplex", cost_flows, threads) {
    idx <- seq_len(nrow(flow)) # or: seq_len(length(start))
    randomisations <- pbapply::pblapply(seq_len(nperm), function(i1) {
        # Yi & Yc
        rand_idx      <- sample(idx)
        new_start     <- start[rand_idx] # new_start n'a pas de nom de col et row contrairement à start (????)
        opt_transport <- wrap_transport(new_start, end, cost, method)
        opt_transport <- as.matrix(opt_transport, dim_mat) # Replace n with N

        # ~ Assez long ci-dessous
        idx_same_count <- lapply(idx, function(i2) which(new_start == start[i2]))
        res_Yci <- vapply(seq_len(length(idx_same_count)), function(i2) {
            mean(vapply(seq_len(length(idx_same_count[[i2]])),
                        function(x) {
                            costTotiCPP(idx_same_count[[i2]][x], opt_transport, cost, "both")
                        },
                        numeric(1L)))
        }, numeric(1L))

        res_Yii <- vapply(idx,
                          function(x) {
                              costTotiCPP(x, opt_transport, cost, "both")
                          },
                          numeric(1L))

        # TODO: To improve below:
        N <- nrow(opt_transport) # TODO: Should be the same as length(start)
        costTot <- vapply(1:N, function(x) {
            costTotiCPP(x, opt_transport, cost, type = "both")
        }, numeric(1L))
        # Due to the convention only positive values are added
        # up. Da = total observed distance (or total cost).
        Ea <- sum(costTot[costTot >= 0]) # Donors

        list(flow = opt_transport,
             idx  = rand_idx,
             Yci = res_Yci,
             Yii = res_Yii,
             #costTot = costTot(opt_transport, cost))) # Pour calculer Ia
             costTot = Ea) #costTotCPP(opt_transport, cost)) # To compute Ia //TODO: Useful costTotCPP???
    }, cl = threads)

    Ycs <- simplify2array(lapply(randomisations, function(x) x[["Yci"]]))
    Ycs <- rowMeans(abs(Ycs))
    Yis <- simplify2array(lapply(randomisations, function(x) x[["Yii"]]))
    Yis <- rowMeans(abs(Yis))
    Y0  <- mean(Ycs) # Should be equal to: mean(Yis)
    Yi  <- vapply(idx,
                  function(x) {
                      costTotiCPP(x, flow, cost, "both")
                  },
                  numeric(1L))

    clust <- ((Yi * Y0) / (Yis * Ycs))
    Ea <- unlist(lapply(randomisations, function(x) x[["costTot"]]))

    list(clust = clust,
         Ea = Ea)
}

#------------------------------------------------------------------------------#
# standardised and dimensionless clustering index (nui) for a donor
#------------------------------------------------------------------------------#
clust_LMX <- function(flow, cost, dim_mat, nperm, start, end,
                             method = "shortsimplex", cost_flows,
                             threads) {
    idx            <- seq_len(nrow(flow)) # or: seq_len(length(start))
    randomisations <- pbapply::pblapply(seq_len(nperm), function(i1) {
        vapply(idx, function(i2) {
            sub_idx   <- idx[idx != i2] ## NEW
            rand_idx  <- sample(sub_idx)
            ##TO DOUBLE CHECK## rand_idx <- insert(rand_idx, i2, i2) ##NEW : insert is in pkgg R.utils à récupérer uniquement fn intéressante en cpp après
            rand_idx  <- append(rand_idx, i2, i2 - 1) # TODO: .insert(i, x)    Insert x at the i^th^ position of, grows vector in RCPP // Attention, at the positin VS after the position
            new_start <- start[rand_idx]
            #dimnames(new_start) <- list(rownames(start), colnames(start))
            opt_transport <- wrap_transport(new_start, end, cost, method)
            opt_transport <- as.matrix(opt_transport, dim_mat)
            costTotiCPP(i2, opt_transport, cost, "both")
        }, numeric(1L))
    }, cl = threads)

    randomisations <- simplify2array(randomisations)
    Dis_bar        <- rowMeans(cbind(randomisations, cost_flows))  # TODO : to check : # Négatif tous pour l'instant. cost_flows = Dis ; il y a (nperm + 1) elements

    clust   <- cost_flows / abs(Dis_bar) # sing_ ... allows to make a donor positive and a receiver negative
    clust_i <- randomisations / abs(Dis_bar)
    prob    <- vapply(seq_len(length(clust)), function(i1) {
        if (clust[i1] >= 0) rpos <- (clust_i[i1, ] > clust[i1]) # Donor
        else                rpos <- (clust_i[i1, ] < clust[i1]) # Receiver
        rpos <- sum(rpos) + 1 # Trick : sum(T, T, F) => 2 ; +1 because Di is +1 away from the other sup or inf Di,rep
        rpos / (nperm + 1)
    }, numeric(1L))

    # Returns:
    list(clust   = clust,
         prob    = prob,
         Dis_all = randomisations) # Needed?
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
plot.sadie <- function(x, y, ..., index = c("Perry", "Li-Madden-Xu"),
                       isoclines = FALSE, resolution = 100,
                       thresholds = c(-1.5, 1.5),
                       point_size = c("radius", "area")) {

    index      <- match.arg(index)
    point_size <- match.arg(point_size)

    idx <- switch(index, "Perry" = "idx_P", "Li-Madden-Xu" = "idx_LMX")
    data_clust <- x$info_clust
    data_loess <- stats::loess(as.formula(paste0(idx, " ~ x * y")), data = data_clust, degree = 2, span = 0.2)
    input_val <- expand.grid(
        x = seq(min(data_clust$x), max(data_clust$x), length = resolution),
        y = seq(min(data_clust$y), max(data_clust$y), length = resolution))
    interpolated   <- predict(data_loess, input_val)
    data_landscape <- data.frame(input_val, z = as.vector(interpolated))

    data_clust$col <- vapply(data_clust[[idx]], function(x) {
        if      (x <  thresholds[1L]) "blue"  # Low values
        else if (x <= thresholds[2L]) "white" # Intermediate values
        else                          "red"   # High values
    }, character(1L))

    gg <- ggplot()
    if (isoclines) {
        gg <- gg + geom_raster(data = data_landscape, aes(x, y, fill = z))
        gg <- gg + geom_contour(data = data_landscape, aes(x, y, z = z),
                                size = 0.6, color = "black")
        # Note that it is not possible in current ggplot2 implementation to use
        # two different colour scales in the same figure.
        gg <- gg + geom_point(data = data_clust,
                              aes_(quote(x), quote(y),
                                   size = call("abs", as.name(idx))),
                              colour = "black", fill = "black", pch = 21)
        switch(point_size,
               "area" = {
                   gg <- gg + scale_size("Absolute\nindex", range = c(0, 10))
               },
               "radius" = {
                   gg <- gg + scale_radius("Absolute\nindex", range = c(0, 10))
               }
        )
        gg <- gg + scale_fill_gradientn("Interpolated\nindex",
                                        colours = terrain.colors(10))
    } else {
        gg <- gg + geom_point(data = data_clust,
                              aes_(quote(x), quote(y), fill = quote(col),
                                   size = call("abs", as.name(idx))),
                              colour = "black",  pch = 21)
        switch(point_size,
            "area" = {
                gg <- gg + scale_size("Absolute\nindex", range = c(0, 10))
            },
            "radius" = {
                gg <- gg + scale_radius("Absolute\nindex", range = c(0, 10))
            }
        )
        gg <- gg + scale_fill_manual("Index\nthreshold",
            breaks = c("red", "white", "blue"),
            labels = c("red"   = paste0("> ", thresholds[2L]),
                       "white" = paste0("Between ", thresholds[1L],
                                        "\nand ", thresholds[2L]),
                       "blue" = paste0("< ", thresholds[1L])),
            values = c("red" = "red", "white" = "white", "blue" = "blue"))
    }
    gg <- gg + theme_bw()
    print(gg)
    invisible(NULL)
}
















