#------------------------------------------------------------------------------#
#' Spatial Analysis by Distance IndicEs (SADIE).
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
#' my_incidence <- clump(my_incidence, unit_size = c(x = 6, y = 6))
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
sadie <- function(data, ...) UseMethod("sadie")

#------------------------------------------------------------------------------#
#' @rdname sadie
#' @method sadie data.frame
#' @export
#------------------------------------------------------------------------------#
sadie.data.frame <- function(data, index = c("Perry", "Li-Madden-Xu", "all"),
                             nperm = 100, rseed = TRUE, seed = 12345,
                             threads = 1, ..., method = "shortsimplex") {

    index <- match.arg(index)
    # data structure:
    # - 1st and 2nd columns: x and y coordinates, respectively.
    # - 3rd column: observed disease intensity data.
    stopifnot(ncol(data) == 3)
    colnames(data) <- c("x", "y", "i")
    # ^ Needed for ggplot2 figures and to simplify the code below.
    cost <- as.matrix(dist(data[c("x", "y")]))

    # More interesting
    if (!rseed) set.seed(seed)

    N <- nrow(data)
    data[["i"]] <- as.numeric(data[["i"]]) # Just in case it's an integer, to have the same nature for data and flat_data

    # Create a homogeneous (flatten) version of data
    flat_data <- rep(mean(data[["i"]]), N)

    # Use optimal transportation algorithm
    opt_transport <- wrap_transport(data[["i"]], flat_data, cost, method)
    opt_transport <- as.matrix(opt_transport, N)

    # Computation the costs
    # Due to the convention only positive values (donor units) are added up
    # to compute Da. Da = total observed distance (or total cost).
    cost_flows <- costTotCPP(opt_transport, cost)
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
                          start = data[["i"]], end = flat_data,
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
                              start = data[["i"]], end = flat_data,
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

    # Returns:
    res <- list(info_clust = info_clust,
                Da = Da, # Da = total observed distance (or total cost).
                Ia = Ia,
                Pa = Pa,
                Ea = Ea,
                summary_idx = summary_idx,
                nperm = nperm,
                rseed = rseed,
                seed = NA_real_) #TODO: to be corrected
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
    sadie.data.frame(mapped_data, index, nperm, rseed, seed, threads, ...,
                     method = method)
}


#==============================================================================#
# Utilities
#==============================================================================#

#------------------------------------------------------------------------------#
# Wrapper for the transport::transport function.
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
#------------------------------------------------------------------------------#
as.matrix.transport <- function(x, dim_mat) {
    as_matrix_transport(x, dim_mat)
}

#------------------------------------------------------------------------------#
# Computation of clustering indices from Perry et al. (1999).
#------------------------------------------------------------------------------#
clust_P <- function(flow, cost, dim_mat, nperm, start, end,
                    method = "shortsimplex", cost_flows, threads) {

    # Randomization procedure:
    idx <- seq_len(nrow(flow)) # or: seq_len(length(start))
    randomizations <- pbapply::pblapply(seq_len(nperm), function(i1) {
        rand_idx       <- sample(idx)
        new_start      <- start[rand_idx]
        opt_transport  <- wrap_transport(new_start, end, cost, method)
        opt_transport  <- as.matrix(opt_transport, dim_mat)
        idx_same_count <- lapply(idx, function(i2) which(new_start == start[i2]))
        res_Yci <- vapply(seq_len(length(idx_same_count)), function(i2) {
            mean(vapply(seq_len(length(idx_same_count[[i2]])),
                function(x) {
                    costTotiCPP(idx_same_count[[i2]][x], opt_transport, cost)
                }, numeric(1L)))
        }, numeric(1L))
        res_Yii <- vapply(idx,
            function(x) {
                costTotiCPP(x, opt_transport, cost)
            }, numeric(1L))
        # Due to the convention, only positive values (donor units) are added up
        # to compute Ea.
        cost_flows_perm <- costTotCPP(opt_transport, cost)
        Ea <- sum(cost_flows_perm[cost_flows_perm >= 0]) # Donors
        list(flow = opt_transport,
             idx  = rand_idx,
             Yci  = res_Yci,
             Yii  = res_Yii,
             cost_flows_perm = Ea)
    }, cl = threads)

    # Compute indices:
    Ycs <- simplify2array(lapply(randomizations, function(x) x[["Yci"]]))
    Ycs <- rowMeans(abs(Ycs))
    Yis <- simplify2array(lapply(randomizations, function(x) x[["Yii"]]))
    Yis <- rowMeans(abs(Yis))
    Y0  <- mean(Ycs) # Must be equal to: mean(Yis)
    Yi  <- vapply(idx,
        function(x) {
            costTotiCPP(x, flow, cost)
        }, numeric(1L))
    clust <- ((Yi * Y0) / (Yis * Ycs))
    Ea <- unlist(lapply(randomizations, function(x) x[["cost_flows_perm"]]))

    # Returns:
    list(clust = clust,
         Ea    = Ea)
}

#------------------------------------------------------------------------------#
# Computation of clustering indices from Li, Madden and Xu (2012).
#------------------------------------------------------------------------------#
clust_LMX <- function(flow, cost, dim_mat, nperm, start, end,
                      method = "shortsimplex", cost_flows, threads) {

    # Randomization procedure:
    idx <- seq_len(nrow(flow)) # or: seq_len(length(start))
    randomizations <- pbapply::pblapply(seq_len(nperm), function(i1) {
        vapply(idx, function(i2) {
            sub_idx       <- idx[idx != i2]
            rand_idx      <- sample(sub_idx)
            rand_idx      <- append(rand_idx, i2, after = i2 - 1)
            new_start     <- start[rand_idx]
            opt_transport <- wrap_transport(new_start, end, cost, method)
            opt_transport <- as.matrix(opt_transport, dim_mat)
            costTotiCPP(i2, opt_transport, cost)
        }, numeric(1L))
    }, cl = threads)

    # Compute indices:
    randomizations <- simplify2array(randomizations)
    Dis_bar        <- rowMeans(cbind(randomizations, cost_flows))
    clust   <- cost_flows / abs(Dis_bar)
    clust_i <- randomizations / abs(Dis_bar)
    prob    <- vapply(seq_len(length(clust)), function(i1) {
        if (clust[i1] >= 0) rpos <- (clust_i[i1, ] > clust[i1]) # Donor
        else                rpos <- (clust_i[i1, ] < clust[i1]) # Receiver
        rpos <- sum(rpos) + 1
        # Trick above:
        # (1) sum(T, T, F) gives 2.
        # (2) +1 b/c clust is +1 away from the nearest lower or upper clust_i.
        rpos / (nperm + 1)
    }, numeric(1L))

    # Returns:
    list(clust   = clust,
         prob    = prob)
}


#==============================================================================#
# Print, summary and plot
#==============================================================================#

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.sadie <- function(x, ...) {
    cat("Spatial Analysis by Distance IndicEs (sadie)\n")
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
#' @export
#------------------------------------------------------------------------------#
plot.sadie <- function(x, ..., index = c("Perry", "Li-Madden-Xu"),
                       isoclines = FALSE, resolution = 100,
                       thresholds = c(-1.5, 1.5),
                       point_size = c("radius", "area")) {

    index      <- match.arg(index)
    point_size <- match.arg(point_size)

    idx <- switch(index, "Perry" = "idx_P", "Li-Madden-Xu" = "idx_LMX")
    data_clust <- x$info_clust
    data_loess <- stats::loess(as.formula(paste0(idx, " ~ x * y")), data = data_clust, degree = 2, span = 0.2) # TODO: No loess if not required
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
















