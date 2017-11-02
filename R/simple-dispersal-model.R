#------------------------------------------------------------------------------#
#' Simple dispersal model
#'
#' This model is used to perform some assessments and test statistical methods
#' implemented in the \code{epiphy} package.
#'
#' @param nfoci The number of (initial) infections.
#' @param xrate Multiplication rate, the number of new infections produced by
#'     each pre-existing infection at each generation.
#' @param ngen The number of generations.
#' @param lambda The rate parameter for the exponential distribution.
#'
#' @references
#'
#' Esker PD, Sparks AH, Antony G, Bates M, Dall' Acqua W, Frank EE, Huebel L,
#' Segovia V, Garrett KA, 2007. Ecology and Epidemiology in R: Modeling
#' dispersal gradients. The Plant Health Instructor.
#' \href{http://dx.doi.org/10.1094/PHI-A-2007-1226-03}{doi:10.1094/PHI-A-2007-1226-03}
#'
#' @examples
#' # To use the "hiden" functions of this simple model:
#' list2env(epiphy:::simple_model, environment())
#'
#' set.seed(12345)
#' foci <- disperse(nfoci = 10, xrate = 30, ngen = 2, lambda = 50)
#' quad <- quadrat(surf_dim = c(x = 1, y = 1), nint = c(x = 90, y = 90))
#' collection <- collect(foci, quad)
#'
#' # To perform analyses:
#' collection$n <- 1
#' my_data <- incidence(collection, mapping(r = incidence))
#' my_data <- clump(my_data, unit_size = c(x = 3, y = 3)) # /!\ avec clump(...), on retire les colonnes non mappées
#' res <- fit_distr(my_data)
#' res
#' summary(res)
#' plot(res)
#'
#' # Multi-data set analyses
#' set.seed(12345)
#' quad <- quadrat(surf_dim = c(x = 1, y = 1), nint = c(x = 90, y = 90))
#' my_data <- list()
#' for (i in 1:30) {
#'     nfoci <- sample(1:100, 1)
#'     foci <- disperse(nfoci = nfoci, xrate = 800, ngen = 1, lambda = 50)
#'     my_data[[i]] <- collect(foci, quad)
#'     print(paste0("set ", i, " done (nfoci = ", nfoci, ")"))
#' }
#' # saveRDS(my_data, "simple_model_data.rds")
#'
#' ## Power law
#' my_data2 <- lapply(my_data, function(x) {
#'     x$n <- 1
#'     x <- incidence(x, mapping(r = incidence))
#'     clump(x, unit_size = c(x = 3, y = 3))
#' })
#'
#' res <- power_law(my_data2)
#' res
#' summary(res)
#' plot(res)
#'
#' ## Spatial hierarchy
#' low <- my_data2
#' high <- lapply(my_data2, level_up)
#'
#' res <- spatial_hier(low, high)
#' res
#' summary(res)
#' plot(res)
#'
#' @keywords internal
#------------------------------------------------------------------------------#
simple_model <- list()

simple_model$disperse <- function(nfoci, xrate, ngen, lambda){
    foci <- matrix(runif(n = 2*nfoci, min = 0, max = 1),
                   nrow = nfoci, ncol = 2,
                   dimnames = list(NULL, c("x", "y")))
    foci <- cbind(foci, gen = 0)
    lapply(seq_len(ngen), function(i1) {
        lapply(seq_len(nfoci), function(i2) {
            focus   <- foci[i2, ]
            newdir  <- runif(n = xrate, min = 0, max = 2*pi)
            newdist <- rexp(n = xrate, rate = lambda)
            newfoci <- cbind(x   = focus[["x"]] + newdist * cos(newdir),
                             y   = focus[["y"]] + newdist * sin(newdir),
                             gen = i1)
            foci <<- rbind(foci, newfoci)
        })
        nfoci <<- nrow(foci)
    })
    structure(as.data.frame(foci), class = c("disperse", "data.frame"))
}

simple_model$quadrat <- function(surf_dim = c(x = 1, y = 1),
                                    nint = c(x = 90, y = 90)) {
    xquad <- surf_dim[["x"]] / nint[["x"]]
    yquad <- surf_dim[["y"]] / nint[["y"]]

    prex <- head(seq(0, surf_dim[["x"]], by = xquad), -1)
    prey <- head(seq(0, surf_dim[["y"]], by = yquad), -1)

    prex <- cbind.data.frame(x = seq_len(length(prex)), x1 = prex)
    prey <- cbind.data.frame(y = seq_len(length(prey)), y1 = prey)

    quad <- expand.grid(x1 = prex[["x1"]], y1 = prey[["y1"]])
    quad <- merge(merge(quad, prex, by = "x1"), prey, by = "y1")
    quad <- quad[c("x", "y", "x1", "y1")]
    quad <- quad[with(quad, order(x, y)), ]
    quad <- cbind(quad, x2 = quad$x1 + xquad, y2 = quad$y1 + yquad)
    structure(quad, class = c("quadrat", "data.frame"))
}

simple_model$collect <- function(disperse, quadrat) {
    quadrat <- do.call(rbind, lapply(seq_len(nrow(quadrat)), function(i1) {
        quad <- quadrat[i1, ]
        inx <- ((quad$x1 < disperse$x) & (disperse$x <= quad$x2))
        iny <- ((quad$y1 < disperse$y) & (disperse$y <= quad$y2))
        quad$count <- sum(inx * iny) # Nice trick!
        quad$incidence <- ifelse(quad$count == 0, 0, 1)
        quad
    }))
    structure(quadrat, class = c("collect", "quadrat", "data.frame"))
}

# Figures
#simple_model$plot.disperse <- function(x, ...) {
#    plot(x, pch = 19, col = "red", xlim = c(0, 1), ylim = c(0, 1), ...)
#}

#simple_model$plot.quadrat <- ...

# simple_model$plot.collect <- function(x, ..., xlim = c(0, 1), ylim = xlim,
#                                       asp = 1) {
#     #if (!add) {
#         plot(0, type = "n", xlim = xlim, ylim = ylim, xlab = "x", ylab = "y",
#              asp = asp, ...)
#     #}
#     x$coords <- rbind(x$coords , x$coords[1, ]) # to close the box
#     x$col    <- pal(length(0:round(maxi)))[round(x[[type]])]
#     polygon(x$coords, border = "grey", col = x$col)
# }




