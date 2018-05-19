#------------------------------------------------------------------------------#
#' @include intensity-classes.R
#------------------------------------------------------------------------------#
NULL

#==============================================================================#
# Aggregation indices
#==============================================================================#

#------------------------------------------------------------------------------#
#' Several aggregation indices.
#'
#' This function can compute different aggregation indices. See "Details"
#' section for more information about the available indices.
#'
#' There are currently three implemented methods to compute indices of
#' aggregation.
#'
#' \code{fisher}: Fisher's index of aggregation. In case of a count, this index
#' corresponds to the ratio of the observed variance to the observed mean, and
#' this is why this index is also known as the variance to mean ratio. For a
#' binary variable, a similar index can be calculated using instead the ratio of
#' the observed variance to the theoretical variance if data follow a binomial
#' law (i.e. a reference distribution for a random pattern of individuals within
#' sampling units).
#'
#' \code{lloyd}: Lloyd's index of patchiness. The value of this index increases
#' with the degree of aggregation. Note that Lloyd's mean crowding can also be
#' returned if \code{type = "mean-crowding"} is provided as parameter.
#'
#' \code{morisita}: Morisita's coefficient of dispersion. This index can be
#' computed for either count or incidence data, but its interpretation can be
#' uncertain.
#'
#' Values of Fisher's and Lloyd's indices can be interpreted as follows:
#' \itemize{
#'     \item index < 1: uniform pattern;
#'     \item index = 1: random pattern;
#'     \item index > 1: aggregated pattern.
#' }
#'
#' The following table gives information about the applicability of the various
#' methods.
#'
#' \tabular{llll}{
#'     \tab count \tab incidence \tab severity \cr
#'     fisher   \tab + \tab + \tab - \cr
#'     lloyd    \tab + \tab - \tab - \cr
#'     morisita \tab + \tab + \tab - \cr
#' }
#' where + means implemented, and -, not implemented (or not possible). At the
#' moment, there is no index of aggregation for severity data.
#'
#' @param x A numeric vector or a \code{count}/\code{incidence} object.
#' @param method The name of the method to be used. "fisher" method is used by
#'     default. See details below.
#' @param flavor Which flavor of this index must be calculated ("count" or
#'     "incidence")?
#' @param n Number of individuals per sampling unit. If \code{n} is provided,
#'     the "incidence" flavor is calculated whatever the value of \code{flavor}.
#'     Note that current implementation only deals with equal size sampling
#'     units.
#' @param ... Additional arguments to be passed to other methods.
#'
#' @examples
#' # Count flavor of Fisher's index:
#' my_fisher_count <- agg_index(aphids$i)
#' my_fisher_count
#'
#' # And incidence flavor of Fisher's index:
#' my_fisher_incidence <- agg_index(tobacco_viruses$i, n = tobacco_viruses$n)
#' my_fisher_incidence
#'
#' # Either standard R or epiphy idioms can be used:
#' identical(my_fisher_count, agg_index(count(aphids)))
#' identical(my_fisher_incidence, agg_index(incidence(tobacco_viruses)))
#'
#' # Lloyd's index (only for count data):
#' agg_index(aphids$i, method = "lloyd")
#' # Lloyd's mean crowding:
#' agg_index(aphids$i, method = "lloyd", type = "mean-crowding")
#'
#' # Count flavor of Morisita's index:
#' agg_index(aphids$i,  method = "morisita")
#' # Incidence flavor of Morisita's index:
#' agg_index(tobacco_viruses$i, n = tobacco_viruses$n, method = "morisita")
#'
#' @seealso \code{\link[vegan]{vegdist}} in \strong{vegan} package.
#'
#' @references
#'
#' Fisher RA. 1925. Statistical methods for research workers. Oliver and Boyd,
#' Edinburgh.
#'
#' Lloyd M. 1967. Mean crowding. The Journal of Animal Ecology 36, 1–30.
#'
#' Morisita M. 1962. I\eqn{\delta}-Index, a measure of dispersion of
#' individuals. Researches on Population Ecology 4, 1–7.
#' \href{http://dx.doi.org/doi:10.1007/BF02533903}{doi:10.1007/BF02533903}
#'
#' Madden LV, Hughes G. 1995. Plant disease incidence: Distributions,
#' heterogeneity, and temporal analysis. Annual Review of Phytopathology 33(1):
#' 529–564.
#' \href{http://dx.doi.org/doi:10.1146/annurev.py.33.090195.002525}{doi:10.1146/annurev.py.33.090195.002525}
#'
#' @export
#------------------------------------------------------------------------------#
agg_index <- function(x, method = c("fisher", "lloyd", "morisita"),
                      flavor = c("count", "incidence"), n = NULL, ...) {
    call   <- match.call() # TODO: Not yet used.
    flavor <- match.arg(flavor)
    method <- match.arg(method)
    switch (method,
        "fisher"   = fisher(x, flavor = flavor, n = n, ...),
        "lloyd"    = lloyd(x, ...),
        "morisita" = morisita(x, flavor = flavor, n = n, ...),
        stop("No such method.")
    )
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.agg_index <- function(x, ...) {
    cat(x$name, ":\n", sep ="")
    if (!is.null(x$flavor)) {
        cat("(Version for ", x$flavor, " data)\n", sep = "")
    }
    cat(formatC(x$index), "\n", sep = "")
}

#------------------------------------------------------------------------------#
# Fisher's index of aggregation.
#------------------------------------------------------------------------------#
fisher <- function(x, ...) UseMethod("fisher")

#------------------------------------------------------------------------------#
fisher.default <- function(x, flavor = c("count", "incidence"), n = NULL, ...) {
    # TODO: Think about the revisited version (section 9.4.6, page 244, Madden
    # et al, 2007), Incidence: D = s_y^2 / s_{bin}^2
    call <- match.call() # TODO: Not yet used.
    stopifnot(is.numeric(x))
    flavor <- match.arg(flavor)
    if (!is.null(n) && !is.na(n)) {
        # Force "incidence" flavor if n is provided.
        flavor <- "incidence"
    }
    x <- na.omit(x)
    N <- length(x) # Number of sampling units.
    switch (flavor,
        "count" = {
            m <- mean(x)
            v <- var(x)
            res <- v / m
        },
        "incidence" = {
            stopifnot(!is.null(n) && !is.na(n))
            n <- unique(n)
            if (length(n) != 1) stop(paste0("Current implementation only deals ",
                                            "with equal size sampling units."))
            if (!(n > 1)) stop(paste0("The number of individuals per sampling ",
                                      " unit ('n') must be > 1."))
            stopifnot(all(x <= n))
            m <- mean(x / n)
            v <- var(x / n)
            # in terms of proportions:
            res <- v / (m * (1 - m) / n)
            # in terms of numbers of diseased individuals per sampling unit:
            #res <- var(d) / (n * m * (1 - m)))
            # TODO: to check (Madden et al 2007, Turechek et al 2011)...
            # Anyway, I prefer "in terms of proportions" because p and v are
            # homogeneous in this case.
        }
    )
    structure(list(index  = res,
                   name   = "Fisher's index of dispersion",
                   flavor = flavor,
                   N      = N,
                   n      = n),
              class = c("fisher", "agg_index"))
}

#------------------------------------------------------------------------------#
fisher.count <- function(x, ...) {
    fisher.default(map_data(x)[["i"]], flavor = "count")
}

#------------------------------------------------------------------------------#
fisher.incidence <- function(x, ...) {
    mapped_data <- map_data(x)
    fisher.default(mapped_data[["i"]], flavor = "incidence",
                   n = mapped_data[["n"]])
}

#------------------------------------------------------------------------------#
# Lloyd's index of patchiness.
#------------------------------------------------------------------------------#
lloyd <- function(x, type, ...) UseMethod("lloyd")

#------------------------------------------------------------------------------#
lloyd.default <- function(x, type = c("patchiness", "mean-crowding"), ...) {
    call <- match.call() # TODO: Not yet used.
    stopifnot(is.numeric(x))
    type <- match.arg(type)
    x <- na.omit(x)
    m <- mean(x)
    v <- var(x)
    mc <- m + ((v / m) - 1)
    switch (type,
            "patchiness"    = {
                name <- "Lloyd's index of patchiness"
                res  <- mc / m
            },
            "mean-crowding" = {
                name <- "Lloyd's mean crowding"
                res  <- mc
            }
    )
    structure(list(index  = res,
                   name   = name,
                   flavor = NULL,
                   N      = length(x), # Number of sampling units.
                   n      = NULL),
              class = c("lloyd", "agg_index"))
}

#------------------------------------------------------------------------------#
lloyd.count <- function(x, type = c("patchiness", "mean-crowding"), ...) {
    lloyd.default(map_data(x)[["i"]], type = type, ...)
}

#------------------------------------------------------------------------------#
# Morisita's coefficient of dispersion.
#------------------------------------------------------------------------------#
morisita <- function(x, ...) UseMethod("morisita")

#------------------------------------------------------------------------------#
morisita.default <- function(x, flavor = c("count", "incidence"), n = NULL, ...) {
    call <- match.call() # TODO: Not yet used.
    stopifnot(is.numeric(x))
    flavor <- match.arg(flavor)
    if (!is.null(n) && !is.na(n)) {
        # Force "incidence" flavor if n is provided.
        flavor <- "incidence"
    }
    x <- na.omit(x)
    N <- length(x) # Number of sampling units.
    tot <- sum(x)  # Total number of individuals over all the sampling units.
    res <- N * sum(x * (x - 1)) / (tot * (tot - 1))
    switch (flavor,
            "count" = {
                # Do no more calculation.
            },
            "incidence" = {
                stopifnot(!is.null(n) && !is.na(n))
                n <- unique(n)
                if (length(n) != 1) stop(paste0("Current implementation only deals ",
                                                "with equal size sampling units."))
                if (!(n > 1)) stop(paste0("The number of individuals per sampling ",
                                          " unit ('n') must be > 1."))
                stopifnot(all(x <= n))
                res <- res * (n - 1/N) / (n - 1)
            }
    )
    structure(list(index  = res,
                   name   = "Morisita's coefficient of dispersion",
                   flavor = flavor,
                   N      = N,
                   n      = n),
              class = c("morisita", "agg_index"))
}

#------------------------------------------------------------------------------#
morisita.count <- function(x, ...) {
    morisita.default(map_data(x)[["i"]], flavor = "count")
}

#------------------------------------------------------------------------------#
morisita.incidence <- function(x, ...) {
    mapped_data <- map_data(x)
    morisita.default(mapped_data[["i"]], flavor = "incidence",
                     n = mapped_data[["n"]])
}

#==============================================================================#
# Stat tests for aggregation indices
#==============================================================================#

#------------------------------------------------------------------------------#
#' Chi-squared test.
#'
#' Performs chi-squared tests for Fisher's aggregation indices (computed with
#' either count or incidence data). If another kind of data is provided, the R
#' standard \code{chisq.test} function is called.
#'
#' Under the null hypothesis for Fisher's aggregation index (index = 1, i.e. a
#' random pattern is observed), (N - 1)*index follows a chi-squared distribution
#' with N - 1 degrees of freedom. N is the number of sampling units.
#'
#' @param x Either the output of the \code{\link{agg_index}} function with
#'     \code{method = "fisher"} as parameter, or another R object. In the latter
#'     case, stats::\code{\link[stats]{chisq.test}} is called.
#' @param ... Further arguments to be passed to
#'     stats::\code{\link[stats]{chisq.test}}.
#'
#' @examples
#' # For incidence data:
#' my_incidence <- incidence(tobacco_viruses)
#' my_fisher <- agg_index(my_incidence, method = "fisher")
#' chisq.test(my_fisher)
#'
#' @seealso \code{\link{calpha.test}}, \code{\link{z.test}}
#'
#' @references
#'
#' For count and incidence data:
#'
#' Madden LV, Hughes G. 1995. Plant disease incidence: Distributions,
#' heterogeneity, and temporal analysis. Annual Review of Phytopathology 33(1):
#' 529–564.
#' \href{http://dx.doi.org/doi:10.1146/annurev.py.33.090195.002525}{doi:10.1146/annurev.py.33.090195.002525}
#'
#' Patil GP, Stiteler WM. 1973. Concepts of aggregation and their
#' quantification: a critical review with some new results and applications.
#' Researches on Population Ecology, 15(1): 238-254.
#'
#' @export
#------------------------------------------------------------------------------#
chisq.test <- function(x, ...) UseMethod("chisq.test")

#------------------------------------------------------------------------------#
#' @rdname chisq.test
#' @method chisq.test default
#' @export
#------------------------------------------------------------------------------#
chisq.test.default <- function(x, ...) stats::chisq.test(x, ...)

#------------------------------------------------------------------------------#
#' @rdname chisq.test
#' @method chisq.test fisher
#' @export
#------------------------------------------------------------------------------#
chisq.test.fisher <- function(x, ...) {
    call  <- match.call() # TODO: Not yet used.
    dname <- deparse(substitute(x))
    N     <- x[["N"]]
    df    <- N - 1
    stat  <- df * x[["index"]]
    pval  <- pchisq(stat, df, lower.tail = FALSE)
    structure(list(
        statistic = c("X-squared" = stat),
        parameter = c("df" = df),
        p.value   = pval,
        method    = paste0("Chi-squared test for (N - 1)*index ",
                           "following a chi-squared distribution (df = N - 1)"),
        data.name = dname
        #observed =,
        #expected =,
        #residuals =,
        #stdres =
    ), class = "htest")
}

#------------------------------------------------------------------------------#
#' Z-test.
#'
#' Performs z-tests for Fisher's aggregation indices (computed with either count
#' or incidence data).
#'
#' For two-sided tests with a confidence level of 95%, if -1.96 <= z <= 1.95,
#' the spatial pattern would be random. If z < -1.96 or z > 1.96, it would be
#' uniform or aggregated, respectively.
#'
#' @param x The output of the \code{\link{agg_index}} function with
#'     \code{method = "fisher"} as parameter.
#' @param alternative A character string specifying the alternative hypothesis.
#'     It must be one of "two.sided" (default), "less" or "greater".
#' @param conf.level The confidence level of the interval.
#' @param ... Not yet implemented.
#'
#' @examples
#' # For incidence data:
#' my_incidence <- incidence(tobacco_viruses)
#' my_fisher <- agg_index(my_incidence, method = "fisher")
#' z.test(my_fisher)
#'
#' @seealso \code{\link{calpha.test}}, \code{\link{chisq.test}}
#'
#' @references
#'
#' For count and incidence data:
#'
#' Moradi-Vajargah M, Golizadeh A, Rafiee-Dastjerdi H, Zalucki MP, Hassanpour M,
#' Naseri B. 2011. Population density and spatial distribution pattern of Hypera
#' postica (Coleoptera: Curculionidae) in Ardabil, Iran. Notulae Botanicae Horti
#' Agrobotanici Cluj-Napoca, 39(2): 42-48.
#'
#' Sun P, Madden LV. 1997. Using a normal approximation to test for the binomial
#' distribution. Biometrical journal, 39(5): 533-544.
#'
#' @export
#------------------------------------------------------------------------------#
z.test <- function(x, ...) UseMethod("z.test")

#------------------------------------------------------------------------------#
#' @rdname z.test
#' @method z.test default
#' @export
#------------------------------------------------------------------------------#
z.test.default <- function(x, ...) {
    # TODO:
    # if (!requireNamespace("snse", quietly = TRUE)) {
    #     stop("snse needed for this function to work. Please install it.",
    #          call. = FALSE)
    # }
    # snse::z.test(x, ...)
    stop("No method implemented for this kind of data.")
}

#------------------------------------------------------------------------------#
#' @rdname z.test
#' @method z.test fisher
#' @export
#------------------------------------------------------------------------------#
z.test.fisher <- function(x, alternative = c("two.sided", "less", "greater"),
                          conf.level = 0.95, ...) {
    call        <- match.call() # TODO: Not yet used.
    dname       <- deparse(substitute(x))
    alternative <- match.arg(alternative)
    N           <- x[["N"]]
    switch (x[["flavor"]],
            "count" = {
                # From Moradi-Vajargah et al. (2011)
                # TODO: To double-check. Try to find the original paper!
                stat <- sqrt(2 * (N - 1) * x[["index"]]) - sqrt(2 * (N - 1) - 1)
            },
            "incidence" = {
                # From Sun and Madden (1997):
                n <- x[["n"]]
                stat <- (x[["index"]] - 1) / (sqrt(2 * (1 - 1/n) / (N - 1)))
            }
    )
    switch(alternative,
           "two.sided" = {
               pval <- 2 * pnorm(abs(stat), lower.tail = FALSE)
               # TODO: Confidence interval
               #alpha <- 1 - conf.level
               #cint <- c(stat * stderr - qnorm((1 - alpha/2)) * stderr,
               #          stat * stderr + qnorm((1 - alpha/2)) * stderr)
           },
           "less" = {
               pval <- pnorm(stat)
               # TODO: Confidence interval
               #cint <- c(-Inf, stat * stderr + qnorm(conf.level) * stderr)
           },
           "greater" = {
               pval <- pnorm(stat, lower.tail = FALSE)
               # TODO: Confidence interval
               #cint <- c(stat * stderr - qnorm(conf.level) * stderr, Inf)
           }
    )
    #TODO: attr(cint, "conf.level") <- conf.level
    structure(list(
        statistic = c("z" = stat),
        #parameter =,
        p.value   = pval,
        #TODO: conf.int  = cint,
        #TODO: estimate =,
        #TODO: null.value = 0,
        alternative = alternative,
        method    = "One-sample z-test",
        data.name = dname
        #observed =,
        #expected =,
        #residuals =,
        #stdres =
    ), class = "htest")
}

#------------------------------------------------------------------------------#
#' C(alpha) test.
#'
#' The C(alpha) test is a test of the binomial distribution against the
#' alternative of the beta-binomial distribution.
#'
#' It is based on calculation of a test statistic, z, that has an asymptotic
#' standard normal distribution under the null hypothesis. It is one-sided (in
#' the way that the alternative is aggregation, not just "non-randomness"), thus
#' with a confidence level of 95%, the null hypothesis is rejected when z >
#' 1.64. When all the sampling units contain the same total number of
#' individuals, n, the test statistic is calculated from:
#'
#' z = (n(N - 1)I - Nn)/(2Nn(n - 1))^(1/2)
#'
#' where N is the number of sampling units, and I, Fisher's index of aggregation
#' for incidence data.
#'
#' @param x The output of the \code{\link{agg_index}} function with
#'     \code{method = "fisher"} as parameter.
#' @param ... Not yet implemented.
#'
#' @examples
#' # For incidence data:
#' my_incidence <- incidence(tobacco_viruses)
#' my_fisher <- agg_index(my_incidence, method = "fisher")
#' calpha.test(my_fisher)
#'
#' @seealso \code{\link{chisq.test}}, \code{\link{z.test}}
#'
#' @references
#'
#' Neyman J. 1959. Optimal asymptotic tests of composite statistical hypotheses.
#' In: Probability and Statistics, 213-234. Wiley, New York.
#'
#' Tarone RE. 1979. Testing the goodness of fit of the binomial distribution.
#' Biometrika, 66(3): 585-590.
#'
#' @export
#------------------------------------------------------------------------------#
calpha.test <- function(x, ...) UseMethod("calpha.test")

#------------------------------------------------------------------------------#
#' @rdname calpha.test
#' @method calpha.test fisher
#' @export
#------------------------------------------------------------------------------#
calpha.test.fisher <- function(x, ...) {
    call  <- match.call() # TODO: Not yet used.
    dname <- deparse(substitute(x))
    N     <- x[["N"]]
    switch (x[["flavor"]],
            "count" = {
                stop("No calpha.test for count data.")
            },
            "incidence" = {
                n <- x[["n"]]
                stat <- ((n * (N - 1) * x[["index"]]) - (N * n)) /
                    sqrt(2 * N * n * (n - 1)) # a one-sided test in the sens that...
            }
    )
    pval <- 2 * pnorm(abs(stat), lower.tail = FALSE)
    structure(list(
        statistic = c("z" = stat),
        #parameter =,
        p.value   = pval,
        method    = "C(alpha) test",
        data.name = dname
        #observed =,
        #expected =,
        #residuals =,
        #stdres =
    ), class = "htest")

}





