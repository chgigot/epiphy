#------------------------------------------------------------------------------#
#' @include intensity-classes.R
#------------------------------------------------------------------------------#
NULL


#       x: Community data, a matrix-like object or a vector.
#
#   index: Diversity index, one of ‘"shannon"’, ‘"simpson"’ or
#         ‘"invsimpson"’.
#
# from vegan

## About agg indexes:
## https://en.wikipedia.org/wiki/Taylor's_law
## cf. http://www.inside-r.org/packages/cran/vegan/docs/vegdist


#==============================================================================#
# Aggregation indices (indices of patchiness)
#==============================================================================#
#' Aggregation indices
#'
#' Compute the eponym indices of aggregation.
#'
#' @param data A vector of integers or an \code{Intensity} object.
#'
#' @name indices
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' @describeIn indices Fisher's index of dispersion, D
#'
#' In case of a binary variable, the index is a ratio of two variances: the
#' observed variance and the theoretical variance if data follow a binomial law,
#' i.e. a reference distribution for a random pattern of diseased individuals
#' within sampling units.
#' Under the null hypothesis (D = 1, random pattern), (N - 1)D follows a
#' chi-squared distribution with N - 1 degrees of freedom. N is the number of
#' quadrats.
#'
#' Count (Moradi-Vajargah et al, 2011: CAPS)... cf Gosme for Incidence
#' V/M < 1: means uniform / REGULAR
#' V/M = 1: means RANDOM
#' V/M > 1: means AGGREGATED / clustered
#'
#' @references
#' Fisher RA. 1925. Statistical methods for research workers. Oliver and Boyd, Edinburgh.
#'
#' @export
#------------------------------------------------------------------------------#
fisher1 <- function(data, confLevel = 0.95) { # confLevel non pris en compte pour l'instant !!!

    # Think about the revisited version (section 9.4.6, page 244, Madden et al, 2007), Incidence
    # D = s_y^2 / s_{bin}^2

    stopifnot(is.count(data) || is.incidence(data))
    call <- match.call()
    d <- data$obs$d # All the raw data
    N <- length(d)
    if (is.Count(data)) {
        m <- mean(d, na.rm = TRUE)
        v <- var(d, na.rm = TRUE)
        D <- v / m
    }
    if (is.incidence(data)) {
        n <- unique(data$obs$n)
        if (length(n) != 1) stop(paste0("Current implementation only deals ",
                                        "with equal size sampling units."))
        m <- mean((d / n), na.rm = TRUE)
        v <- var((d / n), na.rm = TRUE)
        # in terms of proportions:
        D <- v / (m * (1 - m) / n)
        # in terms of numbers of diseased individuals per sampling unit:
        #D <- var(d) / (n * m * (1 - m))) ## À vérifier (Madden et al 2007, Turechek et al 2011)... quoi qu'il en soit, je prefere en treme de proportion car p et v sont homogene dans ce cas-là
    }

    ## Different test are possible:
    # - chi2: OK
    # - z_sm: NOT YET
    #     # Count: (cf. Moradi-Vajargah et al, 2011,...)
    #     --> zValue = sqrt(2 * (N - 1) * D) - sqrt(2 * (N - 1) - 1) # To double check
    #     zValue < -1.96 -----------> UNIFORM spatial distribution
    #     -1.96 >= zValue >= 1.96 --> RANDOM spatial distribution
    #     1.96 < zValue ------------> AGGREGATED spatial distribution
    #     # Incidence: (cf. Madden et al, 2007, p. 243)
    #     --> zValue = ((n * (N - 1) * D) - (N * n)) / sqrt(2 * N * n * (n - 1)) # To double check
    # - C(alpha): NOT YET (cf. Gosme..., Turechek, Madden, 1999 - 2 colonnes)
    #
    # "The test statistic has a standard normal distribution under the null
    # hypothesis and is given by zC(α) = [(nN – 1)D – nN]/[2N(n2 – n)] /2.
    # However, the alternative hypothesis is not just overdispersion but
    # overdispersion described by the beta-binomial. This is a one-sided test,
    # thus, the null hypothesis is rejected when zC(α) > 1.645.
    #
    # ref: Tarone (1979)

    # chi-squared test ---> PENSER à créer une fonction
    C  <- (N - 1) * D
    df <- N - 1
    pvalue <- 1 - pchisq(C, df)
    chisq   <- qchisq(pvalue, df)

    chisqTest <- structure(list(
        method = "Pearson's Chi-squared test",
        #estimate = C,### ???????????
        statistic = c("X-squared" = chisq),
        parameter = c("df" = df),
        p.value = pvalue#,
        #data.name = ??,
        #observed = ??,
        #expected = ??,
        #residuals = ??
    ), class = "htest") # A propos "htest": http://www.inside-r.org/node/219017

    structure(list(
        call = call,
        index = D,
        test = chisqTest
    ), class = "aggIndex")
}


##### TESTS

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
fisher <- function(x, ...) UseMethod("fisher")

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
fisher.default <- function(x, flavor = c("count", "incidence"), n = NULL, ...) {
    # Think about the revisited version (section 9.4.6, page 244, Madden et al, 2007), Incidence
    # D = s_y^2 / s_{bin}^2
    stopifnot(is.numeric(x))
    flavor <- match.arg(flavor)
    #call <- match.call()
    x <- na.omit(x)
    N <- length(x)
    switch (flavor,
        "count" = {
            m <- mean(x)
            v <- var(x)
            D <- v / m
        },
        "incidence" = {
            stopifnot(!is.null(n))
            #stopifnot(n > 1)
            m <- mean(x / n)
            v <- var(x / n)
            # in terms of proportions:
            D <- v / (m * (1 - m) / n)
            # in terms of numbers of diseased individuals per sampling unit:
            #D <- var(d) / (n * m * (1 - m))) ## À vérifier (Madden et al 2007, Turechek et al 2011)... quoi qu'il en soit, je prefere en treme de proportion car p et v sont homogene dans ce cas-là
        }
    )
    attr(D, "flavor") <- flavor
    attr(D, "N") <- N
    attr(D, "n") <- n
    structure(D, class = c("fisher", "index", "numeric"))
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.index <- function(x, ...) cat(x, "\n")

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
fisher.count <- function(x, ...) {
    fisher.default(map_data(x)$r, flavor = "count")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
fisher.incidence <- function(x, ...) {
    mapped_data <- map_data(x)
    n <- unique(mapped_data$n)
    stopifnot(length(n) == 1)
    fisher.default(mapped_data$r, flavor = "incidence", n = n)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
chisq.test <- function(x, ...) UseMethod("chisq.test")

#------------------------------------------------------------------------------#
#' @method chisq.test default
#' @export
#------------------------------------------------------------------------------#
chisq.test.default <- function(x, ...) stats::chisq.test(x, ...)

#------------------------------------------------------------------------------#
#' @method chisq.test fisher
#' @export
#------------------------------------------------------------------------------#
chisq.test.fisher <- function(x, ...) {
    call   <- match.call()
    N      <- attr(x, "N")
    C      <- (N - 1) * x
    df     <- N - 1
    pvalue <- 1 - pchisq(C, df)
    chisq  <- qchisq(pvalue, df)

   # chisq_test <-
    structure(list(
        method = "Pearson's Chi-squared test",
        #estimate = C,### ???????????
        statistic = c("X-squared" = chisq),
        parameter = c("df" = df),
        p.value = pvalue#,
        #data.name = ??,
        #observed = ??,
        #expected = ??,
        #residuals = ??
    ), class = "htest") # A propos "htest": http://www.inside-r.org/node/219017

    # structure(list(
    #     call = call,
    #     index = x,
    #     test = chisq_test
    # ), class = "aggIndex")

}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
z.test <- function(x, ...) UseMethod("z.test")

#------------------------------------------------------------------------------#
#' @method z.test default
#' @export
#------------------------------------------------------------------------------#
z.test.default <- function(x, ...) {
    if (!requireNamespace("snse", quietly = TRUE)) {
        stop("snse needed for this function to work. Please install it.",
             call. = FALSE)
    }
    snse::z.test(x, ...)
}

#------------------------------------------------------------------------------#
#' @source Sun P, Madden LV.
#' @method z.test fisher
#' @export
#------------------------------------------------------------------------------#
z.test.fisher <- function(x, ...) {
    call   <- match.call()
    method <- "One-sample z-test (two-sided)" # Sure????
    N      <- attr(x, "N")
    n      <- attr(x, "n") # n => incidence... to double-check!!!
    m      <- 1
    v      <- 2 * (1 - 1/n) / (N - 1)
    se     <- sqrt(v / n)
    z      <- (x - m) / se
    pvalue <- 2 * pnorm(abs(z), lower.tail = FALSE)
    structure(list(
        method = method,
        #estimate = C,### ???????????
        statistic = c("z" = z),
        #parameter = c("df" = df),
        p.value = pvalue#,
        #data.name = ??,
        #observed = ??,
        #expected = ??,
        #residuals = ??
    ), class = "htest") # A propos "htest": http://www.inside-r.org/node/219017

}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
calpha.test <- function(x, ...) UseMethod("calpha.test")

#------------------------------------------------------------------------------#
#' @source Neyman Tarone
#' @method calpha.test fisher
#' @export
#------------------------------------------------------------------------------#
calpha.test.fisher <- function(x, ...) {
    call   <- match.call()
    method <- "One-sample z-test (one-sided)" # Sure????
    N      <- attr(x, "N")
    n      <- attr(x, "n") # n => incidence... to double-check!!!
    z     <- ((n * (N - 1) * x) - n * N) / sqrt(2 * N * (n^2 - n)) # a one-sided test
    pvalue <- pnorm(abs(z), lower.tail = FALSE) # to check
    structure(list(
        method = method,
        #estimate = C,### ???????????
        statistic = c("z" = z),
        #parameter = c("df" = df),
        p.value = pvalue#,
        #data.name = ??,
        #observed = ??,
        #expected = ??,
        #residuals = ??
    ), class = "htest") # A propos "htest": http://www.inside-r.org/node/219017

}


#------------------------------------------------------------------------------#
#' @describeIn indices Morisita's coefficient of dispersion
#'
#' Interpretation for incidence data remains uncertain, except for low p (Madden and Hughes, 1995)
#'
#' @references
#'
#' Morisita M. 1962. Iδ-Index, a measure of dispersion of
#' individuals. Researches on Population Ecology 4, 1–7.
#' \href{http://dx.doi.org/doi:10.1007/BF02533903}{doi:10.1007/BF02533903}
#'
#' Madden LV, Hughes G. 1995. Plant disease incidence: Distributions,
#' heterogeneity, and temporal analysis. Annual Review of Phytopathology 33 (1):
#' 529–64.
#' \href{http://dx.doi.org/doi:10.1146/annurev.py.33.090195.002525}{doi:10.1146/annurev.py.33.090195.002525}
#'
#' @export
#------------------------------------------------------------------------------#
### To comapre with: dispindmorisita {vegan}, http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/dispindmorisita.html
morisita1 <- function(data) {
    stopifnot(is.count(data) || is.incidence(data))
    call <- match.call()
    d <- data$obs$d # All the raw data
    tot <- sum(d)     # Total number of individuals over all the SU
    N <- length(d)  # Number of sampling units (aka quadrat, aka sample unit)
    I_delta <- N * sum(d * (d - 1)) / (tot * (tot - 1))
    if (is.Count(data)) {
        return(structure(list(call = call, index = I_delta), class = "aggIndex"))
    }
    if (is.Incidence(data)) {# Modified Morisita's index (no more than n possible individuals in a quadrat)
        n <- unique(data@obs$n)
        if (length(n) != 1) stop(paste0("Current implementation only deals ",
                                        "with equal size sampling units."))
        I_B <- I_delta * (n - 1/N) / (n - 1) ### Selon Morisita 1962
        ### On peut écrire I_delta * n / (n - 1) [ou I_delta / (1 - (1/n))] si N grand (Madden et Hughes 1995)
        ### Other remarks from Madden Hughes 1995, Shiyomi Taki and another one to take into account?
        return(structure(list(call = call, index = I_B), class = "aggIndex"))
    }
}

###

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
morisita <- function(x, ...) UseMethod("morisita")

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
morisita.default <- function(x, flavor = c("count", "incidence"), n = NULL, ...) {
    stopifnot(is.numeric(x))
    flavor <- match.arg(flavor)
    #call <- match.call()
    x <- na.omit(x)
    tot <- sum(x)
    N <- length(x)
    Idelta <- N * sum(x * (x - 1)) / (tot * (tot - 1))
    switch (flavor,
        "count" = {
            # Do no more calculation
        },
        "incidence" = {
            stopifnot(!is.null(n))
            stopifnot(n > 1)
            Idelta <- Idelta * (n - 1/N) / (n - 1)
        }
    )
    structure(Idelta, class = c("morisita", "index", "numeric"))
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
morisita.count <- function(x, ...) {
    morisita.default(map_data(x)$r, flavor = "count")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
morisita.incidence <- function(x, ...) {
    mapped_data <- map_data(x)
    n <- unique(mapped_data$n)
    stopifnot(length(n) == 1)
    morisita.default(mapped_data$r, flavor = "incidence", n = n)
}

#------------------------------------------------------------------------------#
#' @describeIn indices Lloyd's index of patchiness or Lloyd's mean crowding
#'
#'   Returns a value between 0 and \code{n}, the size of a sampling unit. A
#'   value lower than 1 indicates a random pattern, while a value greater than 1
#'   indicates an aggregated pattern. The degree of aggregation increases with
#'   the value of the index.
#'
#'   LIP < 1: means REGULAR
#'   LIP = 1: means RANDOM
#'   LIP > 1: means AGGREGATED
#'
#' @references Lloyd M. 1967. Mean crowding. The Journal of Animal Ecology 36,
#'   1–30.
#'
#' @export
#------------------------------------------------------------------------------#
lloyd <- function(x, type) UseMethod("lloyd")

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
lloyd.default <- function(x, type = c("patchiness", "mean-crowding")) {
    stopifnot(is.numeric(x))
    type <- match.arg(type)
    #call <- match.call()
    x <- na.omit(x)
    m <- mean(x)
    stopifnot(m > 1)
    v <- var(x)
    mc <- m + ((v / m) - 1)
    switch (type,
        "mean-crowding" = return(mc),
        "patchiness"    = {
            structure(mc / m, class = c("lloyd", "index", "numeric"))
        }
    )
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
lloyd.count <- function(x, type = c("patchiness", "mean-crowding")) {
    lloyd.default(map_data(x)$r, type)
}

######### BELOW REMOVE

#------------------------------------------------------------------------------#
#' @describeIn indices Iwao's patchiness regression index
#'
#' @references
#' Iwao S. 1968. A new regression method for analyzing the aggregation pattern
#' of animal populations. Researches on Population Ecology 10, 1–20.
#' \href{http://dx.doi.org/10.1007/BF02514729}{doi:10.1007/BF02514729}
#'
#' @export
#------------------------------------------------------------------------------#
iwao <- function(data) {
    stopifnot(is.count(data))
    if (is.numeric(data)) {
        if (!all(is.wholenumber(data))) stop("no whole numbers!")
    }
    return(NULL)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.aggIndex <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\nIndex: ", x$index, "\n\n", sep ="")

    # Below to transfert to summary after // et améliorer le rendu
    if ("test" %in% names(x)) print(x$test)
}

#Campbell, C.L., Madden, L.V., 1990. Spatial aspects of plant disease epidemics: Analysis of spatial pattern. Introduction to Plant Disease Epidemiology 289–328.
#### À remplacer par Madden et al 2007?????

#Madden, L.V., Hughes, G., 1995. Plant disease incidence: distributions, heterogeneity, and temporal analysis. Annual Review of Phytopathology 33, 529–564. doi:10.1146/annurev.py.33.090195.002525


#------------------------------------------------------------------------------#
#' Patchiness/Aggregation/Clustering indices
#'
#' TO DO
#'
#' Different indices are available, and depend on the nature of the data set.
#' Only indices for \code{Count} and \code{Incidence} data are currently
#' implemented.
#'
#' @param data The data set you want to explore.
#' @param method Method implemented.
#'
#' @export
#------------------------------------------------------------------------------#
aggIndices <- function(data, method = "all") {
    stopifnot(is.count(data) || is.incidence(data))
    funs <- list(count = list(
                    fisher = fisher,
                    morisita = morisita,
                    lloyd = lloyd
                ),
                 incidence = list(
                     fisher = fisher,
                     morisita = morisita
                ))
    res <- lapply(funs[class(data) == names(funs)][[1]], function(f) f(data))
    structure(res, class = "aggIndices", call = match.call())
}

#obj$call  <- call
#obj$name  <- c(full = "Fisher's index of dispersion (D)", short = "Fisher's")
#obj$index <- index
#obj$tests <- list(cAlpha, z, chiSq)

#aggIndexes$fisher$index

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.aggIndices <- function(x, ...) {
    cat("\nCall:\n")
    print(attr(x, "call"))
    cat("\nIndices:\n")
    for (idx in names(x)) {
        cat(" - ", idx, ": ", x[[idx]]$index, sep = "")
        if ("test" %in% names(x[[idx]]))
            cat(" (p-value: ", x[[idx]]$test$p.value, ")", sep = "")
        cat("\n")
    }
}











