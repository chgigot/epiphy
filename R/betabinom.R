#------------------------------------------------------------------------------#
#' The beta-binomial distribution.
#'
#' Density, distribution function, quantile function and random generation for
#' the beta-binomial distribution with parameters \code{size}, \code{prob},
#' \code{theta}, \code{shape1}, \code{shape2}. This distribution corresponds to
#' an overdispersed binomial distribution.
#'
#' Be aware that in this implementation \code{theta} = 1 / (\code{shape1} +
#' \code{shape2}). \code{prob} and \code{theta}, or \code{shape1} and
#' \code{shape2} must be specified. if \code{theta} = 0, use *binom family
#' instead.
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param size Number of trials.
#' @param prob Probability of success on each trial.
#' @param theta Aggregation parameter (theta = 1 / (shape1 + shape2)).
#' @param shape1,shape2 Shape parameters.
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are
#'     \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
#'
#' @return
#'
#' \code{dbetabinom} gives the density, \code{pbetabinom} gives the distribution
#' function, \code{qbetabinom} gives the quantile function and \code{rbetabinom}
#' generates random deviates.
#'
#' @examples
#' # Compute P(25 < X < 50) for X following the Beta-Binomial distribution
#' # with parameters size = 100, prob = 0.5 and theta = 0.35:
#' sum(dbetabinom(25:50, size = 100, prob = 0.5, theta = 0.35))
#'
#' # When theta tends to 0, dbetabinom outputs tends to dbinom outputs:
#' sum(dbetabinom(25:50, size = 100, prob = 0.5, theta = 1e-7))
#' sum(dbetabinom(25:50, size = 100, shape1 = 1e7, shape2 = 1e7))
#' sum(dbinom(25:50, size = 100, prob = 0.5))
#'
#' # Example of binomial and beta-binomial frequency distributions:
#' n   <- 15
#' q   <- 0:n
#' p1  <- dbinom(q, size = n, prob = 0.33)
#' p2  <- dbetabinom(q, size = n, prob = 0.33, theta = 0.22)
#' res <- rbind(p1, p2)
#' dimnames(res) <- list(c("Binomial", "Beta-binomial"), q)
#' barplot(res, beside = TRUE, legend.text = TRUE, ylab = "Frequency")
#'
#' # Effect of the aggregation parameter theta on probability density:
#' thetas <- seq(0.001, 2.5, by = 0.001)
#' density1 <- rep(sum(dbinom(25:50, size = 100, prob = 0.5)), length(thetas))
#' density2 <- sapply(thetas, function(theta) {
#'     sum(dbetabinom(25:50, size = 100, prob = 0.5, theta = theta))
#' })
#' plot(thetas, density2, type = "l",
#'      xlab = expression("Aggregation parameter ("*theta*")"),
#'      ylab = "Probability density between 25 and 50 (size = 100)")
#' lines(thetas, density1, lty = 2)
#'
#' @seealso \code{\link[emdbook]{dbetabinom}} in the package \strong{emdbook}
#' where the definition of theta is different.
#'
#' @name BetaBinomial
#' @export
#------------------------------------------------------------------------------#
dbetabinom <- function(x, size, prob, theta, shape1, shape2, log = FALSE) {
    list2env(check_betabinom(prob, theta, shape1, shape2), envir = environment())
    lpmf <- rep(-Inf, length(x))
    if (!all(int <- is.wholenumber(x))) {
        call <- match.call()
        warning(paste0("In ", deparse(call), " : non-integer x = ",
                       eval(call$x)[!int], "\n"), call. = FALSE)
    }
    lpmf[int] <- lchoose(size, x[int]) +
        lbeta(x[int] + shape1, size - x[int] + shape2) -
        lbeta(shape1, shape2)
    if (log) lpmf else exp(lpmf)
}

#------------------------------------------------------------------------------#
#' @rdname BetaBinomial
#' @export
#------------------------------------------------------------------------------#
pbetabinom <- function(q, size, prob, theta, shape1, shape2, lower.tail = TRUE,
                       log.p = FALSE) {
    # TODO: Double-check.
    list2env(check_betabinom(prob, theta, shape1, shape2), envir = environment())
    q <- trunc(q)
    p <- sapply(q, function(x)
        sum(vapply(0:x, dbetabinom, FUN.VALUE = numeric(1L), size = size,
                   shape1 = shape1, shape2 = shape2)))
    if (!lower.tail) p <- 1 - p
    if (log.p) log(p) else p
}

#------------------------------------------------------------------------------#
#' @rdname BetaBinomial
#' @export
#------------------------------------------------------------------------------#
qbetabinom <- function(p, size, prob, theta, shape1, shape2, lower.tail = TRUE,
                       log.p = FALSE) {
    call <- match.call()
    list2env(check_betabinom(prob, theta, shape1, shape2), envir = environment())
    if (log.p) p <- exp(p)
    q <- sapply(p, function(x) {
        if (0 > x | x > 1) return(NaN)
        y <- 0
        while (x > pbetabinom(y, size = size, shape1 = shape1, shape2 = shape2))
            y <- y + 1
        y
    })
    if (lower.tail) q else size - q
}

#------------------------------------------------------------------------------#
#' @rdname BetaBinomial
#' @export
#------------------------------------------------------------------------------#
rbetabinom <- function(n, size, prob, theta, shape1, shape2) {
    list2env(check_betabinom(prob, theta, shape1, shape2), envir = environment())
    rbinom(n, size = size, prob = rbeta(n, shape1 = shape1, shape2 = shape2))
}

#------------------------------------------------------------------------------#
#' @keywords internal
#------------------------------------------------------------------------------#
check_betabinom <- function(prob, theta, shape1, shape2) {#, env) {
    pair1 <- pair2 <- FALSE
    if (!missing(prob) && !missing(theta)) pair1 <- TRUE
    if (!missing(shape1) && !missing(shape2)) pair2 <- TRUE
    if (!xor(pair1, pair2)) {
        stop("(prob, theta) or (shape1, shape2) must be specified.")
    }
    if (pair1) {
        shape1 <- prob / theta
        shape2 <- (1 - prob) / theta
    }
    list(shape1 = shape1, shape2 = shape2)
}



