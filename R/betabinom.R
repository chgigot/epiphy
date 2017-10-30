#------------------------------------------------------------------------------#
#' The Beta-Binomial Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the beta-binomial distribution with parameters \code{size}, \code{prob},
#' \code{theta}, \code{shape1}, \code{shape2}. It corresponds to as an
#' overdispersed binomial distribution.
#'
#' Be careful, in this implementation, \code{theta} = 1 / (\code{shape1} + \code{shape2}).
#' \code{prob} and \code{theta}, or \code{shape1} and \code{shape2} must be specified.
#' if \code{theta} = 0, use *binom family instead.
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param size Number of trials.
#' @param prob Probability of success on each trial.
#' @param theta Aggregation parameter (theta = 1 / (shape1 + shape2)).
#' @param shape1,shape2 Shape parameters.
#' @param log,log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail [Not yet implemented].
#'
#' @return
#'
#' \code{dbetabinom} gives the density, \code{pbetabinom} gives the distribution
#' function, \code{qbetabinom} gives the quantile function and \code{rbetabinom}
#' generates random deviates.
#'
#' @examples
#' # Compute P(25 < X < 50) for X Beta-Binomial(100, 0.5, 0.35)
#' sum(dbetabinom(25:50, size = 100, prob = 0.5, theta = 0.35))
#'
#' # When theta -> 0, dbetabinom -> dbinom
#' sum(dbetabinom(25:50, size = 100, prob = 0.5, theta = 1e-7))
#' sum(dbetabinom(25:50, size = 100, shape1 = 1e7, shape2 = 1e7))
#' sum(dbinom(25:50, size = 100, prob = 0.5))
#'
#' n <- 15
#' q <- 0:n
#' p1 <- dbinom(q, size = n, prob = 0.33)
#' p2 <- dbetabinom(q, size = n, prob = 0.33, theta = 0.22)
#' res <- rbind(p1, p2)
#' dimnames(res) <- list(c("Binomial", "Beta-binomial"), q)
#' barplot(res, beside=TRUE, legend.text = TRUE, ylab = "Frequency")
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

#x <- sapply(seq(0.01,10,by = 0.01), function(x) 1/10^x)
#y1<-sapply(x, function(x) sum(dbinom(25:50, size = 100, prob = 0.5)))
#y2<-sapply(x, function(x) sum(dbetabinom(25:50, size = 100, prob = 0.5, theta = x)))
#plot(x,y2, type = "l");lines(x,y1,lty=2)

### Compute expected beta-binomial probabilities using approach in book from
# Madden et al, 2007
# pbbinom <- function( q, size, prob, theta ) { # q, size belongs to Int
#     for (x in q) {
#         alpha <- prob / theta
#         beta  <- (1-prob)/theta
#         y     <- 1
#         for (i in 0:(size-1))
#             y <- y * ((1-prob+i*theta) / (1+i*theta))
#         if (x!=0) {
#             k <- y
#             for (j in 1:x) {
#                 k <- k * ((size+1-j)/j) * ((alpha-1+j)/(size+beta-j))
#                 y <- y + k
#             }
#         }
#         if ( !exists("val", inherits = FALSE) ) {val <- y} else {val <- c(val, y)} # SI NON EXIST VAL DANS L'ENVIRONNEMENT DE LA FONCTION !!! TRES IMPORTANT
#     }
#     return(val)
# }
#
# dbetabinom2 <- function(x, size, prob, theta) {
#     values <- pbbinom(x, size, prob, theta)
#     (values - c(0, values[1:max(x)])) # Pas propre // temporaire
# }

#------------------------------------------------------------------------------#
#' @rdname BetaBinomial
#' @export
#------------------------------------------------------------------------------#
pbetabinom <- function(q, size, prob, theta, shape1, shape2, lower.tail = TRUE,
                       log.p = FALSE) {
    list2env(check_betabinom(prob, theta, shape1, shape2), envir = environment())
    q <- trunc(q)
    p <- sapply(q, function(x)
        sum(vapply(0:x, dbetabinom, FUN.VALUE = numeric(1L), size = size,
                   shape1 = shape1, shape2 = shape2)))
    if (log.p) log(p) else p
} # TODO: To double check
# TODO: Gerer le lower.tail

#------------------------------------------------------------------------------#
#' @rdname BetaBinomial
#' @export
#------------------------------------------------------------------------------#
qbetabinom <- function(p, size, prob, theta, shape1, shape2, lower.tail = TRUE,
                       log.p = FALSE) {
    list2env(check_betabinom(prob, theta, shape1, shape2), envir = environment())
    if (log.p) p <- 10^p ### TODO: Non non non, c'est exponentiel
    q <- sapply(p, function(x) {
        if (0 > x | x > 1) return(NaN)
        y <- 0
        while (x > pbetabinom(y, size = size, shape1 = shape1, shape2 = shape2))
            y <- y + 1
        y
    })
    q
    # TODO: Gerer le lower.tail
    # Comparer avec le code C de qbinom
    #qbinom(p = 2, size = 10, prob = 0.5)
    #[1] NaN
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
    if (!xor(pair1, pair2)) stop("(prob, theta) or (shape1, shape2) must be specified.")
    if (pair1) {
        shape1 <- prob / theta
        shape2 <- (1 - prob) / theta
    }
    list(shape1 = shape1, shape2 = shape2)
}



