#------------------------------------------------------------------------------#
#' @include utils.R
#' @include intensity-classes.R
#' @include betabinom.R
#' @include mle-factory.R
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' Maximum likelihood fitting of two distributions and goodness-of-fit
#' comparison.
#'
#' Different distributions may be used depending on the kind of provided data.
#' By default, the Poisson and negative binomial distributions are fitted to
#' count data, whereas the binomial and beta-binomial distributions are used
#' with incidence data. Either Randomness assumption (Poisson or binomial
#' distributions) or aggregation assumption (negative binomial or beta-binomial)
#' are made, and then, a goodness-of-fit comparison of both distributions is
#' made using a log-likelihood ratio test.
#'
#' Under the hood, \code{distr_fit} relies on the \code{\link{smle}} utility
#' which is a wrapped around the \code{\link[stats]{optim}} procedure.
#'
#' Note that there may appear warnings about chi-squared goodness-of-fit tests
#' if any expected count is less than 5 (Cochran's rule of thumb).
#'
#' @param data An \code{intensity} object.
#' @param random Distribution to describe random patterns.
#' @param aggregated Distribution to describe aggregated patterns.
#' @param n_est Number of estimated parameters for both distributions.
#' @param ... Additional arguments to be passed to other methods.
#'
#' @return
#'
#' An object of class \code{fit_two_distr}, which is a list containing at least
#' the following components:
#' \tabular{ll}{
#'     \code{call}  \tab The function \code{\link[base]{call}}. \cr
#'     \code{name}  \tab The names of both distributions. \cr
#'     \code{model} \tab The outputs of fitting process for both distributions. \cr
#'     \code{llr}   \tab The result of the log-likelihood ratio test. \cr
#' }
#' Other components can be present such as:
#' \tabular{ll}{
#'     \code{param} \tab A numeric matrix of estimated parameters (that can be
#'                       printed using \code{\link[stats]{printCoefmat}}). \cr
#'     \code{freq}  \tab A data frame or a matrix with the observed and expected
#'                       frequencies for both distributions for the different
#'                       categories. \cr
#'     \code{gof}   \tab Goodness-of-fit tests for both distributions (which are
#'                       typically chi-squared goodness-of-fit tests). \cr
#' }
#'
#' @examples
#' # Simple workflow for incidence data:
#' my_data <- count(arthropods)
#' my_data <- split(my_data, by = "t")[[3]]
#' my_res  <- fit_two_distr(my_data)
#' summary(my_res)
#' plot(my_res)
#'
#' # Simple workflow for incidence data:
#' my_data <- incidence(tobacco_viruses)
#' my_res  <- fit_two_distr(my_data)
#' summary(my_res)
#' plot(my_res)
#'
#' # Note that there are other methods to fit some common distributions.
#' # For example for the Poisson distribution, one can use glm:
#' my_arthropods <- arthropods[arthropods$t == 3, ]
#' my_model <- glm(my_arthropods$i ~ 1, family = poisson)
#' lambda <- exp(coef(my_model)[[1]]) # unique(my_model$fitted.values) works also.
#' lambda
#' # ... or the fitdistr function in MASS package:
#' require(MASS)
#' fitdistr(my_arthropods$i, "poisson")
#'
#' # For the binomial distribution, glm still works:
#' my_model <- with(tobacco_viruses, glm(i/n ~ 1, family = binomial, weights = n))
#' prob <- logit(coef(my_model)[[1]], rev = TRUE)
#' prob
#' # ... but the binomial distribution is not yet recognized by MASS::fitdistr.
#'
#' # Examples featured in Madden et al. (2007).
#' # p. 242-243
#' my_data <- incidence(dogwood_anthracnose)
#' my_data <- split(my_data, by = "t")
#' my_fit_two_distr <- lapply(my_data, fit_two_distr)
#' lapply(my_fit_two_distr, function(x) x$param$aggregated[c("prob", "theta"), ])
#' lapply(my_fit_two_distr, plot)
#'
#' my_agg_index <- lapply(my_data, agg_index)
#' lapply(my_agg_index, function(x) x$index)
#' lapply(my_agg_index, chisq.test)
#'
#' @references
#'
#' Madden LV, Hughes G. 1995. Plant disease incidence: Distributions,
#' heterogeneity, and temporal analysis. Annual Review of Phytopathology 33(1):
#' 529–564.
#' \href{http://dx.doi.org/doi:10.1146/annurev.py.33.090195.002525}{doi:10.1146/annurev.py.33.090195.002525}
#'
#' @export
#------------------------------------------------------------------------------#
fit_two_distr <- function(data, ...) UseMethod("fit_two_distr")

#------------------------------------------------------------------------------#
#' @rdname fit_two_distr
#' @export
#------------------------------------------------------------------------------#
fit_two_distr.default <- function(data, random, aggregated, ...) {
    call            <- match.call()
    name_random     <- attr(random, "name")
    name_aggregated <- attr(aggregated, "name")
    # Fit two distributions to the data and perform a log-likelihood ratio test:
    random     <- random(data)
    aggregated <- aggregated(data)
    llr        <- NULL
    if (!is.null(random$full_output) && !is.null(aggregated$full_output)) {
        llr <- llrtest(random, aggregated)
    }
    structure(list(call  = call,
                   name  = list(random = name_random, aggregated = name_aggregated),
                   model = list(random = random, aggregated = aggregated),
                   llr   = llr),
              class = "fit_two_distr")
}

#------------------------------------------------------------------------------#
#' @rdname fit_two_distr
#' @export
#------------------------------------------------------------------------------#
fit_two_distr.count <- function(data, random = smle_pois,
                                aggregated = smle_nbinom,
                                n_est = c(random = 1, aggregated = 2), ...) {
    call <- match.call()
    mapped_data <- map_data(data)
    ind <- mapped_data[["i"]]
    N <- nrow(mapped_data)
    #--------------------------------------------------------------------------#
    # Try to fit two distributions:
    res <- fit_two_distr.default(data, random = random, aggregated = aggregated,
                             ...)
    # Tweak and complete the result:
    res$call <- call
    res$data_class <- class(data)
    #--------------------------------------------------------------------------#
    # Retrieve and estimate all the parameters:
    param <- list(random = NULL, aggregated = NULL)
    if (!is.null(res$model$random$full_output)) {
        param$random <- coef(summary(res$model$random))
        param$random <- param$random[sort(rownames(param$random)),, drop = FALSE]
    }
    if (!is.null(res$model$aggregated$full_output)) {
        param$aggregated <- coef(summary(res$model$aggregated))
        new_param <- list(
            prob = estimate_param(res$model$aggregated, quote(x2/(x1 + x2)), "norm")
        )
        param$aggregated <- rbind_param(param$aggregated, new_param)
        param$aggregated <- param$aggregated[sort(rownames(param$aggregated)),, drop = FALSE]
    }
    #--------------------------------------------------------------------------#
    # Summary data frame (with observed and theoretical frequencies):
    range_ind <- 0:max(ind)
    freq <- setNames(as.data.frame(table(ind)), c("category", "observed"))
    # To transform a factor into its original numeric values:
    freq$category <- as.numeric(levels(freq$category))[freq$category]
    freq <- merge(x = data.frame(category = range_ind), y = freq,
                  by = "category", all = TRUE)
    freq[is.na(freq)] <- 0
    freq$random <- NA
    if (!is.null(res$model$random$full_output)) {
        freq$random <- N * dpois(x = range_ind,
                                 lambda = coef(res$model$random)[["lambda"]])
    }
    freq$aggregated <- NA
    if (!is.null(res$model$aggregated$full_output)) {
        freq$aggregated <- N * dnbinom(x = range_ind,
                                 mu = coef(res$model$aggregated)[["mu"]],
                                 size = coef(res$model$aggregated)[["k"]])
    }
    #--------------------------------------------------------------------------#
    gof <- list(random = NULL, aggregated = NULL)
    if (!is.null(res$model$random$full_output)) {
        gof$random <- chisq.test2(freq$observed, freq$random,
                                  n_est = n_est[["random"]], rescale.p = TRUE)
    }
    if (!is.null(res$model$aggregated$full_output)) {
        gof$aggregated <- chisq.test2(freq$observed, freq$aggregated,
                                      n_est = n_est[["aggregated"]], rescale.p = TRUE)
    }
    #--------------------------------------------------------------------------#
    res$param <- param
    res$freq  <- freq
    res$gof   <- gof
    res
}

#------------------------------------------------------------------------------#
#' @rdname fit_two_distr
#' @export
#------------------------------------------------------------------------------#
fit_two_distr.incidence <- function(data, random = smle_binom,
                                    aggregated = smle_betabinom,
                                    n_est = c(random = 1, aggregated = 2), ...) {
    call <- match.call()
    mapped_data <- map_data(data)
    ind <- mapped_data[["i"]]
    N   <- nrow(mapped_data)
    n   <- unique(mapped_data[["n"]])
    if (length(n) != 1) stop(paste0("Current implementation only deals ",
                                    "with equal size sampling units."))
    #--------------------------------------------------------------------------#
    # Try to fit two distributions:
    res <- fit_two_distr.default(mapped_data, random = random,
                             aggregated = aggregated, ...)
    # Tweak and complete the result:
    res$call <- call
    res$data_class <- class(data)
    #--------------------------------------------------------------------------#
    # Retrieve and estimate all the parameters:
    param <- list(random = NULL, aggregated = NULL)
    if (!is.null(res$model$random$full_output)) {
        param$random <- coef(summary(res$model$random))
        param$random <- param$random[sort(rownames(param$random)),, drop = FALSE]
    }
    if (!is.null(res$model$aggregated$full_output)) {
        param$aggregated <- coef(summary(res$model$aggregated))
        new_param <- list(
            rho   = estimate_param(res$model$aggregated, quote(x2/(x2 + 1)), "norm"),
            alpha = estimate_param(res$model$aggregated, quote(x1/x2), "norm"),
            beta  = estimate_param(res$model$aggregated, quote((1 - x1)/x2), "norm")
        )
        param$aggregated <- rbind_param(param$aggregated, new_param)
        param$aggregated <- param$aggregated[sort(rownames(param$aggregated)),, drop = FALSE]
    }
    #--------------------------------------------------------------------------#
    # Summary data frame (with observed and theoretical frequencies)
    range_ind <- 0:n
    freq <- setNames(as.data.frame(table(ind)), c("category", "observed"))
    # To transform a factor into its original numeric values:
    freq$category <- as.numeric(levels(freq$category))[freq$category]
    freq <- merge(x = data.frame(category = range_ind), y = freq,
                  by = "category", all = TRUE)
    freq[is.na(freq)] <- 0
    if (!is.null(res$model$random$full_output)) {
        freq$random <- N * dbinom(x = range_ind, size = n,
                                  prob = coef(res$model$random)[["prob"]])
    }
    if (!is.null(res$model$aggregated$full_output)) {
        freq$aggregated <- N * dbetabinom(x = range_ind, size = n,
                                  prob  = coef(res$model$aggregated)[["prob"]],
                                  theta = coef(res$model$aggregated)[["theta"]])
    }
    #--------------------------------------------------------------------------#
    gof <- list(random = NULL, aggregated = NULL)
    if (!is.null(res$model$random$full_output)) {
        gof$random <- chisq.test2(freq$observed, freq$random,
                                  n_est = n_est[["random"]], rescale.p = TRUE)
    }
    if (!is.null(res$model$aggregated$full_output)) {
        gof$aggregated <- chisq.test2(freq$observed, freq$aggregated,
                                      n_est = n_est[["aggregated"]], rescale.p = TRUE)
    }
    #--------------------------------------------------------------------------#
    res$param <- param
    res$freq  <- freq
    res$gof   <- gof
    res
}

#==============================================================================#
# Print, summary and plot
#==============================================================================#

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.fit_two_distr <- function(x, ...) {
    cat("Fitting of two distributions by maximum likelihood\n",
        "for '", x$data_class[1L], "' data.\n",
        "\nParameter estimates:\n",
        "\n(1) ", x$name$random, " (random):\n", sep = "")
    if (!is.null(x$param$random)) {
        print(x$param$random[, 1:2, drop = FALSE])
    } else {
        cat("Maximum likelihood procedure did not succeed.\n")
    }
    cat("\n(2) ", x$name$aggregated, " (aggregated):\n", sep = "")
    if (!is.null(x$param$aggregated)) {
        print(x$param$aggregated[, 1:2, drop = FALSE])
    } else {
        cat("Maximum likelihood procedure did not succeed.\n")
    }
    invisible(x)
}

#------------------------------------------------------------------------------
#' @export
#------------------------------------------------------------------------------
summary.fit_two_distr <- function(object, ...) {
    # TODO: to improve, add a coef function, ...
    res <- object
    structure(res, class = "summary.fit_two_distr")
}

#------------------------------------------------------------------------------
#' @method print summary.fit_two_distr
#' @export
#------------------------------------------------------------------------------
print.summary.fit_two_distr <- function(x, ...) {
    cat("Fitting of two distributions by maximum likelihood\n",
        "for '", x$data_class[1L], "' data.\n",
        "Parameter estimates:\n",
        "\n(1) ", x$name$random, " (random):\n", sep = "")
    if (!is.null(x$param$random)) {
        printCoefmat(x$param$random)
    } else {
        cat("Maximum likelihood procedure did not succeed.\n")
    }
    cat("\n(2) ", x$name$aggregated, " (aggregated):\n", sep = "")
    if (!is.null(x$param$aggregated)) {
        printCoefmat(x$param$aggregated)
    } else {
        cat("Maximum likelihood procedure did not succeed.\n")
    }
    invisible(x)
}

#------------------------------------------------------------------------------
#' @export
#------------------------------------------------------------------------------
plot.fit_two_distr <- function(x, ..., breaks = NULL) {
    if (!is.null(breaks)) {
        stopifnot(is.numeric(breaks))
        stopifnot(is.wholenumber(breaks))
        stopifnot(breaks > 0)
        fac <- rep(0:ceiling(nrow(x$freq)/breaks), each = breaks)[1:nrow(x$freq)]
        tmp <- split(x$freq, fac, drop = TRUE)
        x$freq <- do.call(rbind, lapply(tmp, function(xx) {
            res <- xx[0, ] # Trick to create an "empty" relevant df
            res[1,1] <- paste0(xx[1, 1], "-", xx[nrow(xx), 1]) # changement auto en char pas de probleme
            res[1, 2:ncol(res)] <- colSums(xx[2:ncol(xx)])
            res
        }))
    }
    rownames(x$freq) <- x$freq$category
    freq_long <- setNames(as.data.frame.table(as.matrix(x$freq[-1])),
                          c("Number per sampling unit", "key", "Frequency"))
    freq_long$key <- factor(tocamel(freq_long$key),
                            levels = c("Observed", "Aggregated", "Random"))
    gg <- ggplot(freq_long, aes(x = `Number per sampling unit`, y = Frequency,
                                fill = key))
    gg <- gg + geom_bar(stat = "identity", position = "dodge", color = "black",
                        size = 0.5)
    gg <- gg + scale_fill_manual(values = c(Observed = "black",
                                            Aggregated = "grey",
                                            Random = "white"))
    gg <- gg + theme_bw()
    gg <- gg + theme(legend.justification = c(0.98, 0.98),
                     legend.position = c(0.98, 0.98),
                     legend.background = element_rect(size = 0.5,
                                                      linetype = "solid",
                                                      color = "black"),
                     legend.title = element_blank())
    print(gg)
    invisible(NULL)
}




