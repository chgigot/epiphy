#------------------------------------------------------------------------------#
#' @include utils.R
#' @include intensity-classes.R
#' @include betabinom.R
#' @include mle-factory.R
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' Try to fit distributions to data sets
#'
#' Depending of the nature of the data submitted, different distributions are
#' tryied to be fitted to the data.
#'
#' Two distributions are used. One under the assumption of no aggregation and
#' the other one under the assumption of aggregation.
#'
#' @param x A data of the class \code{Intensity}.
#' @param type The kind of distribution to fit: "random", "aggregated" or "both".
#'
#' @return An object of class \code{DistFitting}, with the following components:
#'
#' \tabular{ll}{
#'     \code{from} \tab Kind of the \code{Intensity} object. \cr
#'     \code{random} \tab Theoretical data under random assumption. \cr
#'     \code{aggregated} \tab Theoretical data under aggregated assumption. \cr
#'     \code{freqTable} \tab Entry 4 \cr
#' }
#'
#' @references
#'
#' @examples
#' \dontrun{
#' data <- Incidence(Cochran1936)
#' fit_distr(data)
#' }
#'
#' # With count data:
#' model   <- glm(aphids$r ~ 1, family = poisson)
#' lambda1 <- exp(coef(model)[[1]])
#' # Or:
#' lambda1 <- unique(model$fitted.values)
#'
#' my_data <- count(aphids)
#' my_fit_distr <- fit_distr(my_data)
#' lambda2 <- my_fit_distr$random$par[[1]]
#'
#' identical(lambda1, lambda2)
#' lambda1 - lambda2
#'
#' require(MASS)
#' fitdistr(aphids$r, "Poisson")
#'
#' #--------
#' # Basic example
#' my_data <- incidence(tomato_t...$field_1929)
#' my_res  <- fit_distr(my_data)
#' my_res
#' summary(my_res)
#' plot(my_res)
#'
#' # Advanced example with another proposed distribution for
#' # the aggregated situation:
#' # L-N-B : logistic-normal-binomial
#' # CLL-N-B: CLL-normal-binomial (Same as LNB, but with CLL link)
#' ll_lnb  <- function(data, param)
#' mle_lnb <- function(data) mle(data, ll_lnb, c(prob = 0.5, sigma = 0.1))
#' my_res2 <- fit_distr(my_data, aggregated = mle_lnb)
#' my_res2
#' summary(res2)
#' plot(res2)
#'
#' #==========
#'
#' data <- tomato_tswv$field_1929[tomato_tswv$field_1929$t == 1, ]
#' data <- clump(incidence(data), unit_size = c(x = 3, y = 3))
#' res <- fit_distr(data)
#' res
#' summary(res)
#' plot(res)
#'
#' @export
#------------------------------------------------------------------------------#
fit_distr <- function(data, ...) UseMethod("fit_distr")

#------------------------------------------------------------------------------#
# fit_distr.default
#------------------------------------------------------------------------------#
fit_distr.default <- function(data, ..., random, aggregated) {
    call        <- match.call()

    # Fit two distributions to the data and perform a log-likelihood ratio test:
    random      <- random(data)
    aggregated  <- aggregated(data)
    test        <- llrtest(random, aggregated)

    structure(list(call = call,
                   model = list(random = random,
                                aggregated = aggregated),
                   test = test),
              class = "fit_distr")
}

#------------------------------------------------------------------------------#
#' @rdname fit_distr
#' @export
#------------------------------------------------------------------------------#
fit_distr.count <- function(data, ..., random = smle_pois,
                            aggregated = smle_nbinom,
                            extra_param, # quotes TODO
                            freq_model) { # functions TOD0
    call <- match.call()
    mapped_data <- map_data(data)
    r <- mapped_data[["r"]]
    N <- nrow(mapped_data)

    # Try to fit two distributions:
    res <- fit_distr.default(data, ..., random = random,
                             aggregated = aggregated)

    # Tweak and complete the result:
    res$call <- call

    # All the parameters:
    param <- list(random     = coef(summary(res$model$random)),
                  aggregated = coef(summary(res$model$aggregated)))
    new_param <- list(
        prob   = estimate_param(res$model$aggregated, quote(x2/(x1 + x2)), type = "norm")
    )
    param$aggregated <- rbind(param$aggregated, do.call(rbind, lapply(new_param, unlist))) # TODO: Not clean

    # Summary data frame (with observed and theoretical frequencies)
    freq <- setNames(as.data.frame(table(r)), c("category", "observed"))
    freq$category <- as.numeric(levels(freq$category))[freq$category] # To transform a factor to its original numeric values:
    freq <- merge(x = data.frame(category = 0:max(r)), y = freq,
                  by = "category", all = TRUE)
    freq[is.na(freq)] <- 0
    freq$random <- N * dpois(x = 0:max(r),
                             lambda = coef(res$model$random)[["lambda"]])
    freq$aggregated <- N * dnbinom(x = 0:max(r) ,
                                   mu = coef(res$model$aggregated)[["mu"]],
                                   size = coef(res$model$aggregated)[["k"]])

    res$param <- param
    res$freq <- freq
    res
}

#------------------------------------------------------------------------------#
#' @rdname fit_distr
#' @export
#------------------------------------------------------------------------------#
fit_distr.incidence <- function(data, ..., random = smle_binom,
                                aggregated = smle_betabinom,
                                extra_param, # quotes TODO
                                freq_model) { # functions TOD0
    call <- match.call()
    mapped_data <- map_data(data)
    r <- mapped_data[["r"]]
    N <- nrow(mapped_data)
    n <- unique(mapped_data[["n"]])
    if (length(n) != 1) stop(paste0("Current implementation only deals ",
                                    "with equal size sampling units."))

    # Try to fit two distributions:
    res <- fit_distr.default(mapped_data, ..., random = random,
                             aggregated = aggregated)

    # Tweak and complete the result:
    res$call <- call

    # All the parameters:
    param <- list(random     = coef(summary(res$model$random)),
                  aggregated = coef(summary(res$model$aggregated)))
    new_param <- list(
        rho   = estimate_param(res$model$aggregated, quote(x2/(x2 + 1)), type = "norm"),
        alpha = estimate_param(res$model$aggregated, quote(x1/x2), type = "norm"),
        beta  = estimate_param(res$model$aggregated, quote((1 - x1)/x2), type = "norm")
    )
    param$aggregated <- rbind(param$aggregated, do.call(rbind, lapply(new_param, unlist))) # TODO: Not clean

    # Summary data frame (with observed and theoretical frequencies)
    freq <- setNames(as.data.frame(table(r)), c("category", "observed"))
    freq$category <- as.numeric(levels(freq$category))[freq$category] # To transform a factor to its original numeric values:
    freq <- merge(x = data.frame(category = 0:n), y = freq,
                  by = "category", all = TRUE)
    freq[is.na(freq)] <- 0
    freq$random <- N * dbinom(x = 0:n, size = n,
                              prob = coef(res$model$random)[["prob"]])
    freq$aggregated <- N * dbetabinom(x = 0:n, size = n,
                              prob  = coef(res$model$aggregated)[["prob"]],
                              theta = coef(res$model$aggregated)[["theta"]])

    res$param <- param
    res$freq <- freq
    res
    # Binomial distribution
    #d <- data$obs$d # All the raw data
    #N <- length(d)
    #n <- unique(data$obs$n)
    #rand$model <- glm(d/n ~ 1, family = binomial, weights = rep(n, N))
}

#==============================================================================#
# Print, summary and plot
#==============================================================================#

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.fit_distr <- function(x, ...) {
    cat_snse <- function(name, x) {
        cat("- ", name, " = ", formatC(x[[1]]), " (\u00b1 ", formatC(x[[2]]),
            ")\n", sep = "")
    }
    cat("Fitting of two distributions\n\nCall:\n")
    print(x$call)
    cat("\nParameter estimates:\nRandom:\n")
    cat_snse("p", x$param$random["prob", 1:2])
    cat("Aggregated:\n")
    cat_snse("theta", x$param$aggregated["prob", 1:2])
    cat_snse("theta", x$param$aggregated["theta", 1:2])
    cat_snse("rho",   x$param$aggregated["rho", 1:2])
    cat_snse("alpha", x$param$aggregated["alpha", 1:2])
    cat_snse("beta",  x$param$aggregated["beta", 1:2])
}

#------------------------------------------------------------------------------
#' @export
#------------------------------------------------------------------------------
summary.fit_distr <- function(object, ...) {
    res <- list(call_orig       = object$call,
                cmat_random     = object$param$random,
                cmat_aggregated = object$param$aggregated)
    structure(res, class = "summary.fit_distr")
}

#------------------------------------------------------------------------------
#' @method print summary.fit_distr
#' @export
#------------------------------------------------------------------------------
print.summary.fit_distr <- function(x, ...) {
    cat("Fitting of two distributions\n\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\nRandom:\n")
    printCoefmat(x$cmat_random)
    cat("\nAggregated:\n")
    printCoefmat(x$cmat_aggregated)
}

#------------------------------------------------------------------------------
#' @export
#------------------------------------------------------------------------------
plot.fit_distr <- function(x, ...) {

    rownames(x$freq) <- x$freq$category
    freq_long <- setNames(as.data.frame.table(as.matrix(x$freq[-1])),
                          c("Number per sampling unit", "key", "Frequency"))
    freq_long$key <- factor(tocamel(freq_long$key),
                            levels = c("Observed", "Aggregated", "Random"))

    gg <- list(
        geom_bar(data = freq_long, aes(x = `Number per sampling unit`,
                                       y = Frequency, fill = key),
                 stat = "identity", position = "dodge", color = "black",
                 size = 0.5),
        scale_fill_manual(values = c(Observed = "black", Aggregated = "grey",
                                     Random = "white")),
        theme_bw(),
        theme(legend.justification = c(0.98, 0.98), legend.position = c(0.98, 0.98),
              legend.background = element_rect(size=.5, linetype="solid", color = "black"),
              legend.title=element_blank())
    )
    ggplot() + gg
}


# fit_distr.incidence <- function(data, ..., progress = TRUE, simulatePValue) { # Extra arguments useful?
    # Add a parameter? gof = c("Chisq", "LLR", "no") ???
    #
    # ...
    #
    # if (!requireNamespace("XNomial", quietly = TRUE) && gof == "LLR") {
    #     warning("Package 'XNomial' needed to use the option 'LLR'. Auto-switch to 'Chisq'.")
    #     gof <- "Chisq"
    # }
    #
    # ...
    #

    ## cf. emdbook for 1 specific function
    ## emdbook::dbetabinom(x, prob, size,  theta, shape1, shape2, log = FALSE)

    # Et si optim n'a pas convergé ? Gérer ce cas ci-dessous.
    # Et puis, plein de warnings!!!
    # gof
    #testB <- with(freq, chisq.test(x = observed, p = binomial, rescale.p = TRUE))
    #testBBD <- with(freq, chisq.test(x = observed, p = betaBinomial, rescale.p = TRUE))
    # Add simulate.p.value = TRUE ???

    ####### Apply the Cochran rule regarding using Xchi2
    # Cochran : 80 % des classes doivent satisfaire la règle des cinq éléments tandis que les autres doivent être non vides.




