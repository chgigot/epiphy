#------------------------------------------------------------------------------#
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

# Print, summary and plot
#print.fit_distr   <- function(x, ...)
#summary.fit_distr <- function(object, ...)
#plot.fit_distr    <- function(x, ...)



#------------------------------------------------------------------------------#
# @describeIn fit_distr
# For incidence data, the two distributions fitted are binomial and beta-binomial
# distributions.
#
# @export
#------------------------------------------------------------------------------#
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

    # Initial checks and data preparation
    # d <- data$obs$d # All the raw data
    # N <- length(d)
    # n <- unique(data$obs$n)
    # if (length(n) != 1) stop(paste0("Current implementation only deals ",
    #                                 "with equal size sampling units."))





    # lp <- predict(rand$model, type = "link", se.fit = TRUE)
    # ## Ci-dessous, on utilise le fait que lp a retourné une liste
    # ## pourrait faire : lp[[residual.scale]] <- NULL
    # ## Et aussi modifyList(...)
    # lp$fit    <- lp$fit[[1]]
    # lp$se     <- lp$se.fit[[1]]
    # lp$zvalue <- lp$fit / lp$se
    # lp$pvalue <- 2 * pnorm(-abs(lp$zvalue))
    #
    # p <- predict(rand$model, type = "response", se.fit = TRUE)
    # p$fit    <- p$fit[[1]] #equivalent to former: pBinom <- fitted(modelBinom)[[1]]
    # p$se     <- p$se.fit[[1]] ## Vérifier la pertinence de cela
    # p$zvalue <- p$fit / p$se
    # p$pvalue <- 2 * pnorm(-abs(p$zvalue))
    #
    # mat <- cbind(c(lp$fit, p$fit),
    #              c(lp$se, p$se),
    #              c(lp$zvalue, p$zvalue),
    #              c(lp$pvalue, p$pvalue))
    # rownames(mat) <- c("logit(p)", "p")
    # colnames(mat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    # rand$par <- mat

    # z value (Pr(>|z|)
    # 2*pnorm(Estimate / Std. Error) = Pr(>|z|)
    # Pourquoi pas t-value??? (contrairement à lm??)
    #pvalue <- 2 * pt(-abs(tvalue), df.r)
    #coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    #dimnames(coef.table) <- list(names(coef.p), c(dn, "t value", "Pr(>|t|)"))


    # # Beta-binomial distribution
    # inits <- c(p$fit, 1 - p$fit)
    # modelBBD <- optimBetaBinom(inits, d, n) # !!!!! si retour NULL, pas d'assignement list possible # optimBetaBinom n'existe plus
    # # use stats4::mle() instead
    # agg <- list()

    # #------------------------------------------------------------------------------#
    #
    # if (!is.null(modelBBD)) { # If the optimization procedure has been successful
    #     agg$model <- modelBBD
    #     alpha <- beta <- p <- theta <- rho <- list()
    #     alpha$fit <- agg$model$par[1]
    #     beta$fit  <- agg$model$par[2]
    #     #varcov <- solve(model$hessian) # var-cov matrix
    #     #se     <- sqrt(diag(varcov))   # standard errors
    #     p$fit     <- alpha$fit / (alpha$fit + beta$fit)
    #     theta$fit <- 1 / (alpha$fit + beta$fit)
    #     rho$fit <- 1 / (alpha$fit + beta$fit + 1)
    #     # theta: index of aggregation.
    #     # Aggregation increases with increasing theta
    #
    #     # Work with the hessian
    #     estVcov  <- solve(agg$model$hessian) # estimator of the asymptotic covariance matrix
    #     estSe    <- sqrt(diag(estVcov))       # solve(mat): compute the inverse of the matrix mat
    #     alpha$se <- estSe[[1]]
    #     beta$se  <- estSe[[2]]
    #     p$se     <- msm::deltamethod(~ x1 / (x1 + x2), agg$model$par, estVcov)
    #     theta$se <- msm::deltamethod(~ 1 / (x1 + x2), agg$model$par, estVcov)
    #     rho$se   <- msm::deltamethod(~ 1 / (x1 + x2 + 1), agg$model$par, estVcov)
    #
    #     #------------------------------------------------------------------------------#
    #
    #     p$zvalue <- p$fit / p$se
    #     p$pvalue <- 2 * pnorm(-abs(p$zvalue))
    #
    #     theta$zvalue <- theta$fit / theta$se
    #     theta$pvalue <- 2 * pnorm(-abs(theta$zvalue))
    #
    #     alpha$zvalue <- alpha$fit / alpha$se
    #     alpha$pvalue <- 2 * pnorm(-abs(alpha$zvalue))
    #
    #     beta$zvalue <- beta$fit / beta$se
    #     beta$pvalue <- 2 * pnorm(-abs(beta$zvalue))
    #
    #     rho$zvalue <- rho$fit / rho$se
    #     rho$pvalue <- 2 * pnorm(-abs(rho$zvalue))



        # # Below should work with bbmle::mle2
        # # Retrieve result matrice (to which we will add extra estimates)
        # param <- coef(summary(modelBBD)) # alpha and beta
        # #mylist <- list(x1 = coef(modelBBD)[[1]], # alpha
        # #               x2 = coef(modelBBD)[[2]]) # beta
        #
        # p     <- estimate_param(modelBBD, quote(x1 / (x1 + x2)), student = FALSE)
        # theta <- estimate_param(modelBBD, quote(1  / (x1 + x2)), student = FALSE)
        # rho   <- estimate_param(modelBBD, quote(1  / (x1 + x2 + 1)), student = FALSE)
        #
        # param <- rbind(param, unlist(p), unlist(theta), unlist(rho))
        # rownames(param) <- c("alpha", "beta", "p", "theta", "rho")




        ## cf. emdbook for 1 specific function
        ## emdbook::dbetabinom(x, prob, size,  theta, shape1, shape2, log = FALSE)
        # mat <- cbind(c(NA, p$fit, theta$fit, rho$fit, alpha$fit, beta$fit),
        #              c(NA, p$se, theta$se, rho$se, alpha$se, beta$se),
        #              c(NA, p$zvalue, theta$zvalue, rho$zvalue, alpha$zvalue, beta$zvalue),
        #              c(NA, p$pvalue, theta$pvalue, rho$pvalue, alpha$pvalue, beta$pvalue))
        # rownames(mat) <- c("logit(p)", "p", "theta", "rho", "alpha", "beta")
        # colnames(mat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
        # agg$par <- mat

    # } else {
    #     agg$model <- list(NULL)
    #     freq$betaBinomial <- NA
    #     agg$par <- list(NULL)
    #     agg$test <- list(NULL)
    # }

    # Et si optim n'a pas convergé ? Gérer ce cas ci-dessous.
    # Et puis, plein de warnings!!!
    # gof
    #testB <- with(freq, chisq.test(x = observed, p = binomial, rescale.p = TRUE))
    #testBBD <- with(freq, chisq.test(x = observed, p = betaBinomial, rescale.p = TRUE))
    # Add simulate.p.value = TRUE ???

    ####### Apply the Cochran rule regarding using Xchi2
    # Cochran : 80 % des classes doivent satisfaire la règle des cinq éléments tandis que les autres doivent être non vides.

    # Returns:

    # structure(list(call = match.call(),
    #                random = list(model = rand$model,
    #                              par = rand$par,
    #                              test = NULL),
    #                aggregated = list(model = agg$model,
    #                                  par = agg$par,
    #                                  test = NULL),
    #                frequency = freq,
    #                test = llrtest(rand$model, agg$model)),
    #           class = "fit_distr") # "cmpDistr" should be better

# }

#------------------------------------------------------------------------------
#' @export
summary.fit_distr <- function(object, ...) {
    cat("\nCall:\n")
    print(object$call)
    cat("\n----------\nRandom:\n")
    printCoefmat(res$random$par)
    cat("\n----------\nAggregated:\n")
    printCoefmat(res$aggregated$par)
    cat("\n----------\n")
    print(object$test)
}

#' @export
print.fit_distr <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    #cat("\nSource: <", x, ">\n", sep = "")
    cat("\nParameters:\n")
    cat("Random (BD distr.):\n")
    cat("  p     = ", formatC(x$random$par["p", 1]),
        " (\u00b1 ", formatC(x$random$par["p", 2]), ")\n", sep = "" )
    cat("Aggregated (BBD distr.):\n")
    cat("  p     = ", formatC(x$aggregated$par["p", 1]),
        " (\u00b1 ", formatC(x$aggregated$par["p", 2]), ")\n", sep = "")
    cat("  theta = ", formatC(x$aggregated$par["theta", 1]),
        " (\u00b1 ", formatC(x$aggregated$par["theta", 2]), ")\n", sep = "")
    cat("  rho   = ", formatC(x$aggregated$par["rho", 1]),
        " (\u00b1 ", formatC(x$aggregated$par["rho", 2]), ")\n", sep = "")
    cat("  alpha = ", formatC(x$aggregated$par["alpha", 1]),
        " (\u00b1 ", formatC(x$aggregated$par["alpha", 2]), ")\n", sep = "")
    cat("  beta  = ", formatC(x$aggregated$par["beta", 1]),
        " (\u00b1 ", formatC(x$aggregated$par["beta", 2]), ")\n", sep = "")
}

#' @export
plot.fit_distr <- function(x, y, ..., plot = TRUE) {
    g <- list(
        #geom_bar(data = (x$frequency %>% tidyr::gather("key", "value", -category)),
        #         aes(x = category, y = value, fill = key),
        #         stat = "identity", position = "dodge")
        geom_bar(data = (x$frequency %>% tidyr::gather("key", "value", -category)),
                 aes(x = category, y = value, fill = key),
                 stat = "identity", position = "dodge", colour = "black", size = 0.5),
        scale_fill_manual(breaks = c("observed", "betaBinomial", "binomial"),
                          values = c(observed = "black", betaBinomial = "grey", binomial = "white"))#,
        #                          guide = FALSE)
        #scale_fill_discrete(guide = FALSE)

        #    ggplot(data = (x$freqTable %>% tidyr::gather("key", "value", -category)),
        #           aes(x = category, y = value, fill = key)) +
        #        geom_bar(stat = "identity", position = "dodge", colour = "black", size = 0.5) +
        #        scale_fill_manual(breaks = c("observed", "betaBinomial", "binomial"),
        #                          values = c(observed = "black", betaBinomial = "grey", binomial = "white"),
        #                          guide = FALSE)
        #scale_fill_discrete(guide = FALSE)
    )
    if (!plot) g else ggplot() + g
}


#plot.fit_distr <- function(x, y, ...,
#                          type = c("all", "distribution", "parameters"),
#                          parameters = c("p-theta", "p-rho", "alpha-beta"),
#                          trend = TRUE) {}







