#' Generate lambda sequence.
#'
#' @param X Input matrix, of dimension \code{nobs} x \code{nvars}; 
#'   each row is  an observation vector. 
#' @param y Response variable, of length \code{nobs}. For \code{family="gaussian"}, 
#'   it should be quantitative; for \code{family="binomial"}, it should be either 
#'   a factor with two levels or a binary vector. 
#' @param weights Observation weights.
#' @param lambda.min.ratio The smallest value for \code{lambda}, as a fraction of 
#'   \code{lambda.max}, the smallest value for which all coefficients are zero. 
#'   The default depends on the sample size \code{nobs} relative to the number 
#'   of variables \code{nvars}. 
#' @param nlambda The number of \code{lambda} values.
#' 
#' @importFrom stats weighted.mean
#' 
setup_lambda <- function(X, y, weights, lambda.min.ratio, nlambda) {
    lambda.max <- get_lambda_max(X, y, weights)

    lambda <- exp(seq(
        from = log(lambda.max),
        to = log(lambda.min.ratio * lambda.max),
        length.out = nlambda
    ))
    lambda
}

get_lambda_max <- function(X, y, weights) {
    rw <- (y - weighted.mean(y, weights)) * weights
    max(abs(crossprod(X, rw)), na.rm = TRUE) / nrow(X)
}
