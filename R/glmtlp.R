#' Fit a GLM with L0, L1, or TLP Penalization
#'
#' @description
#' Fit generalized linear models via penalized maximum likelihood. The
#'   regularization path is computed for the l0, lasso, or truncated lasso
#'   penalty at a grid of values for the regularization parameter \code{lambda}
#'   or \code{kappa}. Fits linear and logistic regression models.
#'
#' @details
#' The sequence of models indexed by \code{lambda} (when \code{penalty = c('l1', 'tlp')})
#'   or \code{kappa} (when \code{penalty = 'l0'}) is fit by the coordinate
#'   descent algorithm.
#'
#'   The objective function for the \code{"gaussian"} family is:
#'   \deqn{1/2 RSS/nobs + \lambda*penalty,} and for the other models it is:
#'   \deqn{-loglik/nobs + \lambda*penalty.}
#'   Also note that, for \code{"gaussian"}, \code{glmtlp} standardizes y to
#'   have unit variance (using 1/(n-1) formula).
#'
#'   ## Details on \code{family} option
#'
#'   \code{glmtlp} currently only supports built-in families, which are specified by a
#'   character string. For all families, the returned object is a regularization
#'   path for fitting the generalized linear regression models, by maximizing the
#'   corresponding penalized log-likelihood. \code{glmtlp(..., family="binomial")}
#'   fits a traditional logistic regression model for the log-odds.
#'
#'   ## Details on \code{penalty} option
#'
#'   The built-in penalties are specified by a character string. For \code{l0}
#'   penalty, \code{kappa} sequence is used for generating the regularization
#'   path, while for \code{l1} and \code{tlp} penalty, \code{lambda} sequence
#'   is used for generating the regularization path.
#'
#' @param X Input matrix, of dimension \code{nobs} x \code{nvars};
#'   each row is  an observation vector.
#' @param y Response variable, of length \code{nobs}. For \code{family="gaussian"},
#'   it should be quantitative; for \code{family="binomial"}, it should be either
#'   a factor with two levels or a binary vector.
#' @param family A character string representing one of the built-in families.
#'   See Details section below.
#' @param penalty A character string representing one of the built-in penalties.
#'   \code{"l0"} represents the \eqn{L_0} penalty, \code{"l1"} represents the
#'   lasso-type penalty (\eqn{L_1} penalty), and \code{"tlp"} represents the
#'   truncated lasso penalty.
#' @param nlambda The number of \code{lambda} values. Default is 100.
#' @param lambda.min.ratio The smallest value for \code{lambda}, as a fraction of
#'   \code{lambda.max}, the smallest value for which all coefficients are zero.
#'   The default depends on the sample size \code{nobs} relative to the number
#'   of variables \code{nvars}. If \code{nobs > nvars}, the default is
#'   \code{0.0001}, and if \code{nobs < nvars}, the default is \code{0.01}.
#' @param lambda A user-supplied \code{lambda} sequence. Typically, users should let
#'   the program compute its own \code{lambda} sequence based on
#'   \code{nlambda} and \code{lambda.min.ratio}. Supplying a value of
#'   \code{lambda} will override this. WARNING: please use this option with care.
#'   \code{glmtlp} relies on warms starts for speed, and it's often faster to
#'   fit a whole path than a single fit. Therefore, provide a decreasing sequence
#'   of \code{lambda} values if you want to use this option. Also, when
#'   \code{penalty = 'l0'}, it is not recommended for the users to supply
#'   this parameter.
#' @param kappa A user-supplied \code{kappa} sequence. Typically, users should
#'   let the program compute its own \code{kappa} sequence based on \code{nvars}
#'   and \code{nobs}. This sequence is used when \code{penalty = 'l0'}.
#' @param tau A tuning parameter used in the TLP-penalized regression models.
#'   Default is  \code{0.3 * sqrt(log(nvars)/nobs)}.
#' @param delta A tuning parameter used in the coordinate majorization descent
#'   algorithm. See Yang, Y., & Zou, H. (2014) in the reference for more detail.
#' @param tol Tolerance level for all iterative optimization algorithms.
#' @param weights Observation weights. Default is 1 for each observation.
#' @param penalty.factor Separate penalty factors applied to each
#'   coefficient, which allows for differential shrinkage. Default is 1
#'   for all variables.
#' @param standardize Logical. Whether or not standardize the input matrix
#'   \code{X}; default is \code{TRUE}.
#' @param dc.maxit Maximum number of iterations for the DC (Difference of
#'   Convex Functions) programming; default is 20.
#' @param cd.maxit Maximum number of iterations for the coordinate descent
#'   algorithm; default is 10^4.
#' @param nr.maxit Maximum number of iterations for the Newton-Raphson method;
#'   default is 500.
#' @param ... Additional arguments.
#' @return An object with S3 class \code{"glmtlp"}.
#'
#' \item{beta}{a \code{nvars x length(kappa)} matrix of
#'   coefficients when \code{penalty = 'l0'}; or a \code{nvars x length(lambda)}
#'   matrix of coefficients when \code{penalty = c('l1', 'tlp')}.}
#' \item{call}{the call that produces this object.}
#' \item{family}{the distribution family used in the model fitting.}
#' \item{intercept}{the intercept vector, of \code{length(kappa)} when
#'   \code{penalty = 'l0'} or \code{length(lambda)} when
#'   \code{penalty = c('l1', 'tlp')}.}
#' \item{lambda}{the actual sequence of \code{lambda} values used. Note that
#'   the length may be smaller than the provided \code{nlambda} due to removal
#'   of saturated values.}
#' \item{penalty}{the penalty type in the model fitting.}
#' \item{penalty.factor}{the penalty factor for each coefficient used in the model fitting.}
#' \item{tau}{the tuning parameter used in the model fitting, available when
#'   \code{penalty = 'tlp'}.}
#'
#' @author Chunlin Li, Yu Yang, Chong Wu
#'   \cr Maintainer: Yu Yang \email{yang6367@umn.edu}
#'
#' @seealso \code{print}, \code{predict}, \code{coef} and \code{plot} methods,
#' and the \code{cv.glmtlp} function.
#'
#' @references Shen, X., Pan, W., & Zhu, Y. (2012).
#'   \emph{Likelihood-based selection and sharp parameter estimation.
#'   Journal of the American Statistical Association, 107(497), 223-232.}
#'   \cr Shen, X., Pan, W., Zhu, Y., & Zhou, H. (2013).
#'   \emph{On constrained and regularized high-dimensional regression.
#'   Annals of the Institute of Statistical Mathematics, 65(5), 807-832.}
#'   \cr Li, C., Shen, X., & Pan, W. (2021).
#'   \emph{Inference for a Large Directed Graphical Model with Interventions.
#'   arXiv preprint arXiv:2110.03805.}
#'   \cr Yang, Y., & Zou, H. (2014).
#'   \emph{A coordinate majorization descent algorithm for l1 penalized learning.
#'   Journal of Statistical Computation and Simulation, 84(1), 84-95.}
#'   \cr Two R package Github: \emph{ncvreg} and \emph{glmnet}.
#'
#' @keywords models regression
#'
#' @examples
#'
#' # Gaussian
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' fit1 <- glmtlp(X, y, family = "gaussian", penalty = "l0")
#' fit2 <- glmtlp(X, y, family = "gaussian", penalty = "l1")
#' fit3 <- glmtlp(X, y, family = "gaussian", penalty = "tlp")
#'
#' # Binomial
#'
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- sample(c(0, 1), 100, replace = TRUE)
#' fit <- glmtlp(X, y, family = "binomial", penalty = "l1")
#' @importFrom stats model.matrix
#' @export glmtlp

glmtlp <- function(X, y, family = c("gaussian", "binomial"), penalty = c("l0", "l1", "tlp"),
                   nlambda = ifelse(penalty == "l0", 50, 100), 
                   lambda.min.ratio = ifelse(nobs < nvars, 5e-2, 1e-3),
                   lambda = NULL, kappa = NULL,
                   tau = 0.3 * sqrt(log(nvars) / nobs), delta = 2.0, tol = 1e-4,
                   weights = NULL, penalty.factor = rep(1.0, nvars), standardize = FALSE,
                   dc.maxit = 20, cd.maxit = 10000, nr.maxit = 20, ...) {
  # Coersion
  this.call <- match.call()
  family <- match.arg(family)
  penalty <- match.arg(penalty)

  # check X
  xdim <- dim(X)
  if (is.null(xdim) | (xdim[2] <= 1)) stop("X should be a matrix with 2 or more columns")
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~ 0 + ., data = X), silent = TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (typeof(X) == "character") stop("X must be a numeric matrix")
  if (typeof(X) == "integer") storage.mode(X) <- "double"

  nobs <- as.integer(xdim[1])
  nvars <- as.integer(xdim[2])

  # check y

  y <- drop(y) # Delete the dimensions of an array which has only one level.
  if (length(y) != nobs) stop(paste("number of observations in y (", length(y), ") not equal to the number of rows of X (", nobs, ")", sep = ""))
  if (!is.double(y)) {
    op <- options(warn = 2) # turn warnings into errors
    on.exit(options(op))
    y <- tryCatch(
      error = function(cond) stop("y must be numeric or able to be coerced to numeric"),
      as.double(y)
    )
    options(op)
  }
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to the model")
  if (family == "binomial" & length(table(y)) > 2) stop("Attemping to use family='binomial' with non-binary data", call. = FALSE)
  if (family == "binomial" & !identical(sort(unique(y)), 0:1)) {
    y <- as.double(y == max(y))
  }

  # check penalty.factor and weights
  penalty.factor <- as.double(penalty.factor)
  if (length(penalty.factor) != nvars) {
    stop(paste("the length of penalty.factor (",
      length(penalty.factor), ") not equal to the number of variables (",
      nvars, ")",
      sep = ""
    ))
  }
  if (is.null(weights)) {
    weights <- rep(1.0, nobs)
  } else if (length(weights) != nobs) {
    stop(paste("number of elements in weights (",
      length(weights), ") not equal to the number of rows of X (",
      nobs, ")",
      sep = ""
    ))
  } else {
    weights <- as.double(weights)
    weights <- weights * nobs / sum(weights)
  }

  ## check/setup lambda and kappa
  if (is.null(lambda)) {
    nlambda <- as.integer(nlambda)
    if (nlambda < 2) stop("nlambda must be at least 2")
    if (lambda.min.ratio >= 1) stop("lambda.min.ratio should be less than 1")
    lambda <- setup_lambda(X, y, weights, lambda.min.ratio, nlambda)
  } else {
    nlambda <- as.integer(length(lambda))
    if (nlambda < 1) stop("the length of input lambda must be at least 1")
    if (any(lambda < 0)) stop("lambdas should be non-negative")
    if (!is.double(lambda)) lambda <- as.double(lambda)
    if (nlambda == 1) {
      lambda_max <- get_lambda_max(X, y, weights)
      if (lambda[0] >= lambda_max) {
        return(get_null_output(
          this.call,
          y, family, penalty, penalty.factor, weights
        ))
      } else {
        lambda <- c(lambda_max, lambda)
      }
    }
    lambda <- sort(lambda, decreasing = TRUE)
  }

  if (penalty == "l0") {
    if (is.null(kappa)) {
      kappa <- 1:min(nvars, max(10, as.integer(nobs / log(nvars))))
      nkappa <- as.integer(length(kappa))
    } else {
      if (any(kappa < 0)) stop("kappa should be non-negative")
      if (!is.integer(kappa)) {
        message("Coercing kappa to integers.")
        kappa <- as.integer(kappa)
      }
      kappa <- sort(kappa, decreasing = FALSE)
      nkappa <- as.integer(length(kappa))
    }
  }

  ## check tau, delta, tol, and maxiters:
  ## may add more on tau, delta, and tol checking
  tau <- as.double(tau)
  delta <- as.double(delta)
  tol <- as.double(tol)
  dc.maxit <- as.integer(dc.maxit)
  cd.maxit <- as.integer(cd.maxit)
  nr.maxit <- as.integer(nr.maxit)

  standardize <- as.integer(standardize)

  ## fit
  if (family == "gaussian") {
    fit <- switch(penalty,
      "l0" = .Call(
        "gaussian_l0",
        y, X, weights, penalty.factor, kappa, lambda, tau, delta,
        standardize, tol, cd.maxit, dc.maxit
      ),
      "l1" = .Call(
        "gaussian_l1",
        y, X, weights, penalty.factor, lambda, delta,
        standardize, tol, cd.maxit
      ),
      "tlp" = .Call(
        "gaussian_tlp",
        y, X, weights, penalty.factor, lambda, tau, delta,
        standardize, tol, cd.maxit, dc.maxit
      )
    )
  } else if (family == "binomial") {
    fit <- switch(penalty,
      "l0" = .Call(
        "logistic_l0",
        y, X, weights, penalty.factor, kappa, lambda, tau, delta,
        standardize, tol, nr.maxit, cd.maxit, dc.maxit
      ),
      "l1" = .Call(
        "logistic_l1",
        y, X, weights, penalty.factor, lambda, delta,
        standardize, tol, nr.maxit, cd.maxit
      ),
      "tlp" = .Call(
        "logistic_tlp",
        y, X, weights, penalty.factor, lambda, tau, delta,
        standardize, tol, nr.maxit, cd.maxit, dc.maxit
      )
    )
  }

  dim(fit[[2]]) <- c(nvars, ifelse(penalty == "l0", nkappa, nlambda))
  varnames <- colnames(X)
  if (is.null(varnames)) varnames <- paste("V", seq(nvars), sep = "") 
  rownames(fit[[2]]) <- varnames
  if (penalty == "l0") {
    colnames(fit[[2]]) <- paste(kappa)
  } else {
    colnames(fit[[2]]) <- lambda_names(lambda)
  }

  out <- structure(list(
    beta = fit[[2]],
    call = this.call,
    family = family,
    intercept = fit[[1]],
    penalty = penalty,
    penalty.factor = penalty.factor,
    weights = weights,
    deviance = fit[[3]]
  ),
  class = "glmtlp"
  )
  if (penalty == "l0") {
    out$lambda <- lambda
    out$kappa <- kappa
    out$tau <- tau
  } else if (penalty == "l1") {
    out$lambda <- lambda
  } else {
    out$lambda <- lambda
    out$tau <- tau
  }
  out
}


get_null_output <- function(this.call, y, family, penalty,
                            penalty.factor, weights) {
  structure(list(
    beta = 0,
    call = this.call,
    family = family,
    intercept = weighted.mean(y, weights),
    penalty = penalty,
    penalty.factor = penalty.factor,
    weights = weights
  ),
  class = "glmtlp"
  )
}