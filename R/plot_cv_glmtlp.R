#' Plot Method for a "cv.glmtlp" Object
#'
#' @description 
#' Plots the cross-validation curve, and the upper and lower standard deviation
#'   curves, as a function of the \code{lambda} or \code{kappa} values. 
#'
#' @details 
#' The generated plot is a \code{ggplot} object, and therefore, the users are able 
#'   to customize the plots following the \code{ggplot2} syntax.
#'
#' @aliases plot.cv.glmtlp
#' @param x Fitted \code{cv.glmtlp} object
#' @param vertical.line Logical. Whether or not include a vertical line 
#'   indicating the position of the index which gives the smallest CV error.
#' @param \dots Additional arguments.
#' 
#' @author Chunlin Li, Yu Yang, Chong Wu 
#'   \cr Maintainer: Yu Yang \email{yang6367@umn.edu}
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
#' @keywords models plot
#'
#' @examples
#' X <- matrix(rnorm(100 * 20), 100, 20)
#' y <- rnorm(100)
#' cv.fit <- cv.glmtlp(X, y, family = "gaussian", penalty = "tlp")
#' plot(cv.fit)
#' plot(cv.fit, vertical.line = FALSE)
#' cv.fit2 <- cv.glmtlp(X, y, family = "gaussian", penalty = "l0")
#' plot(cv.fit2)
#' plot(cv.fit2, vertical.line = FALSE)
#' 
#' data("gau_data")
#' cv.fit <- cv.glmtlp(gau_data$X, gau_data$y, family = "gaussian", penalty = "tlp")
#' plot(cv.fit)
#' 
#' data("bin_data")
#' cv.fit <- cv.glmtlp(bin_data$X, bin_data$y, family = "binomial", penalty = "l1")
#' plot(cv.fit)
#' 
#' @import ggplot2
#' @method plot cv.glmtlp
#' @export
#' @export plot.cv.glmtlp


plot.cv.glmtlp <- function(x, vertical.line=TRUE, ...) {
  if (x$fit$penalty == "l0") {
    xlab <- expression(kappa)
    index <- x$kappa
    min.index <- x$kappa.min
  } else {
    xlab <- expression(Log(lambda))
    index <- log(x$lambda)
    min.index <- log(x$lambda.min)
  }
  if (length(x$cv.mean) <= 1) {
    warning("No plot produced, since the length of cv sequence is <= 1.")
    return ()
  }
  df <- data.frame(matrix(nrow=length(index), ncol = 3))
  df$index <- index
  df$cvm <- x$cv.mean
  df$cvl <- x$cv.mean - x$cv.se
  df$cvu <- x$cv.mean + x$cv.se
  ylab <- "Deviance"
  title <- "Cross-validation Plot"
  g <- ggplot(df, aes_string(x = "index", y = "cvm")) + 
    xlab(xlab) + ylab(ylab) + ggtitle(title) + 
    geom_line(lwd=1, color = "red") + geom_point(pch=19, color = "red") + 
    theme_bw() + theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + 
    geom_errorbar(aes_string(ymin = "cvl", ymax = "cvu"), width=0.02, color = "gray")
  
  if (vertical.line) {
    g <- g + geom_vline(xintercept = min.index, linetype = "dashed")
  }
  if (x$fit$penalty == "l0") {
    g <- g + scale_x_discrete(limits = factor(x$kappa))
  }
  g
}

