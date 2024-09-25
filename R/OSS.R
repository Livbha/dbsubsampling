## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
##
## R main function `OSS`.  `Rcpp`-version by `the package itself`.
##
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------

#' Orthogonal subsampling for big data linear regression (OSS, `Rcpp`-version by `the package itself`)
#'
#' A subsampling method based on orthogonal array for linear model.
#'
#' @param n Subsample size.
#' @param X A matrix or data frame.
#'
#' @return Subsample index.
#'
#' @examples
#' data_numeric_regression["y"] <- NULL
#' X <- as.matrix(data_numeric_regression)
#' OSS(100, X)
#'
#' @references Lin Wang, Jake Elmstedt, Weng Kee Wong & Hongquan Xu (2021)
#' \emph{Orthogonal subsampling for big data linear regression,
#' The Annals of Applied Statistics, 15(3), 1273-1290},
#' \url{https://projecteuclid.org/journals/annals-of-applied-statistics/volume-15/issue-3/Orthogonal-subsampling-for-big-data-linear-regression/10.1214/21-AOAS1462.short?tab=ArticleLink}.
#'
#' @export

scale_neg_pos_1 <- function(X) {
 apply(X,2,function(.col)((.col-mean(.col))/max(max(.col)-mean(.col),mean(.col)-min(.col))))
}
        
OSS <- function(n, X, k){
  X <- scale_neg_pos_1(as.matrix(X)) # need scale
  attributes(X) <- attributes(X)["dim"]
  subindex <- rcppOSS(X = X, n = n, K=k)
  return(subindex)
}

## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
##
## R main function `myArma_OSS`.  `RcppArmadillo`-version by `Zhu`.
##
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------

#' OSS (`RcppArmadillo`-version by `Zhu`)
#'
#' @param n Subsample size.
#' @param X A matrix.
#'
#' @return Subsample index.
#' @export
# @examples
# data_numeric_regression["y"] <- NULL
# X <- as.matrix(data_numeric_regression)
# myR_OSS(X, 100)
myArma_OSS <- function(n, X){
  X <- scale_neg_pos_1(as.matrix(X))
  attributes(X) <- attributes(X)["dim"]
  as.vector(armaOSS(X, n))
}


## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------
##
## R main function `myR_OSS`.  `base-R` by `the package itself`.
##
## ---------------------------------------------------------------------------
## ---------------------------------------------------------------------------

#' Get L2 norm (r-version)
#'
#' Get L2 norm of a matrix or data frame.
#' @param X A matrix or data.frame.
#'
#' @return L2 norm of `X`(every row).
#'
# @examples
# X <- matrix(1:12, 4, 3)
# X <- scale(X)
# rL2norm(X)
rL2norm <- function(X) {
  return(rowSums(X^2))
}

#' Compute loss function for OSS (r-version)
#'
#' @param candi The index of the candidate set.
#' @param last_index The index of the seleted point in last iteration.
#' @param X The whole data.
#' @param norm Norm of the whole data.
#' @param p Numbers of columns of the data.
#'
#' @return Loss of every point in candidate set.
# @examples
# X <- matrix(1:20, 5, 4)
# X <- scale(X)
# norm <- rL2norm(X)
# rComputeLoss(c(1,3,4), 2, X, norm)
rComputeLoss <- function(candi, last_index, X, norm, p = ncol(X)){
  delta <- rowSums(t(apply(X[candi, ], 1, function(.row) sign(.row) == sign(X[last_index,]))))
  loss <- (p - norm[candi]/2  -  norm[last_index]/2 + delta)^2
  return(loss)
}

#' Find t smallest index of a vector.
#'
#' @param loss A vector.
#' @param t A int
#'
#' @return The index of the t smallest element of the vector.
#'
# @examples
# loss <- rnorm(10)
# rbottom_t_index(loss, 3)
rbottom_t_index <- function(loss, t){
  return(which(loss <= sort(loss)[t]))
}


#' OSS with `base R` by `the package itself`.
#'
#' @param n Subsample size.
#' @param X A matrix.
#'
#' @return Subsample index.
#' @export
# @examples
# data_numeric_regression["y"] <- NULL
# X <- as.matrix(data_numeric_regression)
# myR_OSS(X, 100)
myR_OSS <- function(n, X){
  X <- scale_neg_pos_1(as.matrix(X))
  attributes(X) <- attributes(X)["dim"]
  N <- nrow(X)

  index <- numeric(n)
  candi <- 1:N

  norm <- rL2norm(X)
  r <- log(N)/log(n)

  for (i in 1:n) {
    # Initial
    if (i == 1) {
      index[1] <- which.max(norm)
      candi <- candi[-index[1]]
      loss <- rComputeLoss(candi, index[1], X, norm)
      next
    }

    # Election
    tmp <- which.min(loss)
    index[i] <- candi[tmp]
    candi <- candi[-tmp]
    loss <- loss[-tmp]

    # Elimination
    t <- ifelse(N > (n^2), N/i, N/(i^(r-1)))
    if (length(candi) > t) {
      candi <- candi[rbottom_t_index(loss,t)]
      loss <- loss[rbottom_t_index(loss,t)]
    }

    #    if (length(candi) == 0) {
    #      index <- index[1:i]
    #      break
    #    }
    # Update loss
    loss <- loss + rComputeLoss(candi, index[i], X, norm)
  }

  return(index)
}
