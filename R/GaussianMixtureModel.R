#' Gaussian Mixture Model
#'
#' Fit a Gaussian Mixture Model by the EM Algorithm
#'
#' @import mvtnorm
#' @param X an N by p data matrix
#' @param k integer, number of Gaussian components
#' @param initial_mu a k by p matrix, the initial points of the EM Algorithm; each row corresponds to one component
#' @param max_iter integer, the maximum number of EM iterations
#' @param tol scaler or NULL, used as the convergence criteria \cr
#' if tol is scaler, if the distance between centroids in the current iteration and the previous iteration is less than tol, then the algorithm stops early \cr
#' if tol is NULL, then early stopping criterion does not apply
#'
#' @return A list containing mu, Sigma, r and loglikelihood \cr
#' mu: a k by p matrix, fitted centroids; each row corresponds to one component \cr
#' Sigma: a p by p by k array, fitted covariance matrices \cr
#' r: an N by k matrix, fitted responsibility; each row corresponds to one sample; the value (i, j) is the fitted probability of sample i belonging to component j \cr
#'
#' @examples
#' # create samples from a two-component Gaussian Mixture
#' X = matrix(0, nrow = 500, ncol = 2)
#' z = sample(2, 500, replace = TRUE)
#' X[which(z == 1), ] = rnorm(sum(z == 1) * 2, mean = 0, sd = 1)
#' X[which(z == 2), ] = rnorm(sum(z == 2) * 2, mean = 5, sd = 1)
#'
#' # fit a Gaussian Mixture Model to the data
#' gmm = GaussianMixtureModel(X, 2, initial_mu = matrix(c(0, 1, 0, 1), nrow = 2, ncol = 2))
#'
#' # prediction
#' centroids = gmm$mu
#' r = gmm$r
#' cluster = apply(r, 1, which.max)
#'
#' # plot
#' plot(X[, 1], X[, 2], col = c('orange', 'red')[cluster], cex = 0.5)
#' points(x = centroids[, 1], y = centroids[, 2], pch = 10, col = 'blue', cex = 2)
#' legend(x = 'topleft',
#'        legend = c('cluster 1', 'cluster 2', 'centers'),
#'        col = c('orange', 'red', 'blue'),
#'        pch = c(1, 1, 10))
#'
#' @export
GaussianMixtureModel = function(X, k, initial_mu, max_iter = 100, tol = 1e-8) {
  N = dim(X)[1]     # number of samples
  p = dim(X)[2]     # number of features
  loglik_list = rep(0, max_iter)

  # initialization
  pi = rep(1 / k, k)
  mu = initial_mu
  Sigma = array(rep(diag(p), k), dim = c(p, p, k))    # stack k identity matrices on each other

  # EM Algorithm
  for (i in 1:max_iter) {
    mu_prev = mu

    # update responsibility
    r = matrix(0, nrow = N, ncol = k)    # initialize the responsibility
    for (l in 1:k) {
      r[, l] = log(pi[l]) + dmvnorm(X, mean = mu[l, ], sigma = Sigma[, , l], log = TRUE)
    }
    min.logr = apply(r, 1, min)
    r = r - min.logr    # for numerical stability
    r = exp(r)
    r = r / rowSums(r)

    # update N_k
    N_k = colSums(r)

    # update pi
    pi = N_k / N

    # update mu
    for (l in 1:k) {
      mu[l, ] = colSums(r[, l] * X) / N_k[l]
    }

    # update sigma
    for (l in 1:k) {
      s = 0
      for (n in 1:N) {
        d = X[n, ] - mu[l, ]
        s = s + r[n, l] * outer(d, d)
      }
      Sigma[, , l] = s / N_k[l]
    }

    # check for convergence
    if (length(tol) > 0) {
      if (sqrt( sum((mu - mu_prev)^2) ) < tol) {
        cat('converges at iteration', i, '\n')
        break
      }
    }
  }

  return(list(mu = mu,
              Sigma = Sigma,
              r = r))
}
