# Gaussian Mixture Model
#
# Fit a Gaussian Mixture Model by the EM Algorithm
#
#'@import mvtnorm
#'
#'@para X: an N by p data matrix
#'@para k: integer, number of Gaussian components
#'@para initial_mu: a k by p matrix, the initial points of the EM Algorithm; each row corresponds to one component
#'@para max_iter: integer, the maximum number of EM iterations
#'@para tol: scaler, used as the convergence criteria; if the distance between centroids in the current iteration and the previous iteration is less than tol, then the algorithm stops early
#'
#'@return mu: a k by p matrix, fitted centroids; each row corresponds to one component
#'@return Sigma: a p by p by k array, fitted covariance matrices
#'@return r: an N by k matrix, fitted responsibility; each row corresponds to one sample; the value (i, j) is the fitted probability of sample i belonging to component j
#'@return loglikelihood: a numeric vector recording loglikelihood throughout iterations
#'
#'@export
GaussianMixtureModel = function(X, k, initial_mu, max_iter = 100, tol = 1e-8) {
  # function used to compute loglikelihood
  loglikelihood = function(X, pi, mu, Sigma) {
    N = dim(X)[1]
    k = length(pi)
    llk = matrix(0, nrow = N, ncol = k)
    for (l in 1:2) {
      llk[, l] = pi[l] * dmvnorm(X, mean = mu[l, ], sigma = Sigma[, , l])
    }
    value = sum(log(rowSums(llk)))
    return(value)
  }

  N = dim(X)[1]     # number of samples
  p = dim(X)[2]     # number of features
  loglik_list = rep(0, max_iter)

  # initialization
  pi = rep(1 / k, k)
  mu = initial_mu
  Sigma = array(rep(diag(p), k), dim = c(p, p, k))    # stack k identity matrices on each other

  # EM Algorithm
  for (i in 1:max_iter) {
    #browser()
    mu_prev = mu

    # update responsibility
    r = matrix(0, nrow = N, ncol = k)    # initialize the responsibility
    for (l in 1:k) {
      r[, l] = pi[l] * dmvnorm(X, mean = mu[l, ], sigma = Sigma[, , l])
    }
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

    # compute loglikelihood
    loglik_list[i] = loglikelihood(X, pi, mu, Sigma)

    # check for convergence
    if (sqrt( sum((mu - mu_prev)^2) ) < tol) {
      cat('converge at iteration', i)
      loglik_list = loglik_list[1:i]
      break
    }
  }

  return(list(mu = mu,
              Sigma = Sigma,
              r = r,
              loglikelihood = loglik_list))
}
