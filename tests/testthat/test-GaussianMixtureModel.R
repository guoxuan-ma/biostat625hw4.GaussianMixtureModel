test_that("GaussianMixtureModel works", {
  set.seed(2021)
  for (i in 1:20) {
    mean_1 = runif(1, -10, 0)
    mean_2 = runif(1, 5, 15)
    X = matrix(0, nrow = 500, ncol = 2)
    z = sample(2, 500, replace = TRUE)
    X[which(z == 1), ] = rnorm(sum(z == 1) * 2, mean = mean_1, sd = 1)
    X[which(z == 2), ] = rnorm(sum(z == 2) * 2, mean = mean_2, sd = 1)
    gmm = GaussianMixtureModel(X, 2, initial_mu = matrix(c(mean_1, mean_2, mean_1, mean_2), nrow = 2, ncol = 2))
    cluster = apply(gmm$r, 1, which.max)
    acc = max(sum(cluster == z), sum(cluster != z)) / 500

    expect_equal(sum((gmm$mu - c(mean_1, mean_2, mean_1, mean_2))^2) < 0.1, TRUE)
    expect_equal(acc > 0.99, TRUE)
  }
})
