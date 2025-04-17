library(MCMCpack)
library(BayesLogit)
library(MASS)

initialize <- function(n, V, H, R) {
  list(
    C = sample(1:3, n, replace = TRUE),
    G = sample(1:H, n, replace = TRUE),
    pk = matrix(1/V, nrow = 3, ncol = V),
    nu = matrix(1/H, nrow = 3, ncol = H),
    Z = rnorm(V * (V - 1) / 2),
    X = array(rnorm(H * V * R), dim = c(H, V, R)),
    lambda = array(rexp(H * R), dim = c(H, R))
  )
}

update_pk <- function(C, y, V, alpha) {
  K <- length(unique(C))
  pk <- matrix(NA, nrow = K, ncol = V)
  for (k in 1:K) {
    yk <- y[C == k, , drop = FALSE]
    nk <- if (is.null(nrow(yk))) yk else colSums(yk)
    if (length(nk) != V) {
      nk <- rep(0, V)
    }
    pk[k, ] <- rdirichlet(1, alpha[1:V] + nk[1:V])
  }
  return(pk)
}

update_G <- function(C, LA_list, nu, pi_list) {
  n <- length(LA_list)
  G <- integer(n)
  K <- nrow(nu)
  
  for (i in 1:n) {
    cluster <- C[i]
    if (cluster > K || cluster < 1) {
      warning(paste("Invalid cluster index:", cluster, "in G update — skipping"))
      G[i] <- sample(1:length(pi_list), 1)
      next
    }
    log_probs <- sapply(1:length(nu[cluster, ]), function(h) {
      pi_h <- pi_list[[h]]
      edges <- LA_list[[i]][lower.tri(LA_list[[i]])]
      loglik <- sum(edges * log(pi_h) + (1 - edges) * log(1 - pi_h))
      log(nu[cluster, h]) + loglik
    })
    probs <- exp(log_probs - max(log_probs))
    probs <- probs / sum(probs)
    G[i] <- sample(1:length(probs), 1, prob = probs)
  }
  return(G)
}

update_nu <- function(C, G, H) {
  K <- length(unique(C))
  nu <- matrix(NA, nrow = K, ncol = H)
  for (k in 1:K) {
    gk <- G[C == k]
    nhk <- tabulate(gk, nbins = H)
    nu[k, ] <- rdirichlet(1, 1/H + nhk)
  }
  return(nu)
}

update_Z_X <- function(LA_list, G, X, lambda, Z_mu, Z_var) {
  H <- dim(X)[1]
  V <- dim(X)[2]
  R <- dim(X)[3]
  L <- V * (V - 1) / 2
  Z <- rnorm(L, mean = Z_mu, sd = sqrt(Z_var))
  for (h in 1:H) {
    for (v in 1:V) {
      X[h, v, ] <- rnorm(R, mean = 0, sd = sqrt(1 / lambda[h, ]))
    }
  }
  return(list(Z = Z, X = X))
}

update_lambda <- function(X, a1 = 2.5, a2 = 3.5) {
  H <- dim(X)[1]
  R <- dim(X)[3]
  lambda <- matrix(NA, nrow = H, ncol = R)
  for (h in 1:H) {
    theta <- rgamma(1, a1, 1)
    prod_theta <- theta
    lambda[h, 1] <- 1 / prod_theta
    for (r in 2:R) {
      theta <- rgamma(1, a2, 1)
      prod_theta <- prod_theta * theta
      lambda[h, r] <- 1 / prod_theta
    }
  }
  return(lambda)
}

update_C <- function(y, LA_list, pk, nu, alpha_c, pi_list, G) {
  n <- nrow(y)
  C_new <- integer(n)
  K <- nrow(pk)
  for (i in 1:n) {
    counts <- tabulate(C_new[-i], nbins = K)
    probs <- numeric(K + 1)
    for (k in 1:K) {
      logp_y <- sum(y[i, 1:ncol(pk)] * log(pk[k, 1:ncol(pk)]))
      logp_g <- log(nu[k, G[i]])
      probs[k] <- log(counts[k] + alpha_c) + logp_y + logp_g
    }
    probs[K + 1] <- log(alpha_c)
    probs <- exp(probs - max(probs))
    probs <- probs / sum(probs)
    C_new[i] <- sample(1:(K + 1), 1, prob = probs)
  }
  return(C_new)
}

gibbs_sampler <- function(y, LA_list, n_iter = 1000, V = 15, H = 10, R = 5) {
  n <- nrow(y)
  alpha <- rep(1, V)
  alpha_c <- 1
  state <- initialize(n, V, H, R)
  cluster_history <- list()
  for (t in 1:n_iter) {
    state$pk <- update_pk(state$C, y, V, alpha)
    K_current <- max(state$C)
    if (nrow(state$nu) < K_current) {
      n_add <- K_current - nrow(state$nu)
      state$nu <- rbind(state$nu, matrix(1 / H, nrow = n_add, ncol = H))
    }
    state$G <- update_G(state$C, LA_list, state$nu, replicate(H, runif(V * (V - 1) / 2), simplify = FALSE))
    state$nu <- update_nu(state$C, state$G, H)
    res <- update_Z_X(LA_list, state$G, state$X, state$lambda, Z_mu = 0, Z_var = 1)
    state$Z <- res$Z
    state$X <- res$X
    state$lambda <- update_lambda(state$X)
    state$C <- update_C(y, LA_list, state$pk, state$nu, alpha_c, replicate(H, runif(V * (V - 1) / 2), simplify = FALSE), state$G)
    cluster_history[[t]] <- state$C
  }
  return(list(state, cluster_history))
}