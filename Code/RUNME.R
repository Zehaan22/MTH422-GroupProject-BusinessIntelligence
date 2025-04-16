## Loading the algorithm
source("Algorithm.R")

set.seed(42)
#final_state <- gibbs_sampler(y = your_y_data, LA = your_LA_data, n_iter = 1000)


### Simulating the data for testing
simulate_data <- function(n = 50, V = 15, K = 3, H = 5, R = 3, ni = 500) {
  # 1. Simulate cluster-specific monoproduct choice probabilities
  pk_list <- lapply(1:K, function(k) {
    base <- rep(1, V)
    tweak <- sample(1:V, 2)
    base[tweak] <- base[tweak] + runif(2, 3, 5)
    base / sum(base)
  })
  
  # 2. Assign clusters to agencies
  C_true <- sample(1:K, n, replace = TRUE)
  
  # 3. Generate monoproduct counts matrix y
  y <- matrix(0, nrow = n, ncol = V)
  for (i in 1:n) {
    y[i, ] <- rmultinom(1, size = ni, prob = pk_list[[C_true[i]]])
  }
  
  # 4. Simulate H latent eigenmodels (π.h/)
  make_pi_h <- function(V, R) {
    X <- matrix(rnorm(V * R), nrow = V, ncol = R)
    Lambda <- diag(rexp(R, 1))
    S <- X %*% Lambda %*% t(X)
    Z <- rnorm(V * (V - 1) / 2, mean = 0, sd = 1)
    logits <- Z + lower_tri_vec(S)
    1 / (1 + exp(-logits))
  }
  
  lower_tri_vec <- function(M) {
    M[lower.tri(M)]
  }
  
  pi_list <- replicate(H, make_pi_h(V, R), simplify = FALSE)
  
  # 5. Cluster-specific ν_k vectors (mixing probabilities over H components)
  nu_list <- lapply(1:K, function(k) rdirichlet(1, rep(1, H))[1, ])
  
  # 6. Sample G_i and generate co-subscription networks
  G_true <- integer(n)
  LA <- matrix(0, nrow = n, ncol = V * (V - 1) / 2)
  
  for (i in 1:n) {
    k <- C_true[i]
    G_true[i] <- sample(1:H, 1, prob = nu_list[[k]])
    pi <- pi_list[[G_true[i]]]
    LA[i, ] <- rbinom(length(pi), 1, prob = pi)
  }
  
  return(list(y = y, LA = LA, C_true = C_true, G_true = G_true, pk_list = pk_list))
}

### Running the sampler
sim_data <- simulate_data(n = 50, V = 15, K = 3, H = 5, R = 3, ni = 500)

# Run Gibbs Sampler
result <- gibbs_sampler(y = sim_data$y, LA = sim_data$LA, n_iter = 100)
cluster_history <- result[[2]]
result <- result[[1]]
