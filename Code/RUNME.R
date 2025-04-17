## Loading the algorithm
source("Algorithm.R")

## Loading the data
load("data.Rdata")

set.seed(42)

y <- matrix(0,nrow = 200, ncol = 500)
for(i in 1:200){
  y[i,] <- simulated_data$y[[i]]
}

result <- gibbs_sampler(
  y = y, 
  LA = simulated_data$L, n_iter = 1000)

result <- result[[1]]
cluster_history <- result[[2]]

