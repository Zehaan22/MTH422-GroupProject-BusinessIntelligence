## Loading the generated sampels
source("RUNME.R")

plot_cluster_evolution <- function(cluster_history) {
  num_clusters <- sapply(cluster_history, function(C) length(unique(C)))
  plot(num_clusters, type = "l", col = "blue", lwd = 2,
       xlab = "Iteration", ylab = "Number of Clusters",
       main = "Cluster Count Over Iterations")
}


plot_cluster_heatmap <- function(y, cluster_assignments) {
  y_ordered <- y[order(cluster_assignments), ]
  image(t(y_ordered), col = heat.colors(100), axes = FALSE,
        main = "Monoproduct Counts Ordered by Cluster")
  box()
}

plot_confusion_matrix <- function(true_C, estimated_C) {
  tab <- table(True = true_C, Estimated = estimated_C)
  
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    message("⚠️ Only one class in true or estimated labels — using image() instead.")
    image(tab, axes = TRUE, col = cm.colors(100), main = "Confusion Matrix (Fallback)")
    return(invisible(NULL))
  }
  
  heatmap(tab, Rowv = NA, Colv = NA, col = cm.colors(100),
          main = "Confusion Matrix: True vs Estimated Clusters")
}


plot_cluster_sizes <- function(C) {
  cluster_sizes <- table(C)
  barplot(cluster_sizes, col = "skyblue", main = "Final Cluster Sizes",
          xlab = "Cluster", ylab = "Number of Agencies")
}

# Assuming you've stored cluster_history from the sampler loop
plot_cluster_evolution(cluster_history)

# Final cluster assignments
plot_cluster_heatmap(y, result$C)
plot_confusion_matrix(simulated_data$z, result$C)
plot_cluster_sizes(result$C)

