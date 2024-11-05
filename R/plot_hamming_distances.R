#' Calculate and Plot Hamming Distances Between Binary Network Matrices
#'
#' This function calculates the Hamming distance, common edges, and number of edges
#' for each pair of binary network matrices in a specified folder, then visualizes
#' the top N pairs with the lowest Hamming distances.
#'
#' @param binary_network_folder A character string specifying the path to the folder containing binary network CSV files.
#' @param top_n_pairs Integer specifying the number of pairs with the lowest Hamming distances to visualize.
#' @param output_filename Character string specifying the name of the output JPEG file for the plot.
#' @return A data frame containing the calculated distances and metrics for each pair.
#' @importFrom ggplot2 ggplot aes geom_bar geom_text coord_flip labs theme_minimal ggsave
#' @importFrom stats reorder
#' @importFrom readr read_csv cols col_double
#' @export
plot_hamming_distances <- function(binary_network_folder, top_n_pairs, output_filename = "hamming_distances_plot.jpeg") {

  # Get a list of all CSV files in the folder
  file_list <- list.files(binary_network_folder, pattern = "*.csv", full.names = TRUE)

  # Create an empty data frame to store pairwise Hamming distances and additional metrics
  distance_results <- data.frame(Pair = character(), Hamming_Distance = numeric(), Common_Edges = numeric(), Edges_Net1 = numeric(), Edges_Net2 = numeric(), stringsAsFactors = FALSE)

  # Iterate over each pair of files and calculate metrics
  for (i in 1:(length(file_list) - 1)) {
    for (j in (i + 1):length(file_list)) {
      net1 <- as.matrix(suppressMessages(read_csv(file_list[i], col_types = cols(.default = col_double()))))
      net2 <- as.matrix(suppressMessages(read_csv(file_list[j], col_types = cols(.default = col_double()))))

      net1 <- net1[, -1]
      net2 <- net2[, -1]

      # Ensure both matrices have the same dimensions
      if (!is.null(net1) && !is.null(net2) && !all(dim(net1) == dim(net2))) {
        stop("Matrices have different dimensions: ", file_list[i], " and ", file_list[j])
      }

      # Calculate Hamming distance
      hamming_dist <- sum(net1 != net2)

      # Calculate the number of common edges
      common_edges_count <- sum(net1 == 1 & net2 == 1)

      # Calculate the number of edges in each network
      edges_net1 <- sum(net1 == 1)
      edges_net2 <- sum(net2 == 1)

      # Store results
      pair_name <- paste(basename(file_list[i]), "vs", basename(file_list[j]))
      distance_results <- rbind(distance_results, data.frame(Pair = pair_name, Hamming_Distance = hamming_dist, Common_Edges = common_edges_count, Edges_Net1 = edges_net1, Edges_Net2 = edges_net2))
    }
  }

  # Clean up pair names
  distance_results$Pair <- sub("_binary.csv$", "", distance_results$Pair)
  distance_results$Pair <- sub("_binary.csv", "", distance_results$Pair)
  distance_results$Pair <- sub("spiecEasi_", "SE_", distance_results$Pair)
  distance_results$Pair <- sub("vs spiecEasi_", "vs SE_", distance_results$Pair)

  # Select the top pairs with the lowest Hamming distances for visualization
  top_pairs <- distance_results[order(distance_results$Hamming_Distance), ][1:top_n_pairs, ]

  # Plot the Hamming distances with the number of common edges and edges in each network
  plot <- ggplot(top_pairs, aes(x = reorder(Pair, Hamming_Distance), y = Hamming_Distance)) +
    geom_bar(stat = "identity", fill = "#11A0D9") +
    geom_text(aes(label = paste("Common:", Common_Edges, "\nNet1:", Edges_Net1, "\nNet2:", Edges_Net2)), hjust = -0.1, size = 2) +
    coord_flip() +
    labs(title = "Top Pairwise Hamming Distances Between Networks",
         x = "",
         y = "Hamming Distance") +
    theme_minimal()

  # Save the plot as a JPEG file
  ggsave(output_filename, plot = plot, width = 13, height = 6, units = "in")

  # Return the distance results data frame
  return(distance_results)
}
