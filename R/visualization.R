#' Process and Visualize Weighted Microbiome Network
#'
#' This function processes a weighted microbiome network and visualizes it across different thresholds. Each threshold represents a minimum edge weight required for inclusion in the network plot.
#'
#' @param weighted_network A numeric matrix representing the weighted network where rows and columns represent taxa.
#' @param taxa_names A character vector of taxa names corresponding to rows and columns of \code{weighted_network}.
#' @param thresholds A numeric vector specifying thresholds for edge weight inclusion in network visualization.
#' @param show_labels A logical vector indicating whether to display node labels for each threshold. Recycled to match \code{thresholds} length if needed.
#' @param node_colors A character vector specifying colors for nodes at each threshold. Recycled to match \code{thresholds} length if needed.
#' @param edge_colors A character vector specifying colors for edges at each threshold. Recycled to match \code{thresholds} length if needed.
#'
#' @return A data frame containing the edges of the network, with columns for taxa pairs and edge weights.
#'
#' @details The function performs the following steps:
#'   \enumerate{
#'     \item Extracts edges from \code{weighted_network} that meet each specified threshold.
#'     \item Creates and visualizes a static graph layout for each threshold, displaying nodes and edges.
#'     \item Exports each network visualization to a JPEG file.
#'     \item Provides an interactive visualization using \code{visNetwork} for further exploration.
#'   }
#'
#' @import igraph
#' @importFrom graphics par mtext
#' @importFrom grDevices dev.off jpeg
#' @importFrom visNetwork visNetwork visEdges visOptions visLayout
#' @export
process_and_visualize_network <- function(weighted_network, taxa_names, thresholds, show_labels, node_colors, edge_colors) {
  # Ensure parameter arrays are the same length as thresholds or adjust accordingly
  show_labels <- rep(show_labels, length.out = length(thresholds))
  node_colors <- rep(node_colors, length.out = length(thresholds))
  edge_colors <- rep(edge_colors, length.out = length(thresholds))

  # Initialize edge list
  edges <- data.frame(First_Taxa = character(), Second_Taxa = character(), Weight = numeric(), stringsAsFactors = FALSE)

  # Extract upper triangle values (excluding diagonal)
  num_features <- ncol(weighted_network)
  for (i in 1:(num_features - 1)) {
    for (j in (i + 1):num_features) {
      if (weighted_network[i, j] > 0) {
        edges <- rbind(edges, data.frame(First_Taxa = taxa_names[i], Second_Taxa = taxa_names[j], Weight = weighted_network[i, j]))
      }
    }
  }

  # Sort edges by weight
  edges <- edges[order(-edges$Weight), ]

  # Save the edge list as CSV
  write.csv(edges, "edge_list.csv", row.names = FALSE)

  # Set up plot layout dynamically based on the number of thresholds
  plot_rows <- ceiling(sqrt(length(thresholds)))
  plot_cols <- ceiling(length(thresholds) / plot_rows)
  par(mfrow = c(plot_rows, plot_cols), mar = c(4, 4, 4, 4), oma = c(2, 2, 2, 2))

  # Plot each threshold
  for (i in 1:length(thresholds)) {
    threshold <- thresholds[i]
    edges_filtered <- subset(edges, Weight > threshold)

    # Create graph object for visualization
    g <- graph_from_data_frame(edges_filtered, directed = FALSE)

    # Filter nodes that have at least one edge with score > threshold
    nodes_in_edges <- unique(c(edges_filtered$First_Taxa, edges_filtered$Second_Taxa))
    g <- induced_subgraph(g, vids = V(g)[name %in% nodes_in_edges])

    # Set edge weights
    # # E(g)$weight <- edges_filtered$Weight

    # Calculate graph properties
    num_nodes <- vcount(g)
    num_edges <- ecount(g)
    max_degree <- max(degree(g))
    max_degree_nodes <- V(g)[degree(g) == max_degree]$name  # Ensure unique node labels for max degree nodes
    max_degree_nodes_str <- paste(unique(max_degree_nodes), collapse = ", ")

    # Set vertex labels based on user input
    vertex_label <- if (show_labels[i]) V(g)$name else NA

    # Adjust node size based on degree
    vertex_size <- 4

    # Plot the graph with force-directed layout, adjusted node size, and customized edges and nodes
    set.seed(123)  # Fix the layout for reproducibility
    plot(g, layout = layout_with_fr(g), vertex.label = vertex_label, vertex.size = vertex_size, vertex.color = node_colors[i], edge.width = 2, edge.color = edge_colors[i],
         vertex.label.cex = 0.7,          # Use a moderate value, like 1.5 or 2
         vertex.label.family = "Arial",   # Set font family if needed
         main = paste("Microbiome Network with Score > ", threshold), cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
    mtext(side = 1, line = 4, paste("Nodes: ", num_nodes, " Edges: ", num_edges, " Max-Degree: ", max_degree))
  }
  par(mfrow = c(1, 1))  # Reset to default layout

  # Save plot as JPEG with high resolution
  jpeg(filename = paste0("network_", paste(thresholds, collapse = ","), ".jpeg"), width = 3200, height = 2400, res = 300)
  par(mfrow = c(plot_rows, plot_cols), mar = c(4, 4, 4, 4), oma = c(2, 2, 2, 2))
  for (i in 1:length(thresholds)) {
    threshold <- thresholds[i]
    edges_filtered <- subset(edges, Weight > threshold)

    # Create graph object for visualization
    g <- graph_from_data_frame(edges_filtered, directed = FALSE)

    # Filter nodes that have at least one edge with score > threshold
    nodes_in_edges <- unique(c(edges_filtered$First_Taxa, edges_filtered$Second_Taxa))
    g <- induced_subgraph(g, vids = V(g)[name %in% nodes_in_edges])

    # Set edge weights
    E(g)$weight <- edges_filtered$Weight

    # Set vertex labels based on user input
    vertex_label <- if (show_labels[i]) V(g)$name else NA

    # Adjust node size based on degree
    vertex_size <- 4

    # Calculate graph properties for each subplot
    num_nodes <- vcount(g)
    num_edges <- ecount(g)
    max_degree <- max(degree(g))
    max_degree_nodes <- V(g)[degree(g) == max_degree]$name  # Ensure unique node labels for max degree nodes
    max_degree_nodes_str <- paste(unique(max_degree_nodes), collapse = ", ")

    # Plot the graph
    set.seed(123)  # Fix the layout for reproducibility
    plot(g, layout = layout_with_fr(g), vertex.label = vertex_label, vertex.size = vertex_size, vertex.color = node_colors[i],
         vertex.label.cex = 0.7,          # Use a moderate value, like 1.5 or 2
         vertex.label.family = "Arial",   # Set font family if needed
         edge.width = 2, edge.color = edge_colors[i],
         main = paste("Microbiome Network with Score > ", threshold), cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)
    mtext(side = 1, line = 4, paste("Nodes: ", num_nodes, " Edges: ", num_edges, " Max-Degree: ", max_degree))
  }
  par(mfrow = c(1, 1))  # Reset to default layout
  dev.off()

  # Interactive Visualization
  for (i in 1:length(thresholds)) {
    threshold <- thresholds[i]
    edges_filtered <- subset(edges, Weight > threshold)

    # Create nodes and edges data frames for visNetwork
    nodes <- data.frame(id = unique(c(edges_filtered$First_Taxa, edges_filtered$Second_Taxa)),
                        label = ifelse(show_labels[i], as.character(unique(c(edges_filtered$First_Taxa, edges_filtered$Second_Taxa))), NA),
                        value = degree(g)[match(unique(c(edges_filtered$First_Taxa, edges_filtered$Second_Taxa)), V(g)$name)],
                        stringsAsFactors = FALSE)

    edges_vis <- data.frame(from = edges_filtered$First_Taxa,
                            to = edges_filtered$Second_Taxa,
                            value = 1,
                            title = "Edge",
                            stringsAsFactors = FALSE)

    # Create interactive network plot
    visNetwork(nodes, edges_vis, main = paste("Interactive Network with Score >", threshold)) %>%
      visEdges(scaling = list(min = 1, max = 10)) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visLayout(randomSeed = 123)
  }

  # Return edges data frame for further use
  return(edges)
}
