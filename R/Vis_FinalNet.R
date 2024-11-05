#' Plot Network in High Resolution and Interactive Format
#'
#' This function generates a network plot from a binary adjacency matrix using the \code{igraph} and \code{visNetwork} packages.
#' It saves a high-resolution static image of the network and provides an interactive visualization.
#'
#' @param network_final A binary adjacency matrix where values greater than a specified threshold are set to 1, indicating the presence of an edge.
#' @param node_color Character. Color of the nodes in the plot. Default is \code{"skyblue"}.
#' @param edge_color Character. Color of the edges in the plot. Default is \code{"grey"}.
#' @param label_color Character. Color of the node labels in the plot. Default is \code{"black"}.
#'
#' @details The function creates a static plot saved as a high-resolution PNG file and also generates an interactive network plot using \code{visNetwork}.
#' Nodes and edges are colored based on the specified parameters, and labels are positioned within the nodes.
#' The interactive plot highlights nodes and edges upon selection.
#'
#' @return A \code{visNetwork} interactive plot of the network is displayed in the R console.
#'
#' @importFrom igraph graph_from_data_frame layout_with_fr V
#' @importFrom visNetwork visNetwork visNodes visEdges visOptions visLayout
#' @export
# Function to generate and plot network using igraph and save it in high resolution
plot_network <- function(network_final, node_color = "skyblue", edge_color = "grey", label_color = "black") {

  # Generate edge list for the network
  edge_list <- which(network_final == 1, arr.ind = TRUE)
  edges <- data.frame(
    from = rownames(network_final)[edge_list[, 1]],
    to = colnames(network_final)[edge_list[, 2]]
  )

  # Create an igraph object
  network_graph <- graph_from_data_frame(d = edges, directed = FALSE)

  # Set up layout and plot parameters
  set.seed(123)  # Set seed for reproducibility of layout
  layout <- layout_with_fr(network_graph, niter = 500, grid = "nogrid") * 3  # Scale layout

  # Plotting with graphics for a high-resolution static image
  png(filename = "network_plot.png", width = 2000, height = 2000, res = 300)
  plot(network_graph, layout = layout, vertex.color = node_color, vertex.size = 5,
       vertex.label.color = label_color, vertex.label.cex = 0.6, vertex.frame.color = "black",
       edge.color = edge_color, edge.width = 2, main = "Consensus network",
       vertex.label.dist = 0,            # Set label distance to zero to center it within the node
       vertex.label.degree = 0,          # Center the label within the node
       vertex.label.font = 2,            # Optional: adjust font style
       vertex.shape = "circle")          # Use circle shape to better center labels
  dev.off()  # Close the PNG device

  # Create an interactive network plot with visNetwork
  nodes <- data.frame(id = V(network_graph)$name, label = V(network_graph)$name, color = node_color)
  edges <- data.frame(from = edges$from, to = edges$to, color = edge_color)

  visNetwork(nodes, edges, main = "Interactive Network Plot") %>%
    visNodes(color = list(border = "black", highlight = "orange"), font = list(color = label_color)) %>%
    visEdges(color = list(color = edge_color, highlight = "orange")) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visLayout(randomSeed = 123)
}
