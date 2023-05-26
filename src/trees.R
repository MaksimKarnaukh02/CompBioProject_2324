# Load required packages
library(TreeDist)
library(ape)

# Step 1: Read and prepare the input trees
tree_files <- list.files("./Data/trees", pattern = "*.nwk", full.names = TRUE)
trees <- lapply(tree_files, read.tree)

treeNumbers <- c(1:411)
spectrum <- hcl.colors(411, "plasma")
treeCols <- spectrum[treeNumbers]

# Step 2: Compute tree distances
# method = "robinson-foulds"
# dist_matrix <- TreeDistance(trees)
distances <- RobinsonFoulds(trees)

mapping <- cmdscale(distances, k = 2)

par(mar = rep(0, 4))
plot(mapping,
     asp = 1, # Preserve aspect ratio - do not distort distances
     ann = FALSE, axes = FALSE, # Don't label axes: dimensions are meaningless
     col = treeCols, pch = 16
     )

# # Step 3: Create a tree space plot
# tree_space <- cmdscale(dist_matrix, k = 2)  # Perform MDS with 2 dimensions
#
# # Step 4: Plot the tree space
# plot(tree_space, type = "n", xlab = "Dimension 1", ylab = "Dimension 2")
# text(tree_space, labels = tree_files)
#
# # Optional: Customize the plot appearance
# # Add additional plot elements, change colors, labels, etc., based on your preferences
#
# # Save the plot
# png("treespace_plot.png")
# plot(tree_space, type = "n", xlab = "Dimension 1", ylab = "Dimension 2")
# text(tree_space, labels = tree_files)
# dev.off()

