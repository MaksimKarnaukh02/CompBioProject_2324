install.packages("TreeDist")
# Load required packages
library(TreeDist)
library(ape)

# Step 1: Read and prepare the input trees
tree_files <- list.files("/data/antwerpen/208/vsc20886/trees", pattern = "*.nwk", full.names = TRUE)
trees <- lapply(tree_files, read.tree)

# Step 2: Compute tree distances
dist_matrix <- treedist(trees, method = "robinson-foulds")

# Step 3: Create a tree space plot
tree_space <- cmdscale(dist_matrix, k = 2)  # Perform MDS with 2 dimensions

# Step 4: Plot the tree space
plot(tree_space, type = "n", xlab = "Dimension 1", ylab = "Dimension 2")
text(tree_space, labels = tree_files)

# Optional: Customize the plot appearance
# Add additional plot elements, change colors, labels, etc., based on your preferences

# Save the plot
png("treespace_plot.png")
plot(tree_space, type = "n", xlab = "Dimension 1", ylab = "Dimension 2")
text(tree_space, labels = tree_files)
dev.off()