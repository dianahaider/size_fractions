# Load necessary libraries
#install.packages("ape")

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(devtools)
library(ape)
library(vegan)
library(pairwiseAdonis)
library(phyloseq)
install.packages("remotes")
remotes::install_github("Russel88/MicEco")
library(MicEco)

# The root directory of the analysis

directory = "~/Documents/escuela/phd/size_fractions/outputs/02-EUKs/"
infile = "distance_matrix_d1.tsv"
outfile = "PCoA-plt60m16s.png"


# Load the distance matrix
# Assuming the distance matrix is saved in a file named "distance_matrix.tsv"
matrix_data <- read.table(paste(directory,infile,sep = ""), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Convert the matrix data to a 'dist' object
dist_matrix <- as.dist(matrix_data)

# Perform a Principal Coordinates Analysis (PCoA)
pcoa_result <- pcoa(dist_matrix)

# Compute variance percentages
total_variance <- sum(pcoa_result$values$Eigenvalues)
var1 <- round((pcoa_result$values$Eigenvalues[1] / total_variance) * 100, 2)
var2 <- round((pcoa_result$values$Eigenvalues[2] / total_variance) * 100, 2)

xvalues <- seq(from = floor(min(pcoa_result$vectors[,1])),
               to = ceiling(max(pcoa_result$vectors[,1])),
               by = 0.2)

yvalues <- seq(from = floor(min(pcoa_result$vectors[,2])),
               to = ceiling(max(pcoa_result$vectors[,2])),
               by = 0.2)

png(filename = paste(directory,outfile,".png",sep=""))

old_mar <- par("mar")
par(mar = c(old_mar[1], 6, old_mar[3], old_mar[4]))  # Increase the left margin

# Plot the first two PCoA axes with modified settings
plot(pcoa_result$vectors[,1], pcoa_result$vectors[,2], 
     xlab = paste("PCo1 (", var1, "%)", sep = ""),
     ylab = paste("PCo2 (", var2, "%)", sep = ""),
     type = "n", # Set up the plotting area
     cex.lab = 1.8,
     cex.axis = 1.4,
     xaxp = c(floor(min(pcoa_result$vectors[,1])), ceiling(max(pcoa_result$vectors[,1])), 
              (ceiling(max(pcoa_result$vectors[,1])) - floor(min(pcoa_result$vectors[,1]))) / 0.2),
     yaxp = c(floor(min(pcoa_result$vectors[,2])), ceiling(max(pcoa_result$vectors[,2])),
              (ceiling(max(pcoa_result$vectors[,2])) - floor(min(pcoa_result$vectors[,2]))) / 0.2))
# Add gridlines
abline(h = yvalues, col = "gray", lty = "dotted")
abline(v = xvalues, col = "gray", lty = "dotted")

type_colors = c("#4477AA","#EE7733","#BBBBBB", "Black")

# Add points based on the ending character of labels
colors <- ifelse(grepl("AS$", rownames(pcoa_result$vectors)), type_colors[1],
                 ifelse(grepl("AL$", rownames(pcoa_result$vectors)), type_colors[2],
                        ifelse(grepl("ASL", rownames(pcoa_result$vectors)), type_colors[3], "Black"))) # "ALack" as a default, just in case

points(pcoa_result$vectors[,1], pcoa_result$vectors[,2], pch = 21, bg = colors, cex = 1.6)

# Add a legend with circle markers
legend("topright", 
       legend = c("Small", "Large", "DeFr", "Whole"), 
       col = type_colors, 
       pt.bg = type_colors,  # fill color
       pch = 21, 
       cex = 1.6)



dev.off()

# Let's do the PERMANOVA
dist_matrix_rownames <- labels(dist_matrix)
group_labels_time <- ifelse(grepl("BB22.1", dist_matrix_rownames), "Pre-bloom",
                       ifelse(grepl("BB22.2", dist_matrix_rownames), "Pre-bloom",
                              ifelse(grepl("BB22.3", dist_matrix_rownames), "Pre-bloom",
                                     ifelse(grepl("BB22.4", dist_matrix_rownames), "Pre-bloom",
                                            ifelse(grepl("BB22.5", dist_matrix_rownames), "Pre-bloom",
                                                   ifelse(grepl("BB22.6", dist_matrix_rownames), "Pre-bloom",
                                                          ifelse(grepl("BB22.7", dist_matrix_rownames), "Pre-bloom", "ALoom")))))))

group_labels <- ifelse(grepl("AS$", dist_matrix_rownames), "Small",
                       ifelse(grepl("AL", dist_matrix_rownames), "Large",
                              ifelse(grepl("ASL", dist_matrix_rownames), "DeFr", "Whole")))

#permanova_result <- adonis2(as.matrix(dist_matrix) ~ group_labels, permutations = 999)

permanova_result <- adonis2(as.matrix(dist_matrix) ~ group_labels + group_labels_time, 
                            permutations = 999)

print(permanova_result)

adonis_OmegaSq(permanova_result, partial = TRUE)


# Make sure the labels are in a data frame
df_group_labels <- data.frame(group_labels)

# Run pairwise tests

posthoc_results <- pairwise.adonis(dist_matrix, df_group_labels$group_label, p.adjust.m = 'fdr')

print(posthoc_results)

