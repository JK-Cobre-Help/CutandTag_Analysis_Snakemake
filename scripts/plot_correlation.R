#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(corrplot)

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[-length(args)]  # All but last argument are input files
output_dir <- args[length(args)]    # Last argument is the output directory

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Initialize variables
fragCount <- NULL

# Loop through each input file to collect fragment count data
for (file_path in input_files) {
  
  # Extract sample name from file name
  sample_name <- gsub("_bowtie2.fragmentsCount.bin500.bed", "", basename(file_path))
  
  # Read the fragment count data
  fragCountTmp <- read.table(file_path, header = FALSE)
  colnames(fragCountTmp) <- c("chrom", "bin", sample_name)
  
  # Combine with existing data
  if (is.null(fragCount)) {
    fragCount <- fragCountTmp
  } else {
    fragCount <- full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
  }
}

# Add a pseudocount to handle zeros and log-transform the data
fragCount <- fragCount %>%
  mutate(across(-c(chrom, bin), ~ . + 1)) %>%  # Add pseudocount
  mutate(across(-c(chrom, bin), log2))        # Log2 transform

# Remove zero-variance columns
fragCount_filtered <- fragCount %>%
  select(where(~ !is.na(var(.x, na.rm = TRUE)) && var(.x, na.rm = TRUE) > 0))

# If not enough columns for correlation, generate a blank plot
if (ncol(fragCount_filtered) <= 2) {
  output_file <- file.path(output_dir, "fragCount_correlation_plot.pdf")
  pdf(output_file, width = 8, height = 8)
  plot.new()
  text(0.5, 0.5, "Insufficient data for correlation", cex = 1.5)
  dev.off()
  quit(save = "no", status = 0)
}

# Calculate correlation matrix
M <- cor(fragCount_filtered %>% select(-chrom, -bin), use = "complete.obs")

# Dynamically set addrect based on the number of samples
num_samples <- ncol(M)
addrect_value <- ifelse(num_samples >= 3, 3, num_samples)

# Generate correlation plot
output_file <- file.path(output_dir, "fragCount_correlation_plot.pdf")
pdf(output_file, width = 8, height = 8)
corrplot(M, method = "color", outline = TRUE, addgrid.col = "darkgray", order = "hclust",
         addrect = addrect_value, rect.col = "black", rect.lwd = 3, cl.pos = "b", tl.col = "indianred4",
         tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1,
         col = colorRampPalette(c("midnightblue", "white", "darkred"))(100))
dev.off()
