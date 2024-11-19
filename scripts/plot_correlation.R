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
  sample_name <- gsub("_bowtie2.fragmentsCount.bin500.bed", "", basename(file_path))
  fragData <- read.table(file_path, header = FALSE)
  colnames(fragData) <- c("chrom", "bin", sample_name)

  if (is.null(fragCount)) {
    fragCount <- fragData
  } else {
    fragCount <- full_join(fragCount, fragData, by = c("chrom", "bin"))
  }
}

# Replace zero or negative values with 1 and filter zero-variance columns
filtered_fragCount <- fragCount %>%
  mutate(across(-c(chrom, bin), ~ replace(.x, .x <= 0, 1))) %>%
  select(-chrom, -bin) %>%
  select(where(~ var(.x, na.rm = TRUE) > 0))

# Stop if not enough valid data for correlation
if (ncol(filtered_fragCount) < 2) {
  stop("Not enough valid columns for correlation calculation after filtering.")
}

# Calculate correlation matrix
M <- cor(filtered_fragCount %>% log2(), use = "complete.obs")

# Adjust clustering dynamically
addrect_value <- ifelse(ncol(M) >= 3, 3, ncol(M))

# Generate correlation plot
output_file <- file.path(output_dir, "fragCount_correlation_plot.pdf")
pdf(output_file, width = 8, height = 8)
corrplot(M, method = "color", outline = TRUE, addgrid.col = "darkgray", order = "hclust",
         addrect = addrect_value, rect.col = "black", rect.lwd = 3, cl.pos = "b", tl.col = "indianred4",
         tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1,
         col = colorRampPalette(c("midnightblue", "white", "darkred"))(100))
dev.off()
