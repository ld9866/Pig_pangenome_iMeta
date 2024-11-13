#Analysis of core variable gene families

# Loop through all .fa files in the current directory
# Run the primary_transcript.py script to extract the primary transcript for each .fa file
# Note: Adjust the path to the primary_transcript.py script and the working directory as needed

for f in *.fa; do
  python /path/to/OrthoFinder/tools/primary_transcript.py "$f"
done


# Clear the workspace to start fresh
rm(list = ls())

# Load necessary libraries
library(tidyverse)  # Data manipulation and visualization
library(readxl)     # Read Excel files
library(patchwork)  # Combine plots

# Set the working directory (update with your actual path, but keep it private)
setwd("C:/path/to/your/directory")

# Read the gene family data from an Excel file
# Replace "genefamily2.xlsx" with your file name, and make sure to hide private details
df <- read_xlsx("file_name.xlsx", sheet = 1, col_names = TRUE) %>%
  pivot_longer(cols = c("Core", "Pan"), names_to = "group", values_to = "value")

# Separate the data for core and pan gene families
coredf <- df %>% filter(group == "Core")
pandf <- df %>% filter(group == "Pan")

# Create a scatter plot with LOESS smoothing lines for core and pan gene families
p1 <- ggplot(df, aes(x = `Individual Number`, y = value, color = group)) +
  geom_point(shape = 4, size = 6) +  # Scatter plot points
  scale_x_continuous(breaks = seq(2, 28, by = 1)) +  # X-axis with breaks
  scale_y_continuous(labels = scales::number_format(), breaks = scales::pretty_breaks(n = 20)) +  # Y-axis formatting
  geom_smooth(coredf, method = "loess",  # LOESS line for core genes
              mapping = aes(x = `Individual Number`, y = value),
              color = "black",
              se = FALSE) +
  geom_smooth(pandf, method = "loess",  # LOESS line for pan genes
              mapping = aes(x = `Individual Number`, y = value),
              color = "black",
              se = FALSE) +
  scale_color_manual(values = c("#D17431", "#39729C")) +  # Custom colors for the groups
  labs(x = "Number of genomes", y = "Number of gene families", color = "") +  # Labels
  theme_classic()  # Classic theme for a clean look

# Display the plot
p1
