# Clear the environment
rm(list = ls())

# Load required libraries
library(tidyverse)
library(readxl)
library(patchwork)

# Set working directory
setwd("C:\\Users\\xxx\\Desktop\\genefamily")

# Load the data from an Excel file
df <- read_xlsx("GeneFamilyAnalysis.xlsx", sheet = 1, col_names = TRUE) %>%
  pivot_longer(cols = c("Core", "Pan"), names_to = "group", values_to = "value")

# Separate the data into core and pan gene families
coredf <- df %>% filter(group == "Core")
pandf <- df %>% filter(group == "Pan")

# Create the scatter plot with smooth lines
p1 <- ggplot(df, aes(x = `Individual Number`, y = value, color = group)) +
  geom_point(shape = 4, size = 6) +
  scale_x_continuous(breaks = seq(2, 28, by = 1)) +
  scale_y_continuous(labels = scales::number_format(), breaks = scales::pretty_breaks(n = 20)) +
  geom_smooth(data = coredf, method = "loess", mapping = aes(x = `Individual Number`, y = value),
              color = "black", se = FALSE) +
  geom_smooth(data = pandf, method = "loess", mapping = aes(x = `Individual Number`, y = value),
              color = "black", se = FALSE) +
  scale_color_manual(values = c("#D17431", "#39729C")) +
  labs(x = "Number of Genomes", y = "Number of Gene Families", color = "") +
  theme_classic()

# Display the plot
p1
