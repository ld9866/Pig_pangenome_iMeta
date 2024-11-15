# Load required R packages
library(tidyverse)       # For data manipulation and visualization
library(RColorBrewer)    # For color palettes
library(MetBrewer)       # For additional color palettes
library(ggnewscale)      # To use multiple color scales in one plot

# Set working directory (replace with your own path)
setwd("YOUR_DIRECTORY_PATH")

# Load the data from a tab-separated file
df <- read_tsv("YOUR_FILE_NAME.txt") %>%
  select(1, 3, 4, 5, 6) %>%                             # Select specific columns
  mutate(count = Input.number / Background.number) %>%   # Calculate the enrichment ratio
  arrange(desc(count)) %>%                              # Arrange data in descending order of the count
  head(27)                                              # Select the top 27 rows for plotting

# Define the order of factors for plotting
df$layer3 <- factor(df$layer3, levels = df$layer3 %>% unique() %>% rev())

# Data visualization
df %>%
  ggplot(aes(count, layer3)) +
  geom_segment(aes(x = 0, xend = count, y = layer3, yend = layer3, color = group),
               size = 6, show.legend = FALSE, alpha = 0.5) +  # Draw colored segments
  geom_text(aes(x = 0, y = layer3, label = layer3),           # Add text labels
            vjust = 0.5, hjust = 0, size = 2.8, color = "black") +
  scale_color_manual(values = c("#3B9AB2", "#00A08A")) +      # Custom color palette
  new_scale_color() +                                         # Use a new color scale
  geom_point(aes(size = Input.number, color = FDR, fill = FDR), shape = 19) + # Points representing enrichment
  scale_fill_gradientn(colors = met.brewer("Cassatt1")) +     # Gradient color scale for FDR
  scale_color_gradientn(colors = met.brewer("Cassatt1")) +    # Gradient color for points
  guides(size = guide_legend(title = "Input.number")) +       # Legend for point size
  facet_grid(group ~ ., scale = "free_y") +                   # Facet by group, free y-axis scaling
  scale_x_continuous(expand = c(0, 0.005)) +                  # Adjust x-axis scale
  theme_classic() +                                           # Classic theme
  coord_cartesian(clip = "off") +                             # Allow labels to go outside the plotting area
  theme(
    axis.title = element_blank(),                             # Remove axis titles
    axis.text.y = element_blank(),                            # Remove y-axis text
    axis.text.x = element_text(color = "black", face = "bold"), # Customize x-axis text
    axis.ticks.y = element_blank(),                           # Remove y-axis ticks
    panel.background = element_blank(),                       # Remove panel background
    strip.background = element_rect(fill = "grey85", colour = NA), # Facet background
    panel.spacing.y = unit(0, "cm"),                          # Adjust panel spacing
    legend.key = element_blank(),                             # Remove legend key background
    legend.title = element_text(color = "black", size = 9),   # Customize legend title
    legend.text = element_text(color = "black", size = 8),    # Customize legend text
    legend.spacing.x = unit(0.1, 'cm'),                       # Adjust legend spacing
    legend.key.width = unit(0.5, 'cm'),                       # Set legend key width
    legend.key.height = unit(0.5, 'cm'),                      # Set legend key height
    legend.background = element_blank(),                      # Remove legend background
    legend.box.margin = margin(1, 1, 1, 1),                   # Set legend margin
    legend.position = c(0.92, 0.78),                          # Position legend
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")            # Set plot margins
  )
