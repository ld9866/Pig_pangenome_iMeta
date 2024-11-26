# Clear workspace to start fresh
rm(list = ls())

# Load necessary libraries
library(data.table)  # For reading large datasets efficiently
library(psych)       # For visualizing data relationships
library(vegan)       # For performing RDA (Redundancy Analysis)

# Load genetic data (final raw file with genetic information)
gen <- fread("final.raw", header = TRUE, check.names = FALSE, data.table = FALSE)

# Clean the genetic data by removing unnecessary columns (first 6 columns)
gen_clean <- gen[, c(-1:-6)]
gen_mat <- as.matrix(gen_clean)  # Convert data to matrix form
rownames(gen_mat) <- gen[, 1]   # Set row names from the first column of the gen data

# Check dimensions and missing values
dim(gen_mat)
sum(is.na(gen_mat))

# Impute missing values using the most frequent value in each column (here, it's commented out, as imputation isn't needed)
gen.imp <- apply(gen_mat, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

# Confirm no missing values after imputation
sum(is.na(gen.imp))

# Load environmental data (csv file)
env <- read.csv("env.csv")

# Convert the individual column to character format for consistency
env$individual <- as.character(env$individual)

# Check if the row names in genetic data match the individual IDs in environmental data
identical(rownames(gen_mat), env[,1])

# Visualize the relationship between environmental variables (bio variables)
library(png)  # For saving plots as images
file_path <- "1.png"
png(file_path, width = 8200, height = 8200, units = "px", res = 300)
pairs.panels(env[, 5:35], scale = TRUE)  # Scatter plot matrix of selected environmental variables
dev.off()

# Subset environmental data to select relevant predictors for RDA
pred <- subset(env, select = -c(bio1, bio4, bio5, bio6, bio8, bio9, bio13, bio14, bio15, bio16, bio17, bio18, soc15cm, tmin2021.2040, prec2021.2040))

# Further subset for a focused set of variables
pred1 <- pred[, 5:20]
colnames(pred1) <- c("bio2", "bio3", "bio7", "bio10", "bio11", "bio12", "bio19", 
                     "Sum_of_UV-B_Radiation_of_Highest_Quarter", "elev_units", "bdod15cm", 
                     "cfvo15cm", "clay15cm", "nitrogen60cm", "ocd60cm", "phd15cm", "silt15cm")

# Visualize relationships between selected environmental predictors
pdf("2.pdf")
pairs.panels(pred1, scale = TRUE)
dev.off()

# Perform Redundancy Analysis (RDA) using environmental predictors
rda_532model <- rda(gen_mat ~ ., data = pred1, scale = TRUE)

# Check RDA results and perform further analysis
RsquareAdj(rda_532model)  # Adjusted R^2 for the model
summary(eigenvals(rda_532model, model = "constrained"))

# Create a scree plot for the RDA model
pdf("screeplot.pdf")
screeplot(rda_532model)
dev.off()

# Check for multicollinearity using VIF (Variance Inflation Factor)
vif.cca(rda_532model)

# Create a color mapping for the ecotype variable and plot the RDA results
levels(env$ecotype) <- c("ED", "AD", "EW", "AW")
eco <- env$ecotype
bg <- c("#ff7f00", "#1f78b4", "#ffff33", "#33a02c")
bg <- rep(bg, length.out = length(eco))

# Plot RDA with points representing species and sites
pdf("rda1.pdf")
plot(rda_532model, type = "n", scaling = 3)  # Empty plot, setting up scaling
points(rda_532model, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3)
points(rda_532model, display = "sites", pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = bg)
text(rda_532model, scaling = 3, display = "bp", col = "#0868ac", cex = 1)
legend("bottomright", legend = levels(eco), bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bg)
dev.off()

# Create an alternative RDA plot with specific axes (1st and 3rd)
pdf("rda2.pdf")
plot(rda_532model, type = "n", scaling = 3, choices = c(1, 3))
points(rda_532model, display = "species", pch = 20, cex = 0.7, col = "gray32", scaling = 3, choices = c(1, 3))
points(rda_532model, display = "sites", pch = 21, cex = 1.3, col = "gray32", scaling = 3, bg = bg, choices = c(1, 3))
text(rda_532model, scaling = 3, display = "bp", col = "#0868ac", cex = 1, choices = c(1, 3))
legend("bottomright", legend = levels(eco), bty = "n", col = "gray32", pch = 21, cex = 1, pt.bg = bg)
dev.off()

# Extract and visualize loadings from RDA
load.rda <- scores(rda_532model, choices = c(1:3), display = "species")

# Plot histograms of loadings for RDA components 1, 2, and 3
pdf("load.rda.pdf")
hist(load.rda[, 1], main = "Loadings on RDA1")
hist(load.rda[, 2], main = "Loadings on RDA2")
hist(load.rda[, 3], main = "Loadings on RDA3")
dev.off()

# Detect outliers based on z-scores
outliers <- function(x, z) {
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}

# Identify outliers in RDA loadings
cand1 <- outliers(load.rda[, 1], 3)
cand2 <- outliers(load.rda[, 2], 3)
cand3 <- outliers(load.rda[, 3], 3)
ncand <- length(cand1) + length(cand2) + length(cand3)

# Combine and label outliers
cand1 <- cbind.data.frame(rep(1, times = length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2, times = length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3, times = length(cand3)), names(cand3), unname(cand3))

# Create a combined data frame of outliers
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# Create a matrix to store correlations between outlier SNPs and environmental variables
foo <- matrix(nrow = ncand, ncol = 16)
colnames(foo) <- c("bio2", "bio3", "bio7", "bio10", "bio11", "bio12", "bio19", 
                   "Sum_of_UV-B_Radiation_of_Highest_Quarter", "elev_units", "bdod15cm", 
                   "cfvo15cm", "clay15cm", "nitrogen60cm", "ocd60cm", "phd15cm", "silt15cm")

# Correlate each SNP with the environmental variables and add to the 'foo' matrix
for (i in 1:length(cand$snp)) {
  nam <- cand[i, 2]
  snp.gen <- gen_mat[, nam]
  foo[i, ] <- apply(pred1, 2, function(x) cor(x, snp.gen))
}

# Add correlations to the outlier data frame
cand <- cbind.data.frame(cand, foo)

# Filter out duplicated SNPs
cand <- cand[!duplicated(cand$snp),]

# Determine the strongest predictor for each outlier SNP
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i, 20] <- names(which.max(abs(bar[4:19])))  # strongest predictor
  cand[i, 21] <- max(abs(bar[4:19]))  # correlation value
}

# Rename the new columns for clarity
colnames(cand)[20] <- "predictor"
colnames(cand)[21] <- "correlation"

# Output the results to a CSV file
write.csv(cand, file = "cand.csv")

# Ensure that the columns listed for "keep" and "remove" are correctly identified in the dataset.
# Keep the specified variables
# Remove the specified variables