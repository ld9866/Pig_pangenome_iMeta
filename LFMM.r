# Convert VCF to BED format using PLINK
plink --bfile 04.geno0.1maf0.05mind0.5 --recode vcf --out 04.geno0.1maf0.05mind0.5

# Impute missing genotypes using Beagle
java -jar beagle.01Mar24.d36.jar gt=04.geno0.1maf0.05mind0.5.vcf out=04.geno0.1maf0.05mind0.5.beagle

# Convert back from VCF to BED format using PLINK
plink --vcf 04.geno0.1maf0.05mind0.5.beagle.vcf.gz --make-bed --out 553

# Convert SV to SNP format and adjust .bim file
awk '{print "chr" $1 "\t" "chr"$1"_"$4 "\t" $3 "\t" $4 "\t" "A" "\t" "T"}' 553.bim > chr.bim
cp 553.bed chr.bed
cp 553.fam chr.fam

# Generate raw data file for LFMM format conversion
plink --bfile chr --recodeA --out final
cut -d " " -f 7- final.raw | sed '1d' > final.lfmm

# Copy necessary R script files to the working directory
cp fn-landgen1.R landgen1.R lib.R ./

# Read and process genotype data for LFMM analysis
x = "final1"
Y <- data.table::fread(paste(x, ".raw", sep = ""), na.strings = "NA", header = T)
Y <- Y[, -c(1:6)]
fwrite(Y, paste(x, ".lfmm", sep = ""), col.names = F, row.names = F, sep = "\t", na = "9")
Y <- LEA::lfmm2geno(paste(x, ".lfmm", sep = ""))

# Read environmental data (ecotype)
X <- read.table("eco.txt", h = T, stringsAsFactors = F)
X <- as.matrix(X[, c(5:20)])

# PCA analysis using LEA
pc <- LEA::pca(paste(x, ".lfmm", sep = ""), scale = TRUE)
tw <- LEA::tracy.widom(pc)
print(tw$pvalues[1:10])

# Scree plot of PCA eigenvalues
pdf("1.pdf")
plot(tw$eigenvalues, main = "Scree plot", ylab = "Eigenvalues", xlab = "PCs", t = "b")
dev.off()

# Perform SNP-based clustering (snmf) and plot results
obj.snmf <- LEA::snmf(Y, K = 1:20, entropy = T, ploidy = 2, CPU = 48, project = "new")
pdf("2.pdf")
plot(obj.snmf, pch = 16, col = "blue")
dev.off()

# Barplot of the inferred genetic clusters
K <- 16
pdf("3.pdf")
barplot(t(Q(obj.snmf, K = K)), col = c("forestgreen", "blue", "pink", "black"), border = NA, las = 2)
box(lwd = 2)
dev.off()

# LFMM ridge regression for association testing
mod.lfmm <- lfmm::lfmm_ridge(Y = Y, X = X, K = K)

# Perform LFMM tests for associations with environmental variables
pv <- lfmm::lfmm_test(Y = Y, X = X, lfmm = mod.lfmm, calibrate = "gif")
pvalues <- pv$calibrated.pvalue

# Save p-values and plot histograms
write.table(pvalues, file = "pvalues.txt", sep = "\t")

# Histogram of p-values for each environmental variable
pdf("4.pdf")
par(mfrow = c(3, 3), mar = c(3.5, 3.5, 3, 0.5), mgp = c(2.5, 0.8, 0))
for (i in 1:ncol(pvalues)) {
  hist(pvalues[, i], breaks = 20, xlab = "p-values", xlim = c(0, 1), main = paste("Env. variable ", i, sep = ""), col = "darkgray", border = "darkgray")
}
dev.off()

# Q-Q plot for p-values
pdf("5.pdf")
qqplot(rexp(length(pvalues), rate = log(10)), -log10(pvalues), xlab = "Expected quantile", ylab = expression(-log[10] ~ p - values), pch = 19, cex = 0.4)
abline(0, 1)
dev.off()

# Control False Discovery Rate (FDR) with q-values
qobj <- qvalue::qvalue(pvalues[, 1])
pdf("6.pdf")
hist(qobj)
dev.off()
pdf("7.pdf")
plot(qobj)
dev.off()

# Apply q-value calculation for each environmental variable
my_qvalue <- function(x) {
  q <- qvalue::qvalue(x)
  return(q$qvalues)
}

qvalues <- apply(pvalues, 2, my_qvalue)

# Q-value histograms for each environmental variable
pdf("8.pdf")
par(mfrow = c(3, 3), mar = c(3.5, 3.5, 3, 0.5), mgp = c(2.5, 0.8, 0))
for (i in 1:ncol(qvalues)) {
  hist(qvalues[, i], breaks = 20, xlab = "q-values", xlim = c(0, 1), main = paste("Env. variable ", i, sep = ""), col = "darkgray", border = "darkgray")
}
dev.off()

# Writing final results with SNP information and q-values
map <- read.table(paste("chr1", ".map", sep = ""))
qvalues <- as.data.frame(qvalues)
qvalues$SNP <- as.character(map$V2)
qvalues$CHR <- map$V1
qvalues$BP <- map$V4

write.table(qvalues, file = "qvalues.txt", sep = "\t")

# Decide on a q-value cutoff and filter results
lfmm_qvalcut <- function(qvalues, cutoff) {
  res <- list()
  e <- qvalues[, 1:(ncol(qvalues) - 3)]
  s <- qvalues[, (ncol(qvalues) - 2):ncol(qvalues)]
  
  for (i in 1:ncol(e)) {
    j <- which(e[, i] < cutoff)
    if (length(j) == 0) {
      t <- data.frame("ENV" = colnames(e)[i], "qval" = NA, "SNP" = NA, "CHR" = NA, "BP" = NA)
      res[[i]] <- t
    } else {
      t <- data.frame("ENV" = colnames(e)[i], "qval" = e[j, i], "SNP" = s[j, 1], "CHR" = s[j, 2], "BP" = s[j, 3])
      res[[i]] <- t
    }
  }
  res <- do.call(rbind.data.frame, res)
  res <- na.omit(res)
  return(res)
}

# Filter results based on q-value cutoff
lfmm.res <- lfmm_qvalcut(qvalues = qvalues, cutoff = 0.01)
write.table(lfmm.res, file = "lfmm.res.txt", sep = "\t")

# Generate Manhattan plots for results
manh <- pvalues[, c(17:19, 1:16)]
write.table(manh, file = "manh.txt", sep = "\t")

# Create Manhattan plots for each environmental variable
pdf("9.pdf")
par(mfrow = c(3, 1), oma = c(1, 1, 1, 1), mar = c(4, 5, 3, 1), mgp = c(2.5, 0.8, 0))
for (i in 4:ncol(manh)) {
  manh_i <- manh[, c(1:3, i)]
  SNPs <- lfmm.res$SNP[lfmm.res$ENV == colnames(manh_i)[4]]
  colnames(manh_i)[4] <- "P"
  qqman::manhattan(manh_i, genomewideline = F, suggestiveline = F, main = paste("LFMM - Env. variable ", colnames(manh)[i], sep = ""), highlight = SNPs, col = c(alpha("gray70", 0.7), alpha("gray30", 0.7)), xlab = "Chromosome", cex.lab = 1.5, cex.axis = 1)
  box()
}
dev.off()

# Repeat plotting for the other files as needed
