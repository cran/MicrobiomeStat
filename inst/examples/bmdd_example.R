# BMDD (Bimodal Dirichlet Distribution) Example
# Demonstrates the use of BMDD method for zero-inflated microbiome data analysis
# integrated within the MicrobiomeStat package

library(MicrobiomeStat)
library(phyloseq)  # BiocManager::install("phyloseq") if needed

# ============================================================================
# Example 1: Basic BMDD Usage with Simulated Data
# ============================================================================

cat("\n=== Example 1: Basic BMDD Usage ===\n")

# Simulate microbiome count data
set.seed(123)
m <- 50  # number of taxa
n <- 30  # number of samples
W <- matrix(rpois(m * n, lambda = 100), nrow = m, ncol = n)
rownames(W) <- paste0('Taxon', 1:m)
colnames(W) <- paste0('Sample', 1:n)

# Fit bimodal Dirichlet distribution
cat("\nFitting BMDD model...\n")
bmdd.obj <- bmdd(W = W, type = "count", trace = TRUE)

# View results structure
cat("\nBMDD results structure:\n")
str(bmdd.obj)

# Key outputs:
cat("\nKey outputs:\n")
cat("- gamma: bimodality indicators (", length(bmdd.obj$gamma), " taxa)\n")
cat("- pi: mixing proportions (", length(bmdd.obj$pi), " samples)\n")
cat("- beta: posterior Dirichlet parameters (", dim(bmdd.obj$beta)[1], "x", dim(bmdd.obj$beta)[2], ")\n")
cat("- alpha0: mode 0 parameters (", length(bmdd.obj$alpha$alp0), " taxa)\n")
cat("- alpha1: mode 1 parameters (", length(bmdd.obj$alpha$alp1), " taxa)\n")
cat("- Method used:", bmdd.obj$method, "\n")

# ============================================================================
# Example 2: Extracting Posterior Mean Compositions
# ============================================================================

cat("\n=== Example 2: Posterior Mean Compositions ===\n")

# Extract posterior parameters
beta <- bmdd.obj$beta

# Calculate posterior mean compositions
post.mean <- t(t(beta) / colSums(beta))

# Verify compositions sum to 1
all_sum_to_one <- all(abs(colSums(post.mean) - 1) < 1e-10)
cat("All compositions sum to 1:", all_sum_to_one, "\n")

# Display first few taxa for first 3 samples
cat("\nPosterior mean compositions (first 5 taxa, first 3 samples):\n")
print(post.mean[1:5, 1:3])

# ============================================================================
# Example 3: Generating Posterior Samples for Downstream Analysis
# ============================================================================

cat("\n=== Example 3: Generating Posterior Samples ===\n")

# Function to handle zero values in generated samples
zero.fun <- function(X) {
  X <- t(apply(X, 1, function(x) {
    if(all(x == 0)) {
      x[x == 0] <- min(X[X != 0])
    } else {
      x[x == 0] <- min(x[x != 0])
    }
    return(x)
  }))
  return(X)
}

# Generate K posterior samples per sample
K <- 100
cat("Generating", K, "posterior samples per observation...\n")

beta_rep <- beta[, rep(1:n, K)]
X <- matrix(rgamma(m * n * K, beta_rep, 1), m)
X <- t(t(X) / colSums(X))

# Handle any zeros
if(any(X == 0)) {
  cat("Handling zero values in generated samples...\n")
  X <- zero.fun(X)
  X <- t(t(X) / colSums(X))
}

# Set names
colnames(X) <- paste0('Sample', rep(1:n, K), '_Rep', rep(1:K, each = n))
rownames(X) <- rownames(beta)

cat("Generated posterior samples dimensions:", dim(X), "\n")
cat("These samples can be used for differential abundance testing\n")

# ============================================================================
# Example 4: Real Microbiome Data with LinDA Integration
# ============================================================================

cat("\n=== Example 4: Real Data Analysis with LinDA ===\n")

# Load example dataset (included in package)
data(phy)
cat("Loaded example phyloseq object\n")

# Filter function: remove low-prevalence taxa and low-depth samples
otu_filter <- function(feature.dat, prev = 0.1, dep = 1000){
  idx <- apply(feature.dat, 1, function(x) sum(x > 0) > (ncol(feature.dat) * prev))
  idx2 <- colSums(feature.dat) > dep
  return(feature.dat[idx, idx2])
}

# Extract and filter data
feature.dat <- as.data.frame(as.matrix(phyloseq::otu_table(phy)))
meta.dat <- as.data.frame(as.matrix(phyloseq::sample_data(phy)))
meta.dat$grp <- as.factor(meta.dat$grp)

cat("Original data dimensions:", dim(feature.dat), "\n")

feature.dat <- otu_filter(feature.dat)
meta.dat <- meta.dat[colnames(feature.dat), ]

m <- nrow(feature.dat)
n <- ncol(feature.dat)

cat("Filtered data dimensions:", m, "taxa x", n, "samples\n")

# Fit BMDD model
cat("\nFitting BMDD model to real data...\n")
bmdd.obj <- bmdd(W = feature.dat, type = 'count', trace = TRUE)

# Posterior mean
beta <- bmdd.obj$beta
post.mean <- t(t(beta) / colSums(beta))

cat("\nGenerating 100 posterior samples for LinDA analysis...\n")

# Generate 100 posterior samples per observation
K <- 100
beta_exp <- beta[, rep(1:n, K)]
X <- matrix(rgamma(m * n * K, beta_exp, 1), m)
X <- t(t(X) / colSums(X))

if(any(X == 0)) {
  X <- zero.fun(X)
  X <- t(t(X) / colSums(X))
}

colnames(X) <- paste0('sample', 1:(n * K))
rownames(X) <- rownames(beta)

# ============================================================================
# Differential abundance analysis with LinDA
# ============================================================================

cat("\nPerforming differential abundance analysis with LinDA...\n")

# Apply to BMDD-imputed proportions with 100 replicates per sample
id <- factor(rep(1:n, K))
grp <- rep(meta.dat$grp, K)
Z <- cbind.data.frame(grp, id)
rownames(Z) <- colnames(X)

cat("\nLinDA with BMDD-imputed data (accounting for imputation uncertainty)...\n")
linda.bmdd.obj <- linda(feature.dat = X, meta.dat = Z,
                        formula = "~ grp + (1 | id)",
                        feature.dat.type = "proportion")

# Apply to original count matrix (for comparison)
cat("\nLinDA with original count data...\n")
linda.obj <- linda(feature.dat = feature.dat, meta.dat = meta.dat,
                   formula = "~ grp",
                   feature.dat.type = "count")

# Display results
cat("\n=== Results Summary ===\n")
cat("\nTop 10 differentially abundant taxa (BMDD + LinDA):\n")
top_taxa_bmdd <- head(linda.bmdd.obj$output[[1]][order(linda.bmdd.obj$output[[1]]$padj), ], 10)
print(top_taxa_bmdd[, c("baseMean", "log2FoldChange", "pvalue", "padj")])

cat("\nTop 10 differentially abundant taxa (Original + LinDA):\n")
top_taxa_orig <- head(linda.obj$output[[1]][order(linda.obj$output[[1]]$padj), ], 10)
print(top_taxa_orig[, c("baseMean", "log2FoldChange", "pvalue", "padj")])

# ============================================================================
# Example 5: Method Selection and Performance
# ============================================================================

cat("\n=== Example 5: Method Selection ===\n")

# Auto-select method (uses NLopt if available)
cat("\nAuto method selection:\n")
result_auto <- bmdd(W, type = "count", method = "auto", trace = FALSE)
cat("Method used:", result_auto$method, "\n")

# Force R implementation
cat("\nForcing R implementation:\n")
result_r <- bmdd(W, type = "count", method = "R", trace = FALSE)
cat("Method used:", result_r$method, "\n")

# Try NLopt implementation
cat("\nTrying NLopt implementation (requires NLopt library):\n")
result_nlopt <- tryCatch({
  bmdd(W, type = "count", method = "nlopt", trace = FALSE)
}, error = function(e) {
  cat("NLopt not available:", e$message, "\n")
  NULL
})

if (!is.null(result_nlopt)) {
  cat("Method used:", result_nlopt$method, "\n")
  cat("Note: NLopt provides 50-90x speedup for large datasets\n")
}

cat("\n=== BMDD Integration Complete ===\n")
cat("\nFor more information:\n")
cat("- ?bmdd for function documentation\n")
cat("- ?bmdd.nlopt for NLopt implementation details\n")
cat("- data(phy) for example dataset\n")
cat("\nReference:\n")
cat("Zhou, H., Chen, J., & Zhang, X. (2025). BMDD: A Probabilistic Framework\n")
cat("for Accurate Imputation of Zero-inflated Microbiome Sequencing Data.\n")
cat("PLOS Computational Biology.\n")
cat("https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013124\n")
