library(data.table)
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]       # Input file name
outfile <- args[2]      # Output prefix
threads <- as.numeric(args[3])  # Number of threads for fread 
thres <- as.numeric(args[4])    # Base P-value threshold (e.g., 0.05)
k_size <- as.numeric(args[5])   # k-mer size (e.g., 31)

# Load only the 12th column (p_wald) from the input file using multiple threads
kmers <- fread(
  infile, 
  nThread = threads, 
  select = 12,            
  col.names = "p_wald"    
)

# Ensure no missing or invalid values in p_wald
kmers <- kmers[!is.na(p_wald) & p_wald > 0 & p_wald <= 1]

# Total number of k-mers
n_total <- nrow(kmers)

# Calculate effective number of tests based on k-mer size
n_effective <- n_total / k_size

# Calculate adjusted thresholds
bh_cutoff <- thres / n_effective  # More relaxed threshold
bf_cutoff <- thres / n_total      # Traditional Bonferroni

# Compute Benjamini-Hochberg (BH) adjusted P-values
kmers[, p_wald_bh := p.adjust(p_wald, method = "BH")]

# Calculate cutoffs for Manhattan plots
log_bh_relaxed <- -log10(bh_cutoff)
log_bf_standard <- -log10(bf_cutoff)
log_sig_kmers <- -log10(max(kmers[p_wald_bh < thres, p_wald], na.rm = TRUE))

# Print results
cat("Total k-mers:", n_total, "\n")
cat("Effective number of tests:", n_effective, "\n")
cat("Relaxed threshold (α/effective tests):", bh_cutoff, "\n")
cat("Standard Bonferroni threshold (α/total tests):", bf_cutoff, "\n")
cat("-log10(Relaxed threshold):", log_bh_relaxed, "\n")
cat("-log10(Bonferroni threshold):", log_bf_standard, "\n")
cat("-log10(BH significant k-mers threshold):", log_sig_kmers, "\n")

# Save results
fwrite(
  data.table(
    Total_Kmers = n_total,
    Effective_Tests = n_effective,
    Relaxed_Threshold = bh_cutoff,
    Bonferroni_Threshold = bf_cutoff,
    Log10_Relaxed = log_bh_relaxed,
    Log10_Bonferroni = log_bf_standard,
    Log10_Sig_BH = log_sig_kmers
  ), 
  file = paste0(outfile, "_thresholds.txt"), 
  sep = "\t", 
  col.names = TRUE
)
