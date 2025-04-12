library(data.table)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]       # Input file name
outfile <- args[2]    # Output prefix
threads <- as.numeric(args[3])  # Number of threads for fread 
thres <- as.numeric(args[4])      # P-value significance threshold. i.e: 0.05

# Load only the 12th column (p_wald) from the input file using multiple threads
kmers <- fread(
  infile, 
  nThread = threads, 
  select = 12,            # Load only the 12th column
  col.names = "p_wald"    # Assign column name
)

# Ensure no missing or invalid values in p_wald
kmers <- kmers[!is.na(p_wald) & p_wald > 0 & p_wald <= 1]

# Number of markers (rows)
n <- nrow(kmers)

# Compute Benjamini-Hochberg (BH) adjusted P-values
kmers[, p_wald_bh := p.adjust(p_wald, method = "BH")]

# Calculate -log10(BH cutoff)
bhCutoff <- -log10(max(kmers[p_wald_bh < thres, p_wald], na.rm = TRUE))

# Calculate -log10(Bonferroni cutoff)
bfCutoff <- -log10(thres / n)

# Print results
cat("-log10(BH cutoff):", bhCutoff, "\n")
cat("-log10(BF cutoff):", bfCutoff, "\n")

# Optional: Save results to output file (if needed)
outfile <- paste0(outfile, "_results.txt")
fwrite(
  data.table(BH_Cutoff = bhCutoff, BF_Cutoff = bfCutoff), 
  file = outfile, 
  sep = "\t", 
  col.names = TRUE
)
