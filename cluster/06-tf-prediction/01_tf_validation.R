# ============================================================
# Analyze DeepTFactor predictions vs functional annotation
# ============================================================

library(tidyverse)

# --- File paths (adjust as needed) ---
fasta_file <- "/mnt/data/alfredvar/rlopezt/DeepFactor1/DeepFactorV1/deeptfactor/result/derLaeGenome_namesDlasi_v2.fasta.functional_note.proteins.fasta"
pred_file  <- "/mnt/data/alfredvar/rlopezt/DeepFactor1/DeepFactorV1/deeptfactor/result/prediction_result.txt"

# -------------------------------------------------------
# 1. Parse the FASTA headers to get ID and annotation
# -------------------------------------------------------
fasta_lines <- readLines(fasta_file)
header_lines <- fasta_lines[grepl("^>", fasta_lines)]

# Extract sequence ID and whether function is unknown
fasta_df <- tibble(header = header_lines) %>%
  mutate(
    sequence_ID = str_extract(header, "(?<=^>)\\S+"),
    annotation  = ifelse(grepl('function unknown', header, ignore.case = TRUE),
                         "unknown", "annotated")
  )

cat("===== PROTEOME SUMMARY =====\n")
cat("Total proteins in FASTA:", nrow(fasta_df), "\n")
print(table(fasta_df$annotation))
cat("\n")

# -------------------------------------------------------
# 2. Read DeepTFactor predictions
# -------------------------------------------------------
pred_df <- read_tsv(pred_file, show_col_types = FALSE)

cat("===== DEEPTFACTOR SUMMARY =====\n")
cat("Total predictions:", nrow(pred_df), "\n")
cat("TF (True):", sum(pred_df$prediction == TRUE),  "\n")
cat("Non-TF (False):", sum(pred_df$prediction == FALSE), "\n\n")

# -------------------------------------------------------
# 3. Merge and build 2x2 table
# -------------------------------------------------------
merged <- inner_join(pred_df, fasta_df, by = "sequence_ID")

merged <- merged %>%
  mutate(
    TF_status = ifelse(prediction == TRUE, "TF", "non-TF"),
    annotation = factor(annotation, levels = c("annotated", "unknown"))
  )

cat("===== 2x2 CONTINGENCY TABLE =====\n")
contingency <- table(TF_status = merged$TF_status,
                     Annotation = merged$annotation)
print(contingency)
cat("\n")

# Also show proportions
cat("===== PROPORTIONS (row-wise) =====\n")
print(round(prop.table(contingency, margin = 1) * 100, 2))
cat("\n")

# Optional: Fisher's exact test
cat("===== FISHER'S EXACT TEST =====\n")
print(fisher.test(contingency))