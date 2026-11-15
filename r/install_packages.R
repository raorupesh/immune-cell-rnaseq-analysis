# Run once. Takes 10 to 20 minutes the first time.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages(c("ggplot2", "ggrepel", "pheatmap", "rmarkdown", "knitr"))

BiocManager::install(c(
  "DESeq2",
  "airway",
  "apeglm",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
), ask = FALSE, update = FALSE)

# Check that every package loads
pkgs <- c("DESeq2", "airway", "apeglm", "clusterProfiler", "org.Hs.eg.db",
          "enrichplot", "ggplot2", "ggrepel", "pheatmap")
ok <- sapply(pkgs, requireNamespace, quietly = TRUE)
print(ok)
if (all(ok)) message("All packages installed.") else message("Missing: ", paste(pkgs[!ok], collapse = ", "))
