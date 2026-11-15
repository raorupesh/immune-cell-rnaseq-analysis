# Numbers log

Every number in the README comes from this file, and every entry here comes from a real run of the notebooks. The source column names the notebook and section that printed it. All runs reproduced the same values after deleting `data/` and rerunning from scratch.

| Metric | Value | Source (notebook, section) | Date |
|---|---|---|---|
| Cells before QC | 2,700 | 01, QC | 2026-11-15 |
| Cells after QC | 2,638 | 01, QC | 2026-11-15 |
| Genes after filtering | 13,714 | 01, QC | 2026-11-15 |
| Highly variable genes | 1,838 | 01, HVG | 2026-11-15 |
| Leiden clusters at resolution 0.9 | 8 | 01, Clustering | 2026-11-15 |
| Cell types labeled | 8 (CD4 T 1,152; CD14 mono 479; B 341; CD8 T 302; NK 158; FCGR3A mono 157; DC 36; platelets 13) | 01, Labels | 2026-11-15 |
| Genes kept after low-count filter | 16,139 | 02, Filtering | 2026-11-15 |
| DE genes, padj < 0.05, `~ cell + dex` | 4,081 | 02, Results | 2026-11-15 |
| Up / down | 2,190 / 1,891 | 02, Results | 2026-11-15 |
| Genes with raw p < 0.05 | 5,488 (about 807 expected by chance if no gene had an effect) | 02, Rigor R1 | 2026-11-15 |
| DE genes, padj < 0.05, `~ dex` only | 2,773 | 02, Rigor R2 | 2026-11-15 |
| Permutation: median DE genes over null labelings | 0 (14 null labelings) | 02, Rigor R3 | 2026-11-15 |
| Permutation: max DE genes over null labelings | 20 | 02, Rigor R3 | 2026-11-15 |
| FKBP5 / DUSP1 / PER1 / TSC22D3 / CRISPLD2 log2FC (apeglm) and padj | 3.97 (7.3e-26) / 2.94 (1.2e-124) / 3.17 (1.8e-81) / 3.11 (3.2e-19) / 2.60 (1.6e-45) | 02, Sanity check | 2026-11-15 |
| GO BP terms enriched (padj < 0.05, \|LFC\| > 1 genes) | 306 | 02, GO enrichment | 2026-11-15 |
| Logistic regression CV accuracy (mean ± sd) | 0.960 ± 0.009 (balanced 0.943 ± 0.011) | 03, Classifiers | 2026-11-15 |
| Random forest CV accuracy (mean ± sd) | 0.959 ± 0.009 (balanced 0.915 ± 0.041) | 03, Classifiers | 2026-11-15 |
| Resolutions tested | 0.3, 0.5, 0.7, 0.9, 1.2, 1.5, 2.0 (5 to 23 clusters) | 03, Resolution sweep | 2026-11-15 |
| ARI range across seeds | 0.949 to 1.000 (10 seeds; 9 give 8 clusters) | 03, Seeds | 2026-11-15 |
| ARI vs reference with 10 / 20 / 30 / 40 PCs | 0.541 / 0.668 / 0.884 / 1.000 | 03, Number of PCs | 2026-11-15 |
