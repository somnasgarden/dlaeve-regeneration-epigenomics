###############################################################################
# PART 3: Hierarchical Clustering
###############################################################################

load("Part2_DEGs/data/Part2_objects.RData")

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(cluster)

# ── Output directories ───────────────────────────────────────────────────────
dir.create("Part3_Hierarchical/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("Part3_Hierarchical/data",  recursive = TRUE, showWarnings = FALSE)

# ── 1. Hierarchical clustering on GENES ──────────────────────────────────────
gene_dist <- as.dist(1 - cor(t(deg_scaled), method = "pearson"))
gene_hclust <- hclust(gene_dist, method = "complete")

# ── 2. Hierarchical clustering on SAMPLES ────────────────────────────────────
sample_dist <- dist(t(deg_scaled), method = "euclidean")
sample_hclust <- hclust(sample_dist, method = "complete")

# ── 3. Full heatmap ──────────────────────────────────────────────────────────
pheatmap(deg_scaled,
         cluster_rows = gene_hclust, cluster_cols = sample_hclust,
         annotation_col = annotation_col,
         show_rownames = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Hierarchical Clustering: All DEGs",
         fontsize = 8, fontsize_col = 6,
         filename = "Part3_Hierarchical/plots/hierarchical_full.png",
         width = 14, height = 10)
pheatmap(deg_scaled,
         cluster_rows = gene_hclust, cluster_cols = sample_hclust,
         annotation_col = annotation_col,
         show_rownames = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Hierarchical Clustering: All DEGs",
         fontsize = 8, fontsize_col = 6,
         filename = "Part3_Hierarchical/plots/hierarchical_full.pdf",
         width = 14, height = 10)

# ── 4. Silhouette to pick k ─────────────────────────────────────────────────
sil_scores <- sapply(2:10, function(k) {
  clusters <- cutree(gene_hclust, k = k)
  sil <- silhouette(clusters, gene_dist)
  mean(sil[, 3])
})

sil_df <- data.frame(K = 2:10, Silhouette = sil_scores)
optimal_k_hc <- sil_df$K[which.max(sil_df$Silhouette)]
cat("Optimal k (silhouette):", optimal_k_hc, "\n\n")

p_sil <- ggplot(sil_df, aes(x = K, y = Silhouette)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_point(data = sil_df[sil_df$K == optimal_k_hc, ], color = "red", size = 5) +
  annotate("text", x = optimal_k_hc, y = max(sil_scores) + 0.01,
           label = paste("Optimal k =", optimal_k_hc), color = "red") +
  labs(title = "Silhouette: Hierarchical Clustering", x = "k", y = "Avg Silhouette Width") +
  scale_x_continuous(breaks = 2:10) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave("Part3_Hierarchical/plots/silhouette_hierarchical.png", p_sil, width = 8, height = 6, dpi = 300)
ggsave("Part3_Hierarchical/plots/silhouette_hierarchical.pdf", p_sil, width = 8, height = 6)

# ── 5. Cut dendrogram ────────────────────────────────────────────────────────
n_clusters <- optimal_k_hc   # ← CHANGE if needed

gene_clusters_hc <- cutree(gene_hclust, k = n_clusters)
cat("Cluster sizes:\n")
print(table(gene_clusters_hc))

# ── 6. Gapped heatmap ───────────────────────────────────────────────────────
gene_order <- order(gene_clusters_hc)
deg_scaled_ordered <- deg_scaled[gene_order, ]
clusters_ordered <- gene_clusters_hc[gene_order]

cluster_annotation_row <- data.frame(
  Cluster = factor(clusters_ordered),
  row.names = rownames(deg_scaled_ordered)
)
gap_positions <- which(diff(clusters_ordered) != 0)

pheatmap(deg_scaled_ordered,
         cluster_rows = FALSE, cluster_cols = sample_hclust,
         annotation_row = cluster_annotation_row, annotation_col = annotation_col,
         gaps_row = gap_positions,
         show_rownames = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = paste0("Hierarchical Clusters (k=", n_clusters, ")"),
         fontsize = 8, fontsize_col = 6,
         filename = "Part3_Hierarchical/plots/hierarchical_gapped.png",
         width = 14, height = 10)
pheatmap(deg_scaled_ordered,
         cluster_rows = FALSE, cluster_cols = sample_hclust,
         annotation_row = cluster_annotation_row, annotation_col = annotation_col,
         gaps_row = gap_positions,
         show_rownames = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = paste0("Hierarchical Clusters (k=", n_clusters, ")"),
         fontsize = 8, fontsize_col = 6,
         filename = "Part3_Hierarchical/plots/hierarchical_gapped.pdf",
         width = 14, height = 10)

# ── 7. Centroids ─────────────────────────────────────────────────────────────
cluster_centroids_hc <- do.call(rbind, lapply(1:n_clusters, function(k) {
  genes_in_k <- names(gene_clusters_hc[gene_clusters_hc == k])
  colMeans(deg_scaled[genes_in_k, , drop = FALSE])
}))
rownames(cluster_centroids_hc) <- paste0("Cluster_", 1:n_clusters)

pheatmap(cluster_centroids_hc,
         cluster_rows = FALSE, cluster_cols = sample_hclust,
         annotation_col = annotation_col,
         show_rownames = TRUE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Hierarchical Cluster Centroids",
         fontsize = 10, fontsize_col = 7,
         filename = "Part3_Hierarchical/plots/hierarchical_centroids.png",
         width = 14, height = 6)
pheatmap(cluster_centroids_hc,
         cluster_rows = FALSE, cluster_cols = sample_hclust,
         annotation_col = annotation_col,
         show_rownames = TRUE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Hierarchical Cluster Centroids",
         fontsize = 10, fontsize_col = 7,
         filename = "Part3_Hierarchical/plots/hierarchical_centroids.pdf",
         width = 14, height = 6)

# ── Save data ────────────────────────────────────────────────────────────────
hc_assignments <- data.frame(Gene = names(gene_clusters_hc), Cluster = gene_clusters_hc)
write.csv(hc_assignments, "Part3_Hierarchical/data/hierarchical_cluster_assignments.csv", row.names = FALSE)
write.csv(sil_df, "Part3_Hierarchical/data/silhouette_scores.csv", row.names = FALSE)

save(dds, vsd, normalized_counts, df,
     deg_union, deg_expression, deg_scaled, annotation_col, all_degs,
     gene_hclust, sample_hclust, gene_dist,
     gene_clusters_hc, n_clusters, cluster_centroids_hc,
     file = "Part3_Hierarchical/data/Part3_objects.RData")

cat("\nCheck Part3_Hierarchical/plots/\n")
cat("═══ Part 3 complete ═══\n")