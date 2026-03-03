###############################################################################
# PART 4: K-means Clustering & Comparison
###############################################################################

load("Part3_Hierarchical/data/Part3_objects.RData")

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(cluster)
library(mclust)
library(reshape2)

# ── Output directories ───────────────────────────────────────────────────────
dir.create("Part4_Kmeans/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("Part4_Kmeans/data",  recursive = TRUE, showWarnings = FALSE)

# ── 1. Elbow Method ──────────────────────────────────────────────────────────
cat("Elbow analysis (K=1 to 15)...\n")
set.seed(42)
max_k <- 15

wss <- sapply(1:max_k, function(k) {
  kmeans(deg_scaled, centers = k, nstart = 25, iter.max = 100)$tot.withinss
})

elbow_df <- data.frame(K = 1:max_k, WSS = wss)

p_elbow <- ggplot(elbow_df, aes(x = K, y = WSS)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 3) +
  labs(title = "Elbow Method", x = "K", y = "Within-Cluster SS") +
  scale_x_continuous(breaks = 1:max_k) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("Part4_Kmeans/plots/elbow_plot.png", p_elbow, width = 8, height = 6, dpi = 300)
ggsave("Part4_Kmeans/plots/elbow_plot.pdf", p_elbow, width = 8, height = 6)

# ── 2. Silhouette ────────────────────────────────────────────────────────────
cat("Silhouette analysis (K=2 to 15)...\n")
set.seed(42)

sil_scores_km <- sapply(2:max_k, function(k) {
  km <- kmeans(deg_scaled, centers = k, nstart = 25, iter.max = 100)
  ss <- silhouette(km$cluster, dist(deg_scaled))
  mean(ss[, 3])
})

sil_df_km <- data.frame(K = 2:max_k, Silhouette = sil_scores_km)
optimal_k_km <- sil_df_km$K[which.max(sil_df_km$Silhouette)]
cat("Optimal K (silhouette):", optimal_k_km, "\n\n")

p_sil <- ggplot(sil_df_km, aes(x = K, y = Silhouette)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 3) +
  geom_point(data = sil_df_km[sil_df_km$K == optimal_k_km, ], color = "red", size = 5) +
  annotate("text", x = optimal_k_km, y = max(sil_scores_km) + 0.01,
           label = paste("Optimal K =", optimal_k_km), color = "red") +
  labs(title = "Silhouette: K-means", x = "K", y = "Avg Silhouette Width") +
  scale_x_continuous(breaks = 2:max_k) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("Part4_Kmeans/plots/silhouette_kmeans.png", p_sil, width = 8, height = 6, dpi = 300)
ggsave("Part4_Kmeans/plots/silhouette_kmeans.pdf", p_sil, width = 8, height = 6)

# ── 3. Apply K-means ─────────────────────────────────────────────────────────
k_final <- optimal_k_km   # ← CHANGE if needed
cat("K-means with K =", k_final, "...\n\n")

set.seed(42)
km_result <- kmeans(deg_scaled, centers = k_final, nstart = 25, iter.max = 100)
gene_clusters_km <- km_result$cluster

cat("K-means cluster sizes:\n")
print(table(gene_clusters_km))

# ── 4. Gapped heatmap ───────────────────────────────────────────────────────
km_order <- order(gene_clusters_km)
deg_scaled_ordered_km <- deg_scaled[km_order, ]
clusters_ordered_km <- gene_clusters_km[km_order]

cluster_annotation_km <- data.frame(
  Cluster = factor(clusters_ordered_km),
  row.names = rownames(deg_scaled_ordered_km)
)
gap_positions_km <- which(diff(clusters_ordered_km) != 0)

pheatmap(deg_scaled_ordered_km,
         cluster_rows = FALSE, cluster_cols = sample_hclust,
         annotation_row = cluster_annotation_km, annotation_col = annotation_col,
         gaps_row = gap_positions_km,
         show_rownames = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = paste0("K-means Clustering (K=", k_final, ")"),
         fontsize = 8, fontsize_col = 6,
         filename = "Part4_Kmeans/plots/kmeans_gapped.png",
         width = 14, height = 10)
pheatmap(deg_scaled_ordered_km,
         cluster_rows = FALSE, cluster_cols = sample_hclust,
         annotation_row = cluster_annotation_km, annotation_col = annotation_col,
         gaps_row = gap_positions_km,
         show_rownames = FALSE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = paste0("K-means Clustering (K=", k_final, ")"),
         fontsize = 8, fontsize_col = 6,
         filename = "Part4_Kmeans/plots/kmeans_gapped.pdf",
         width = 14, height = 10)

# ── 5. Centroids ─────────────────────────────────────────────────────────────
cluster_centroids_km <- km_result$centers

pheatmap(cluster_centroids_km,
         cluster_rows = FALSE, cluster_cols = sample_hclust,
         annotation_col = annotation_col,
         show_rownames = TRUE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "K-means Cluster Centroids",
         fontsize = 10, fontsize_col = 7,
         filename = "Part4_Kmeans/plots/kmeans_centroids.png",
         width = 14, height = 6)
pheatmap(cluster_centroids_km,
         cluster_rows = FALSE, cluster_cols = sample_hclust,
         annotation_col = annotation_col,
         show_rownames = TRUE, show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "K-means Cluster Centroids",
         fontsize = 10, fontsize_col = 7,
         filename = "Part4_Kmeans/plots/kmeans_centroids.pdf",
         width = 14, height = 6)

# ── 6. Cluster sizes ─────────────────────────────────────────────────────────
size_df <- data.frame(Cluster = factor(1:k_final), Count = as.integer(table(gene_clusters_km)))

p_size <- ggplot(size_df, aes(x = Cluster, y = Count, fill = Cluster)) +
  geom_col() + geom_text(aes(label = Count), vjust = -0.5) +
  labs(title = "K-means Cluster Sizes", x = "Cluster", y = "Genes") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none")
ggsave("Part4_Kmeans/plots/kmeans_cluster_sizes.png", p_size, width = 8, height = 6, dpi = 300)
ggsave("Part4_Kmeans/plots/kmeans_cluster_sizes.pdf", p_size, width = 8, height = 6)

# ── 7. Compare Hierarchical vs K-means ───────────────────────────────────────
cat("\n══ Hierarchical vs K-means ══\n")

if (n_clusters == k_final) {
  ari <- adjustedRandIndex(gene_clusters_hc, gene_clusters_km)
  cat(sprintf("ARI: %.3f\n", ari))
  
  confusion <- table(Hierarchical = gene_clusters_hc, Kmeans = gene_clusters_km)
  print(confusion)
  
  conf_melt <- melt(confusion)
  conf_melt$Hierarchical <- as.factor(conf_melt$Hierarchical)
  conf_melt$Kmeans <- as.factor(conf_melt$Kmeans)
  
  p_compare <- ggplot(conf_melt, aes(x = Kmeans, y = Hierarchical, fill = value)) +
    geom_tile() +
    geom_text(aes(label = value), color = "white", size = 5) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(title = sprintf("Hierarchical vs K-means (ARI=%.3f)", ari),
         x = "K-means Cluster", y = "Hierarchical Cluster", fill = "Genes") +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave("Part4_Kmeans/plots/method_comparison.png", p_compare, width = 8, height = 6, dpi = 300)
  ggsave("Part4_Kmeans/plots/method_comparison.pdf", p_compare, width = 8, height = 6)
} else {
  cat("Different k values — skipping ARI. Compare heatmaps visually.\n")
}

# ── 8. Save data ─────────────────────────────────────────────────────────────
km_assignments <- data.frame(Gene = names(gene_clusters_km), Kmeans_Cluster = gene_clusters_km)
write.csv(km_assignments, "Part4_Kmeans/data/kmeans_cluster_assignments.csv", row.names = FALSE)

combined <- data.frame(
  Gene         = names(gene_clusters_km),
  Hierarchical = gene_clusters_hc[names(gene_clusters_km)],
  Kmeans       = gene_clusters_km
)
write.csv(combined, "Part4_Kmeans/data/combined_cluster_assignments.csv", row.names = FALSE)

for (k in 1:k_final) {
  genes_k <- names(gene_clusters_km[gene_clusters_km == k])
  writeLines(genes_k, paste0("Part4_Kmeans/data/kmeans_cluster_", k, "_genes.txt"))
  cat(sprintf("Cluster %d: %d genes\n", k, length(genes_k)))
}

write.csv(elbow_df, "Part4_Kmeans/data/elbow_scores.csv", row.names = FALSE)
write.csv(sil_df_km, "Part4_Kmeans/data/silhouette_scores.csv", row.names = FALSE)

save(dds, vsd, normalized_counts, df,
     deg_union, deg_expression, deg_scaled, annotation_col, all_degs,
     gene_hclust, sample_hclust, gene_dist,
     gene_clusters_hc, n_clusters, cluster_centroids_hc,
     km_result, gene_clusters_km, k_final, cluster_centroids_km,
     file = "Part4_Kmeans/data/Part4_objects.RData")

cat("\n═══ Part 4 complete. All clustering done! ═══\n")