# bioinform project demo script
#Description:
# This script simulates patient-level gene expression data,
# performs exploratory analysis, and generates a barplot of totals.
# Set seed for reproducibility (same random numbers each run
#Create patient IDs and groups
set.seed(123)
patient_id <- paste0("P", sprintf("%02d", 1:10))
group <- c("Responder","Responder","Non-Responder","Responder","Non-Responder",
           "Non-Responder","Responder","Non-Responder","Responder","Non-Responder")

# 2. Simulate genes (10 patients)
G1 <- round(rnorm(10, mean = 6, sd = 2), 2)
G2 <- round(rnorm(10, mean = 7, sd = 2.5), 2)
G3 <- round(rnorm(10, mean = 5.5, sd = 1.8), 2)
G4 <- round(rnorm(10, mean = 4.8, sd = 2.2), 2)
G5 <- round(rnorm(10, mean = 8, sd = 3), 2)
G6 <- round(rnorm(10, mean = 3.5, sd = 1.2), 2)
G7 <- round(rnorm(10, mean = 6.5, sd = 2), 2)
G8 <- round(rnorm(10, mean = 5, sd = 1.5), 2)

# 3. Build expr and verify it exists
expr <- data.frame(patient_id = patient_id,
                   group = group,
                   G1 = G1, G2 = G2, G3 = G3, G4 = G4,
                   G5 = G5, G6 = G6, G7 = G7, G8 = G8,
                   stringsAsFactors = FALSE)

if (!exists("expr") || !is.data.frame(expr)) stop("ERROR: 'expr' not created correctly.")
cat("Created expr: dimensions =", dim(expr), "\n")
print(head(expr))

# 4. Simple diagnostics
gene_cols <- grep("^G", names(expr), value = TRUE)
cat("Gene columns detected:", paste(gene_cols, collapse = ", "), "\n\n")

# 5. Structure & summary
str(expr)
print(summary(expr[ , gene_cols]))

# 6. Per-gene overall means and top 2
overall_means <- colMeans(expr[ , gene_cols])
cat("\nOverall gene means:\n"); print(overall_means)
cat("\nTop 2 genes by mean:\n"); print(sort(overall_means, decreasing = TRUE)[1:2])

# 7. Responder vs Non-Responder means
group_means <- sapply(expr[ , gene_cols], function(g) tapply(g, expr$group, mean))
cat("\nGroup means (Responder / Non-Responder):\n")
print(t(group_means))
cat("\nGenes higher in Responders:\n")
print(names(which(group_means["Responder", ] > group_means["Non-Responder", ])))

# 8. Totals and HighExpr label
expr$total <- rowSums(expr[ , gene_cols])
expr$HighExpr <- ifelse(expr$total > median(expr$total), "HighExpr", "LowExpr")
cat("\nCounts HighExpr by group:\n"); print(table(expr$group, expr$HighExpr))
cat("\nRow with highest total:\n"); print(expr[which.max(expr$total), c("patient_id","group","total")])

# 9. Save a barplot to disk (avoids GUI/display problems)
png(filename = "patient_totals.png", width = 900, height = 400)
barplot(expr$total, names.arg = expr$patient_id, las = 2,
        main = "Per-patient total expression", ylab = "Total expression")
dev.off()
cat("\nSaved barplot as: patient_totals.png (in working directory)\n")


