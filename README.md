title: "BioinformHer R Project"
author: "Garima Rana"
date: 29 Sept 2025
output: html_document
In this project, we simulate gene expression data for 10 patients, divided into two groups: Responders and Non-Responders.
The goal is to explore gene-level statistics, compare expression patterns across groups, and generate a visualization of patient totals.
Methods: 

Step 1: Create patient IDs and groups
patient_id <- paste0("P", sprintf("%02d", 1:10))
group <- c("Responder","Responder","Non-Responder","Responder","Non-Responder",
           "Non-Responder","Responder","Non-Responder","Responder","Non-Responder")
           
Step 2: Simulate gene expression data
G1 <- round(rnorm(10, mean = 6, sd = 2), 2)
G2 <- round(rnorm(10, mean = 7, sd = 2.5), 2)
G3 <- round(rnorm(10, mean = 5.5, sd = 1.8), 2)
G4 <- round(rnorm(10, mean = 4.8, sd = 2.2), 2)
G5 <- round(rnorm(10, mean = 8, sd = 3), 2)
G6 <- round(rnorm(10, mean = 3.5, sd = 1.2), 2)
G7 <- round(rnorm(10, mean = 6.5, sd = 2), 2)
G8 <- round(rnorm(10, mean = 5, sd = 1.5), 2)

Step 3: Build final dataframe
expr <- data.frame(patient_id = patient_id,
                   group = group,
                   G1 = G1, G2 = G2, G3 = G3, G4 = G4,
                   G5 = G5, G6 = G6, G7 = G7, G8 = G8,
                   stringsAsFactors = FALSE)

head(expr)

Results:

Data structure and summary
gene_cols <- grep("^G", names(expr), value = TRUE)
str(expr)
summary(expr[, gene_cols])

Overall gene means
overall_means <- colMeans(expr[, gene_cols])
overall_means
sort(overall_means, decreasing = TRUE)[1:2] # Top 2 genes

Group-wise comparison
group_means <- sapply(expr[, gene_cols], function(g) tapply(g, expr$group, mean))
t(group_means)

# Genes higher in Responders
names(which(group_means["Responder", ] > group_means["Non-Responder", ]))
Patient totals and high expression
expr$total <- rowSums(expr[, gene_cols])
expr$HighExpr <- ifelse(expr$total > median(expr$total), "HighExpr", "LowExpr")

table(expr$group, expr$HighExpr)
expr[which.max(expr$total), c("patient_id","group","total")]

Visualization
barplot(expr$total, names.arg = expr$patient_id, las = 2,
        main = "Per-patient total expression", ylab = "Total expression")
        
Discussion
Some genes have overall higher expression than others.
Certain genes show higher means in Responders compared to Non-Responders, suggesting potential biomarkers.
The HighExpr/LowExpr classification indicates variation in total expression burden among patients.

Conclusion
This analysis demonstrates how simulated gene expression data can be explored to identify patterns across patient groups. The pipeline includes:
Data simulation
Summary statistics
Group comparisons
Visualization
