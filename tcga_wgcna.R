library(GSVA)
library(tidyverse)
library(WGCNA)
library(ggvenn)

set.seed(123)

tcga_expression <- readRDS("input/TCGA-GBM_expression_TPM.rds")
tcga_clinical <- readRDS("input/TCGA-GBM_clinical_survival.rds")

genecards_mes_genes <- read.csv(
  "input/signatures/GeneCards_MES_signature.csv",
  check.names = FALSE
)$gene_symbol

genecards_mes_genes <- unique(
  toupper(na.omit(genecards_mes_genes))
)

common_samples <- intersect(
  colnames(tcga_expression),
  rownames(tcga_clinical)
)

tcga_expression <- tcga_expression[
  ,
  common_samples,
  drop = FALSE
]

tcga_clinical <- tcga_clinical[
  common_samples,
  ,
  drop = FALSE
]

mes_gene_set <- list(
  MES_transition = intersect(
    genecards_mes_genes,
    rownames(tcga_expression)
  )
)

mes_ssgsea_parameter <- ssgseaParam(
  as.matrix(tcga_expression),
  mes_gene_set
)

mes_score_matrix <- gsva(mes_ssgsea_parameter)

tcga_mes_metadata <- data.frame(
  sample = colnames(tcga_expression),
  mes_score = as.numeric(
    mes_score_matrix[
      "MES_transition",
      colnames(tcga_expression)
    ]
  ),
  row.names = colnames(tcga_expression)
)

tcga_mes_metadata$mes_group <- ifelse(
  tcga_mes_metadata$mes_score >=
    median(tcga_mes_metadata$mes_score, na.rm = TRUE),
  "High",
  "Low"
)

tcga_mes_metadata$mes_group <- factor(
  tcga_mes_metadata$mes_group,
  levels = c("Low", "High")
)

gene_variance <- apply(
  tcga_expression,
  1,
  var,
  na.rm = TRUE
)

variable_genes <- names(
  gene_variance[
    gene_variance >=
      quantile(gene_variance, 0.92, na.rm = TRUE)
  ]
)

wgcna_expression <- as.data.frame(
  t(
    tcga_expression[
      variable_genes,
      ,
      drop = FALSE
    ]
  )
)

quality_check <- goodSamplesGenes(
  wgcna_expression,
  verbose = 0
)

if (!quality_check$allOK) {
  wgcna_expression <- wgcna_expression[
    quality_check$goodSamples,
    quality_check$goodGenes,
    drop = FALSE
  ]
}

tcga_mes_metadata <- tcga_mes_metadata[
  rownames(wgcna_expression),
  ,
  drop = FALSE
]

mes_traits <- model.matrix(
  ~ 0 + mes_group,
  data = tcga_mes_metadata
)

colnames(mes_traits) <- levels(
  tcga_mes_metadata$mes_group
)

rownames(mes_traits) <- rownames(
  tcga_mes_metadata
)

enableWGCNAThreads()

soft_threshold_results <- pickSoftThreshold(
  wgcna_expression,
  powerVector = 1:20,
  RsquaredCut = 0.85,
  verbose = 0
)

soft_power <- 10

adjacency_matrix <- adjacency(
  wgcna_expression,
  power = soft_power
)

tom_matrix <- TOMsimilarity(
  adjacency_matrix
)

tom_dissimilarity <- 1 - tom_matrix

gene_tree <- hclust(
  as.dist(tom_dissimilarity),
  method = "average"
)

dynamic_modules <- cutreeDynamic(
  dendro = gene_tree,
  distM = tom_dissimilarity,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = 50,
  verbose = 0
)

dynamic_module_colors <- labels2colors(
  dynamic_modules
)

merge_results <- mergeCloseModules(
  wgcna_expression,
  dynamic_module_colors,
  cutHeight = 0.4,
  verbose = 0
)

merged_module_colors <- merge_results$colors
names(merged_module_colors) <- colnames(wgcna_expression)

merged_module_eigengenes <- orderMEs(
  merge_results$newMEs
)

module_trait_correlations <- WGCNA::cor(
  merged_module_eigengenes,
  mes_traits,
  method = "spearman",
  use = "p"
)

module_trait_pvalues <- corPvalueStudent(
  module_trait_correlations,
  nrow(wgcna_expression)
)

module_trait_labels <- paste0(
  round(module_trait_correlations, 2),
  "\n(",
  signif(module_trait_pvalues, 3),
  ")"
)

dim(module_trait_labels) <- dim(
  module_trait_correlations
)

par(mfrow = c(1, 2))

plot(
  soft_threshold_results$fitIndices$Power,
  -sign(
    soft_threshold_results$fitIndices$slope
  ) *
    soft_threshold_results$fitIndices$SFT.R.sq,
  xlab = "Soft threshold power",
  ylab = "Scale-free topology model fit, signed R²",
  type = "n",
  main = "Scale independence"
)

text(
  soft_threshold_results$fitIndices$Power,
  -sign(
    soft_threshold_results$fitIndices$slope
  ) *
    soft_threshold_results$fitIndices$SFT.R.sq,
  labels = soft_threshold_results$fitIndices$Power
)

abline(
  h = 0.85,
  lty = 2
)

plot(
  soft_threshold_results$fitIndices$Power,
  soft_threshold_results$fitIndices$mean.k.,
  xlab = "Soft threshold power",
  ylab = "Mean connectivity",
  type = "n",
  main = "Mean connectivity"
)

text(
  soft_threshold_results$fitIndices$Power,
  soft_threshold_results$fitIndices$mean.k.,
  labels = soft_threshold_results$fitIndices$Power
)

par(mfrow = c(1, 1))

plotDendroAndColors(
  gene_tree,
  merged_module_colors,
  groupLabels = "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)

labeledHeatmap(
  Matrix = module_trait_correlations,
  xLabels = colnames(mes_traits),
  yLabels = colnames(merged_module_eigengenes),
  ySymbols = colnames(merged_module_eigengenes),
  textMatrix = module_trait_labels,
  colors = colorRampPalette(
    c("#6a87b7", "white", "#e2af61")
  )(50),
  setStdMargins = FALSE,
  zlim = c(-1, 1),
  main = "Module–trait relationships"
)

turquoise_module_genes <- names(
  merged_module_colors[
    merged_module_colors == "turquoise"
  ]
)

mes_degs <- read.csv(
  "input/TCGA-GBM_MES_DEGs.csv",
  check.names = FALSE
)$gene_symbol

mes_degs <- unique(
  toupper(na.omit(mes_degs))
)

candidate_genes <- sort(
  Reduce(
    intersect,
    list(
      turquoise_module_genes,
      mes_degs,
      genecards_mes_genes
    )
  )
)

fig1e_candidate_venn <- ggvenn(
  list(
    WGCNA_module = turquoise_module_genes,
    MES_DEGs = mes_degs,
    MES_signature = genecards_mes_genes
  ),
  show_percentage = FALSE,
  fill_alpha = 0.5,
  stroke_color = NA
)

logistic_results <- map_dfr(
  candidate_genes,
  function(gene) {
    model_data <- data.frame(
      os = tcga_clinical$os,
      expression = as.numeric(
        tcga_expression[
          gene,
          rownames(tcga_clinical)
        ]
      )
    ) %>%
      drop_na()

    logistic_model <- glm(
      os ~ expression,
      data = model_data,
      family = binomial
    )

    coefficient_table <- summary(
      logistic_model
    )$coefficients

    beta <- coefficient_table[
      "expression",
      "Estimate"
    ]

    standard_error <- coefficient_table[
      "expression",
      "Std. Error"
    ]

    data.frame(
      gene = gene,
      OR = exp(beta),
      lower_95_CI = exp(
        beta - 1.96 * standard_error
      ),
      upper_95_CI = exp(
        beta + 1.96 * standard_error
      ),
      p_value = coefficient_table[
        "expression",
        "Pr(>|z|)"
      ]
    )
  }
)

significant_logistic_results <- logistic_results %>%
  filter(p_value < 0.05) %>%
  arrange(p_value)

fig1f_logistic <- significant_logistic_results %>%
  mutate(
    gene = forcats::fct_reorder(
      gene,
      OR
    )
  ) %>%
  ggplot(
    aes(
      x = OR,
      y = gene
    )
  ) +
  geom_vline(
    xintercept = 1,
    linetype = 2
  ) +
  geom_errorbar(
    aes(
      xmin = lower_95_CI,
      xmax = upper_95_CI
    ),
    orientation = "y",
    width = 0.2
  ) +
  geom_point(
    size = 2
  ) +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Odds ratio",
    y = NULL
  )