library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(harmony)
library(AUCell)
library(edgeR)
library(Matrix)
library(ggplot2)
library(ggrepel)

set.seed(123)

sample_directories <- list.dirs(
  "input/GSE162631",
  recursive = FALSE,
  full.names = TRUE
)

sample_ids <- sub(
  "_.*$",
  "",
  basename(sample_directories)
)

count_matrices <- setNames(
  lapply(
    sample_directories,
    function(sample_directory) {
      Read10X(
        sample_directory,
        gene.column = 2
      )
    }
  ),
  sample_ids
)

seurat_objects <- Map(
  function(count_matrix, sample_id) {
    object <- CreateSeuratObject(
      counts = count_matrix,
      project = sample_id,
      min.cells = 3,
      min.features = 200
    )

    object$sample <- sample_id
    object
  },
  count_matrices,
  sample_ids
)

sc_data <- merge(
  seurat_objects[[1]],
  y = seurat_objects[-1],
  add.cell.ids = sample_ids
)

sc_data <- JoinLayers(sc_data)

sc_data$percent_mt <- PercentageFeatureSet(
  sc_data,
  pattern = "^MT-"
)

sc_data$log10_genes_per_umi <- log10(
  sc_data$nFeature_RNA
) / log10(
  sc_data$nCount_RNA
)

sc_data <- subset(
  sc_data,
  subset =
    nCount_RNA > 200 &
    nCount_RNA < 50000 &
    nFeature_RNA > 200 &
    log10_genes_per_umi > 0.8 &
    percent_mt < 20
)

single_cell_experiment <- as.SingleCellExperiment(
  sc_data
)

single_cell_experiment <- scDblFinder(
  single_cell_experiment,
  samples = "sample"
)

sc_data$doublet_class <- single_cell_experiment$scDblFinder.class

sc_data <- subset(
  sc_data,
  subset = doublet_class == "singlet"
)

sc_data <- NormalizeData(
  sc_data,
  verbose = FALSE
)

s_genes <- CaseMatch(
  search = cc.genes$s.genes,
  match = rownames(sc_data)
)

g2m_genes <- CaseMatch(
  search = cc.genes$g2m.genes,
  match = rownames(sc_data)
)

sc_data <- CellCycleScoring(
  sc_data,
  s.features = s_genes,
  g2m.features = g2m_genes
)

sc_data <- SCTransform(
  sc_data,
  vars.to.regress = c(
    "S.Score",
    "G2M.Score",
    "percent_mt"
  ),
  verbose = FALSE
)

sc_data <- RunPCA(
  sc_data,
  assay = "SCT",
  npcs = 50,
  seed.use = 123,
  verbose = FALSE
)

sc_data <- RunHarmony(
  sc_data,
  group.by.vars = "sample",
  reduction.use = "pca",
  dims.use = 1:30,
  plot_convergence = FALSE
)

sc_data <- FindNeighbors(
  sc_data,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

sc_data <- FindClusters(
  sc_data,
  resolution = 0.8,
  verbose = FALSE
)

sc_data <- RunUMAP(
  sc_data,
  reduction = "harmony",
  dims = 1:30,
  seed.use = 123,
  verbose = FALSE
)

source(
  "input/ScType/gene_sets_prepare.R"
)

source(
  "input/ScType/sctype_score.R"
)

sctype_gene_sets <- gene_sets_prepare(
  "input/ScType/ScTypeDB_brain.xlsx",
  "Brain"
)

sct_scaled_expression <- GetAssayData(
  sc_data,
  assay = "SCT",
  layer = "scale.data"
)

sctype_scores <- sctype_score(
  scRNAseqData = sct_scaled_expression,
  scaled = TRUE,
  gs = sctype_gene_sets$gs_positive,
  gs2 = sctype_gene_sets$gs_negative
)

cluster_ids <- as.character(
  sc_data$seurat_clusters
)

cluster_annotations <- bind_rows(
  lapply(
    unique(cluster_ids),
    function(cluster_id) {
      cluster_cells <- colnames(sc_data)[
        cluster_ids == cluster_id
      ]

      cluster_scores <- sort(
        rowSums(
          sctype_scores[
            ,
            cluster_cells,
            drop = FALSE
          ]
        ),
        decreasing = TRUE
      )

      data.frame(
        cluster = cluster_id,
        cell_type = names(cluster_scores)[1],
        score = cluster_scores[1],
        cell_count = length(cluster_cells)
      )
    }
  )
)

cluster_annotations$cell_type[
  cluster_annotations$score <
    cluster_annotations$cell_count / 4
] <- "Unknown"

sc_data$cell_type <- cluster_annotations$cell_type[
  match(
    cluster_ids,
    cluster_annotations$cluster
  )
]

fig2a_umap_clusters <- DimPlot(
  sc_data,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

fig2b_cell_types <- DimPlot(
  sc_data,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
)

cell_composition <- sc_data@meta.data %>%
  count(
    sample,
    cell_type,
    name = "cell_count"
  ) %>%
  group_by(sample) %>%
  mutate(
    fraction = cell_count / sum(cell_count)
  ) %>%
  ungroup()

fig2c_cell_composition <- ggplot(
  cell_composition,
  aes(
    x = sample,
    y = fraction,
    fill = cell_type
  )
) +
  geom_col(
    width = 0.8
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Cell fraction",
    fill = "Cell type"
  )

gbm_cells <- subset(
  sc_data,
  subset = cell_type == "GBM Cells"
)

DefaultAssay(gbm_cells) <- "RNA"

gbm_cells <- NormalizeData(
  gbm_cells,
  verbose = FALSE
)

gbm_cells <- FindVariableFeatures(
  gbm_cells,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

gbm_cells <- ScaleData(
  gbm_cells,
  features = VariableFeatures(gbm_cells),
  verbose = FALSE
)

gbm_cells <- RunPCA(
  gbm_cells,
  features = VariableFeatures(gbm_cells),
  npcs = 30,
  seed.use = 123,
  verbose = FALSE
)

gbm_cells <- RunHarmony(
  gbm_cells,
  group.by.vars = "sample",
  reduction.use = "pca",
  dims.use = 1:30,
  plot_convergence = FALSE
)

gbm_cells <- FindNeighbors(
  gbm_cells,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

gbm_cells <- FindClusters(
  gbm_cells,
  resolution = 0.8,
  verbose = FALSE
)

gbm_cells <- RunUMAP(
  gbm_cells,
  reduction = "harmony",
  dims = 1:30,
  seed.use = 123,
  verbose = FALSE
)

gbm_mes_signature <- read.csv(
  "input/signatures/Verhaak_GBM_MES_signature.csv",
  check.names = FALSE
)$gene_symbol

gbm_mes_signature <- unique(
  na.omit(gbm_mes_signature)
)

gbm_mes_signature <- CaseMatch(
  search = gbm_mes_signature,
  match = rownames(gbm_cells)
)

gbm_mes_gene_set <- list(
  GBM_MES = gbm_mes_signature
)

gbm_expression <- GetAssayData(
  gbm_cells,
  assay = "RNA",
  layer = "data"
)

gbm_rankings <- AUCell_buildRankings(
  gbm_expression,
  plotStats = FALSE
)

gbm_auc <- AUCell_calcAUC(
  gbm_mes_gene_set,
  gbm_rankings
)

gbm_cells$gbm_mes_score <- as.numeric(
  getAUC(gbm_auc)["GBM_MES", ]
)

gbm_mes_metadata <- gbm_cells@meta.data %>%
  rownames_to_column("cell") %>%
  group_by(sample) %>%
  mutate(
    rank = rank(
      gbm_mes_score,
      ties.method = "first"
    ),
    gbm_mes_group = if_else(
      rank > floor(n() / 2),
      "High GBM MES",
      "Low GBM MES"
    )
  ) %>%
  ungroup()

gbm_cells$gbm_mes_group <- gbm_mes_metadata$gbm_mes_group[
  match(
    colnames(gbm_cells),
    gbm_mes_metadata$cell
  )
]

gbm_cells$gbm_mes_group <- factor(
  gbm_cells$gbm_mes_group,
  levels = c(
    "Low GBM MES",
    "High GBM MES"
  )
)

fig2d_gbm_mes_groups <- DimPlot(
  gbm_cells,
  reduction = "umap",
  group.by = "gbm_mes_group"
)

pseudobulk_metadata <- gbm_cells@meta.data %>%
  rownames_to_column("cell") %>%
  count(
    sample,
    gbm_mes_group,
    name = "cell_count"
  )

valid_samples <- pseudobulk_metadata %>%
  filter(
    cell_count >= 5
  ) %>%
  count(sample) %>%
  filter(
    n == 2
  ) %>%
  pull(sample)

gbm_cells <- subset(
  gbm_cells,
  subset = sample %in% valid_samples
)

gbm_cells$pseudobulk_group <- interaction(
  gbm_cells$sample,
  gbm_cells$gbm_mes_group,
  sep = "__",
  drop = TRUE
)

gbm_counts <- GetAssayData(
  gbm_cells,
  assay = "RNA",
  layer = "counts"
)

pseudobulk_design_matrix <- sparse.model.matrix(
  ~ 0 + pseudobulk_group,
  data = gbm_cells@meta.data
)

colnames(
  pseudobulk_design_matrix
) <- levels(
  gbm_cells$pseudobulk_group
)

pseudobulk_counts <- gbm_counts %*%
  pseudobulk_design_matrix

pseudobulk_counts <- round(
  as.matrix(pseudobulk_counts)
)

pseudobulk_annotation <- tibble(
  pseudobulk_group = colnames(
    pseudobulk_counts
  )
) %>%
  separate(
    pseudobulk_group,
    into = c(
      "sample",
      "gbm_mes_group"
    ),
    sep = "__",
    remove = FALSE
  )

pseudobulk_annotation$sample <- factor(
  pseudobulk_annotation$sample
)

pseudobulk_annotation$gbm_mes_group <- factor(
  pseudobulk_annotation$gbm_mes_group,
  levels = c(
    "Low GBM MES",
    "High GBM MES"
  )
)

pseudobulk_dge <- DGEList(
  counts = pseudobulk_counts
)

pseudobulk_design <- model.matrix(
  ~ sample + gbm_mes_group,
  data = pseudobulk_annotation
)

expressed_genes <- filterByExpr(
  pseudobulk_dge,
  design = pseudobulk_design
)

pseudobulk_dge <- pseudobulk_dge[
  expressed_genes,
  ,
  keep.lib.sizes = FALSE
]

pseudobulk_dge <- calcNormFactors(
  pseudobulk_dge,
  method = "TMM"
)

pseudobulk_dge <- estimateDisp(
  pseudobulk_dge,
  pseudobulk_design
)

pseudobulk_fit <- glmQLFit(
  pseudobulk_dge,
  pseudobulk_design,
  robust = TRUE
)

gbm_mes_coefficient <- grep(
  "^gbm_mes_groupHigh",
  colnames(pseudobulk_design)
)

pseudobulk_test <- glmQLFTest(
  pseudobulk_fit,
  coef = gbm_mes_coefficient
)

gbm_mes_de_results <- topTags(
  pseudobulk_test,
  n = Inf
)$table %>%
  rownames_to_column("gene") %>%
  arrange(FDR) %>%
  mutate(
    regulation = case_when(
      FDR < 0.05 & logFC > 0.25 ~ "Up",
      FDR < 0.05 & logFC < -0.25 ~ "Down",
      TRUE ~ "Not significant"
    ),
    negative_log10_fdr = -log10(
      FDR + 1e-300
    )
  )

fig2e_gbm_mes_volcano <- ggplot(
  gbm_mes_de_results,
  aes(
    x = logFC,
    y = negative_log10_fdr,
    color = regulation
  )
) +
  geom_point(
    size = 1,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = c(
      -0.25,
      0.25
    ),
    linetype = 2
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2
  ) +
  theme_classic() +
  labs(
    x = "log2 fold change",
    y = "-log10(FDR)",
    color = NULL
  )