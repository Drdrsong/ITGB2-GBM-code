library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(harmony)
library(AUCell)
library(edgeR)
library(Matrix)

set.seed(123)

sample_directories <- sort(
  list.dirs(
    "input/GSE162631",
    recursive = FALSE,
    full.names = TRUE
  )
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

sc_data$log10_genes_per_umi <- log10(
  sc_data$nFeature_RNA
) / log10(
  sc_data$nCount_RNA
)

sc_data$mito_ratio <- PercentageFeatureSet(
  sc_data,
  pattern = "^MT-"
) / 100

sc_data <- subset(
  sc_data,
  subset =
    nCount_RNA > 200 &
    nCount_RNA < 50000 &
    nFeature_RNA > 200 &
    log10_genes_per_umi > 0.8 &
    mito_ratio < 0.2
)

count_matrix <- LayerData(
  sc_data,
  assay = "RNA",
  layer = "counts"
)

sc_data <- CreateSeuratObject(
  counts = count_matrix[
    rowSums(count_matrix > 0) >= 10,
    ,
    drop = FALSE
  ],
  meta.data = sc_data@meta.data
)

sc_data <- NormalizeData(
  sc_data,
  verbose = FALSE
)

single_cell_experiment <- as.SingleCellExperiment(
  sc_data
)

single_cell_experiment <- scDblFinder(
  single_cell_experiment,
  samples = "sample"
)

sc_data <- as.Seurat(
  single_cell_experiment
)

sc_data <- subset(
  sc_data,
  subset = scDblFinder.class == "singlet"
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

sc_data$cell_cycle_difference <-
  sc_data$S.Score - sc_data$G2M.Score

sc_data <- SCTransform(
  sc_data,
  vars.to.regress = c(
    "S.Score",
    "G2M.Score",
    "cell_cycle_difference",
    "mito_ratio"
  ),
  verbose = FALSE
)

sc_data <- RunPCA(
  sc_data,
  assay = "SCT",
  features = VariableFeatures(sc_data),
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

source("input/ScType/gene_sets_prepare.R")
source("input/ScType/sctype_score.R")

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

gbm_cells <- subset(
  sc_data,
  subset = cell_type == "GBM Cells"
)

DefaultAssay(gbm_cells) <- "RNA"

gbm_cells <- NormalizeData(
  gbm_cells,
  assay = "RNA",
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
  getAUC(gbm_auc)[
    "GBM_MES",
    colnames(gbm_cells)
  ]
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
      "high_GBM_MES",
      "low_GBM_MES"
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
    "low_GBM_MES",
    "high_GBM_MES"
  )
)

min_cells_per_group <- 5

pseudobulk_metadata <- gbm_cells@meta.data %>%
  as.data.frame() %>%
  count(
    sample,
    gbm_mes_group,
    name = "cell_count"
  ) %>%
  pivot_wider(
    names_from = gbm_mes_group,
    values_from = cell_count,
    values_fill = 0
  )

valid_samples <- pseudobulk_metadata %>%
  filter(
    low_GBM_MES >= min_cells_per_group,
    high_GBM_MES >= min_cells_per_group
  ) %>%
  pull(sample)

if (length(valid_samples) == 0) {
  valid_samples <- pseudobulk_metadata %>%
    filter(
      low_GBM_MES > 0,
      high_GBM_MES > 0
    ) %>%
    pull(sample)
}

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

pseudobulk_counts <- round(
  as.matrix(
    gbm_counts %*%
      pseudobulk_design_matrix
  )
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

valid_samples <- pseudobulk_annotation %>%
  distinct(
    sample,
    gbm_mes_group
  ) %>%
  count(sample) %>%
  filter(n == 2) %>%
  pull(sample)

keep_columns <- pseudobulk_annotation$sample %in%
  valid_samples

pseudobulk_counts <- pseudobulk_counts[
  ,
  keep_columns,
  drop = FALSE
]

pseudobulk_annotation <- pseudobulk_annotation[
  keep_columns,
  ,
  drop = FALSE
]

pseudobulk_annotation$sample <- factor(
  pseudobulk_annotation$sample
)

pseudobulk_annotation$gbm_mes_group <- factor(
  pseudobulk_annotation$gbm_mes_group,
  levels = c(
    "low_GBM_MES",
    "high_GBM_MES"
  )
)

pseudobulk_dge <- DGEList(
  counts = pseudobulk_counts
)

pseudobulk_dge <- calcNormFactors(
  pseudobulk_dge,
  method = "TMM"
)

pseudobulk_design <- model.matrix(
  ~ sample + gbm_mes_group,
  data = pseudobulk_annotation
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

gbm_mes_coefficient <- which(
  colnames(pseudobulk_design) ==
    "gbm_mes_grouphigh_GBM_MES"
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
  arrange(FDR)
