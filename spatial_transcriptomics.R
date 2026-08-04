library(tidyverse)
library(Seurat)
library(AUCell)

set.seed(123)

spatial_data <- Load10X_Spatial(
  data.dir = "input/GSE194329/T1",
  filename = "filtered_feature_bc_matrix.h5"
)

DefaultAssay(spatial_data) <- "Spatial"

spatial_data <- NormalizeData(
  spatial_data,
  verbose = FALSE
)

spatial_data <- FindVariableFeatures(
  spatial_data,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

spatial_data <- ScaleData(
  spatial_data,
  features = rownames(spatial_data),
  verbose = FALSE
)

spatial_data <- RunPCA(
  spatial_data,
  npcs = 50,
  seed.use = 123,
  verbose = FALSE
)

spatial_data <- FindNeighbors(
  spatial_data,
  reduction = "pca",
  dims = 1:20,
  verbose = FALSE
)

spatial_data <- FindClusters(
  spatial_data,
  resolution = 0.8,
  verbose = FALSE
)

spatial_data <- RunUMAP(
  spatial_data,
  reduction = "pca",
  dims = 1:20,
  seed.use = 123,
  verbose = FALSE
)

source("input/ScType/gene_sets_prepare.R")
source("input/ScType/sctype_score.R")

sctype_gene_sets <- gene_sets_prepare(
  "input/ScType/ScTypeDB_brain.xlsx",
  "Brain"
)

scaled_spatial_expression <- GetAssayData(
  spatial_data,
  assay = "Spatial",
  layer = "scale.data"
)

sctype_scores <- sctype_score(
  scRNAseqData = scaled_spatial_expression,
  scaled = TRUE,
  gs = sctype_gene_sets$gs_positive,
  gs2 = sctype_gene_sets$gs_negative
)

cluster_ids <- as.character(
  spatial_data$seurat_clusters
)

cluster_annotations <- bind_rows(
  lapply(
    sort(unique(cluster_ids)),
    function(cluster_id) {
      cluster_spots <- colnames(spatial_data)[
        cluster_ids == cluster_id
      ]

      annotation_scores <- sort(
        rowSums(
          sctype_scores[
            ,
            cluster_spots,
            drop = FALSE
          ]
        ),
        decreasing = TRUE
      )

      data.frame(
        cluster = cluster_id,
        annotation = names(annotation_scores)[1],
        score = annotation_scores[1],
        spot_count = length(cluster_spots)
      )
    }
  )
)

cluster_annotations$annotation[
  cluster_annotations$score <
    cluster_annotations$spot_count / 4
] <- "Unknown"

spatial_data$spatial_annotation <- cluster_annotations$annotation[
  match(
    cluster_ids,
    cluster_annotations$cluster
  )
]

gbm_mes_signature <- read.csv(
  "input/signatures/Verhaak_GBM_MES_signature.csv",
  check.names = FALSE
)$gene_symbol

gbm_mes_signature <- unique(
  na.omit(gbm_mes_signature)
)

gbm_mes_signature <- CaseMatch(
  search = gbm_mes_signature,
  match = rownames(spatial_data)
)

spatial_expression <- GetAssayData(
  spatial_data,
  assay = "Spatial",
  layer = "data"
)

spatial_rankings <- AUCell_buildRankings(
  spatial_expression,
  plotStats = FALSE,
  verbose = FALSE
)

spatial_auc <- AUCell_calcAUC(
  list(GBM_MES = gbm_mes_signature),
  spatial_rankings,
  aucMaxRank = ceiling(
    0.05 * nrow(spatial_rankings)
  )
)

spatial_data$gbm_mes_score <- as.numeric(
  getAUC(spatial_auc)[
    "GBM_MES",
    colnames(spatial_data)
  ]
)

itgb2_gene <- CaseMatch(
  search = "ITGB2",
  match = rownames(spatial_data)
)

fig3a_spatial_clusters <- SpatialDimPlot(
  spatial_data,
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  image.alpha = 1
)

fig3b_spatial_annotation <- SpatialDimPlot(
  spatial_data,
  group.by = "spatial_annotation",
  label = FALSE,
  image.alpha = 1
)

fig3c_spatial_mes_score <- SpatialFeaturePlot(
  spatial_data,
  features = "gbm_mes_score",
  image.alpha = 1,
  combine = FALSE
)[[1]] +
  scale_fill_gradient(
    low = "#6a87b7",
    high = "#e2af61"
  )

fig3d_umap_clusters <- DimPlot(
  spatial_data,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
)

fig3e_umap_annotation <- DimPlot(
  spatial_data,
  reduction = "umap",
  group.by = "spatial_annotation",
  label = TRUE,
  repel = TRUE
)

fig3f_itgb2_umap <- FeaturePlot(
  spatial_data,
  features = itgb2_gene,
  reduction = "umap",
  order = TRUE,
  combine = FALSE
)[[1]] +
  scale_color_gradient(
    low = "#6a87b7",
    high = "#e2af61"
  )

fig3g_itgb2_spatial <- SpatialFeaturePlot(
  spatial_data,
  features = itgb2_gene,
  image.alpha = 1,
  combine = FALSE
)[[1]] +
  scale_fill_gradient(
    low = "#6a87b7",
    high = "#e2af61"
  )
