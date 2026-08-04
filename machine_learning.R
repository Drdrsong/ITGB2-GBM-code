library(tidyverse)
library(glmnet)
library(caret)
library(klaR)
library(Boruta)
library(UpSetR)
library(e1071)
library(randomForest)
library(ipred)

set.seed(123)

machine_learning_genes <- intersect(
  significant_logistic_results$gene,
  rownames(tcga_expression)
)

machine_learning_x <- as.data.frame(
  t(
    tcga_expression[
      machine_learning_genes,
      rownames(tcga_mes_metadata),
      drop = FALSE
    ]
  ),
  check.names = FALSE
)

machine_learning_y <- tcga_mes_metadata[
  rownames(machine_learning_x),
  "mes_group"
]

machine_learning_y <- factor(
  machine_learning_y,
  levels = c("Low", "High")
)

complete_samples <- complete.cases(
  machine_learning_x,
  machine_learning_y
)

machine_learning_x <- machine_learning_x[
  complete_samples,
  ,
  drop = FALSE
]

machine_learning_y <- machine_learning_y[
  complete_samples
]

run_repeated_lasso <- function(
    x,
    y,
    nfolds = 5,
    iterations = 1000
) {
  y_numeric <- as.integer(y == "High")
  model_signatures <- character(iterations)

  for (iteration in seq_len(iterations)) {
    set.seed(iteration)

    lasso_model <- cv.glmnet(
      x = as.matrix(x),
      y = y_numeric,
      family = "binomial",
      alpha = 1,
      nfolds = nfolds
    )

    coefficients <- coef(
      lasso_model,
      s = "lambda.min"
    )

    selected_genes <- rownames(coefficients)[
      as.vector(coefficients != 0)
    ]

    selected_genes <- setdiff(
      selected_genes,
      "(Intercept)"
    )

    model_signatures[iteration] <- paste(
      sort(selected_genes),
      collapse = "|"
    )
  }

  model_signatures <- model_signatures[
    model_signatures != ""
  ]

  model_frequency <- sort(
    table(model_signatures),
    decreasing = TRUE
  )

  frequency_data <- data.frame(
    model = names(model_frequency),
    frequency = as.integer(model_frequency),
    stringsAsFactors = FALSE
  )

  selected_genes <- strsplit(
    frequency_data$model[1],
    "\\|"
  )[[1]]

  list(
    selected_genes = selected_genes,
    frequency = frequency_data
  )
}

lasso_results <- run_repeated_lasso(
  x = machine_learning_x,
  y = machine_learning_y,
  nfolds = 5,
  iterations = 1000
)

lasso_genes <- lasso_results$selected_genes

feature_sizes <- seq(
  2,
  ncol(machine_learning_x)
)

svm_control <- rfeControl(
  functions = caretFuncs,
  method = "cv",
  number = 5
)

svm_model <- rfe(
  x = machine_learning_x,
  y = machine_learning_y,
  sizes = feature_sizes,
  rfeControl = svm_control,
  method = "svmRadial",
  tuneLength = 5
)

svm_genes <- predictors(
  svm_model
)

bayesian_control <- rfeControl(
  functions = nbFuncs,
  method = "cv",
  number = 2
)

bayesian_model <- rfe(
  x = machine_learning_x,
  y = machine_learning_y,
  sizes = feature_sizes,
  rfeControl = bayesian_control
)

bayesian_genes <- predictors(
  bayesian_model
)

boruta_model <- Boruta(
  x = machine_learning_x,
  y = machine_learning_y,
  doTrace = 0
)

boruta_model <- TentativeRoughFix(
  boruta_model
)

boruta_genes <- getSelectedAttributes(
  boruta_model,
  withTentative = FALSE
)

random_forest_control <- rfeControl(
  functions = rfFuncs,
  method = "cv",
  number = 2
)

random_forest_model <- rfe(
  x = machine_learning_x,
  y = machine_learning_y,
  sizes = feature_sizes,
  rfeControl = random_forest_control
)

random_forest_genes <- predictors(
  random_forest_model
)

bagged_tree_control <- rfeControl(
  functions = treebagFuncs,
  method = "cv",
  number = 5
)

bagged_tree_model <- rfe(
  x = machine_learning_x,
  y = machine_learning_y,
  sizes = feature_sizes,
  rfeControl = bagged_tree_control
)

bagged_tree_genes <- predictors(
  bagged_tree_model
)

lvq_control <- trainControl(
  method = "repeatedcv",
  number = 10
)

lvq_grid <- expand.grid(
  size = 6,
  k = c(1, 2)
)

lvq_model <- train(
  x = machine_learning_x,
  y = machine_learning_y,
  method = "lvq",
  preProcess = "scale",
  trControl = lvq_control,
  tuneGrid = lvq_grid
)

lvq_importance <- varImp(
  lvq_model,
  scale = FALSE
)

if (ncol(lvq_importance$importance) == 1) {
  lvq_scores <- lvq_importance$importance[, 1]
} else {
  lvq_scores <- apply(
    lvq_importance$importance,
    1,
    max
  )
}

lvq_genes <- names(
  lvq_scores[
    lvq_scores > 0.5
  ]
)

selected_features <- list(
  LASSO = unique(lasso_genes),
  SVM = unique(svm_genes),
  Bayesian = unique(bayesian_genes),
  Boruta = unique(boruta_genes),
  Random_Forest = unique(random_forest_genes),
  Bagged_Trees = unique(bagged_tree_genes),
  LVQ = unique(lvq_genes)
)

consensus_gene_frequency <- table(
  unlist(selected_features)
)

consensus_hub_genes <- sort(
  names(
    consensus_gene_frequency[
      consensus_gene_frequency >= 3
    ]
  )
)

fig1g_lasso <- lasso_results$frequency %>%
  slice_head(n = 20) %>%
  mutate(
    model = forcats::fct_reorder(
      model,
      frequency
    )
  ) %>%
  ggplot(
    aes(
      x = model,
      y = frequency
    )
  ) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Frequency",
    title = "LASSO"
  )

plot(
  svm_model,
  type = "b",
  main = "SVM"
)

plot(
  bayesian_model,
  type = "b",
  main = "Bayesian model"
)

plot(
  boruta_model,
  las = 2,
  main = "Boruta"
)

plot(
  random_forest_model,
  type = "b",
  main = "Random Forest"
)

plot(
  bagged_tree_model,
  type = "b",
  main = "Bagged Trees"
)

plot(
  lvq_importance,
  main = "LVQ"
)

fig1n_upset_data <- fromList(
  selected_features
)

upset(
  fig1n_upset_data,
  nsets = 7,
  keep.order = TRUE,
  order.by = "freq"
)