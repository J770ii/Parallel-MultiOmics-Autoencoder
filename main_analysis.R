# ==============================================================================
# Project: A parallel latent-space framework for interpretable multi-omics 
#          integration in microbiome studies
# Description: This script implements the PLAE framework for integrating 
#              Metagenomics and Functional Metagenomics (iHMP IBD dataset).
# ==============================================================================

# 0. Environment Setup
Sys.setenv(KERAS_BACKEND = "tensorflow")
library(reticulate)
# use_virtualenv("D:/MultiOmics_Project/env/r-keras", required = TRUE) 
library(keras3)
tf <- reticulate::import("tensorflow")
library(curatedMetagenomicData)
library(SummarizedExperiment)
library(dplyr)
library(pROC)
library(ggplot2)
library(cowplot)

results_dir <- "results/real_data/"
if(!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# 1. Data Loading (iHMP IBD Dataset)
cat("Loading real data from iHMP...\n")
tse_meta <- sampleMetadata |> filter(study_name == "HMP_2019_ibdmdb") |> returnSamples("relative_abundance", rownames = "short")
tse_metab <- sampleMetadata |> filter(study_name == "HMP_2019_ibdmdb") |> returnSamples("pathway_abundance", rownames = "short")

X_meta_raw  <- t(assay(tse_meta))
X_metab_raw <- t(assay(tse_metab))
y_raw       <- colData(tse_meta)$study_condition
subject_ids <- colData(tse_meta)$host_subject_id

common_samples <- intersect(rownames(X_meta_raw), rownames(X_metab_raw))
X_meta_raw  <- X_meta_raw[common_samples, , drop = FALSE]
X_metab_raw <- X_metab_raw[common_samples, , drop = FALSE]
y_raw       <- y_raw[match(common_samples, rownames(X_meta_raw))]
subject_ids <- subject_ids[match(common_samples, rownames(X_meta_raw))]
y <- as.numeric(y_raw == "IBD")

# 2. MMD Loss Function
mmd_loss <- function(x, y) {
  x_sq <- tf$reduce_sum(tf$square(x), axis = 1L, keepdims = TRUE)
  y_sq <- tf$reduce_sum(tf$square(y), axis = 1L, keepdims = TRUE)
  xy <- tf$matmul(x, y, transpose_b = TRUE)
  xx <- tf$matmul(x, x, transpose_b = TRUE)
  yy <- tf$matmul(y, y, transpose_b = TRUE)
  dxx <- x_sq - 2 * xx + tf$transpose(x_sq)
  dyy <- y_sq - 2 * yy + tf$transpose(y_sq)
  dxy <- x_sq - 2 * xy + tf$transpose(y_sq)
  gamma <- 1.0 / 16.0
  kxx <- tf$exp(-gamma * dxx); kyy <- tf$exp(-gamma * dyy); kxy <- tf$exp(-gamma * dxy)
  tf$reduce_mean(kxx) + tf$reduce_mean(kyy) - 2 * tf$reduce_mean(kxy)
}

# 3. Patient-Level 5-Fold CV (No Data Leakage)
set.seed(42)
patients <- unique(subject_ids)
n_folds <- 5
fold_patients <- split(patients, cut(seq_along(patients), n_folds, labels = FALSE))

auc_results <- numeric(n_folds)
oof_y_test <- c(); oof_prob_test <- c()

for (fold in 1:n_folds) {
  test_p <- fold_patients[[fold]]; train_p <- setdiff(patients, test_p)
  train_idx <- which(subject_ids %in% train_p); test_idx <- which(subject_ids %in% test_p)
  
  # Feature Selection (Train only)
  p_meta <- apply(X_meta_raw[train_idx, , drop = FALSE], 2, function(col) wilcox.test(col ~ y[train_idx])$p.value)
  p_metab <- apply(X_metab_raw[train_idx, , drop = FALSE], 2, function(col) wilcox.test(col ~ y[train_idx])$p.value)
  top_meta <- order(p_meta)[1:min(200, length(p_meta))]; top_metab <- order(p_metab)[1:min(100, length(p_metab))]
  
  X_tr_m <- X_meta_raw[train_idx, top_meta, drop = FALSE]; X_te_m <- X_meta_raw[test_idx, top_meta, drop = FALSE]
  X_tr_p <- X_metab_raw[train_idx, top_metab, drop = FALSE]; X_te_p <- X_metab_raw[test_idx, top_metab, drop = FALSE]
  
  # Min-Max Scaling (Train only)
  min_m <- apply(X_tr_m, 2, min); max_m <- apply(X_tr_m, 2, max)
  X_tr_m <- sweep(sweep(X_tr_m, 2, min_m, "-"), 2, max_m - min_m, "/"); X_te_m <- sweep(sweep(X_te_m, 2, min_m, "-"), 2, max_m - min_m, "/")
  min_p <- apply(X_tr_p, 2, min); max_p <- apply(X_tr_p, 2, max)
  X_tr_p <- sweep(sweep(X_tr_p, 2, min_p, "-"), 2, max_p - min_p, "/"); X_te_p <- sweep(sweep(X_te_p, 2, min_p, "-"), 2, max_p - min_p, "/")
  
  # Build Model
  input_meta <- layer_input(shape = ncol(X_tr_m), name = "meta_input")
  enc_meta <- input_meta |> layer_dense(32, activation = "relu", kernel_regularizer = regularizer_l2(0.01)) |> layer_dropout(0.3) |> layer_dense(16, name = "latent_meta")
  dec_meta <- enc_meta |> layer_dense(32, activation = "relu") |> layer_dense(ncol(X_tr_m), activation = "sigmoid", name = "recon_meta")
  
  input_metab <- layer_input(shape = ncol(X_tr_p), name = "metab_input")
  enc_metab <- input_metab |> layer_dense(32, activation = "relu", kernel_regularizer = regularizer_l2(0.01)) |> layer_dropout(0.3) |> layer_dense(16, name = "latent_metab")
  dec_metab <- enc_metab |> layer_dense(32, activation = "relu") |> layer_dense(ncol(X_tr_p), activation = "linear", name = "recon_metab")
  
  latent_concat <- layer_concatenate(list(enc_meta, enc_metab))
  pheno_out <- latent_concat |> layer_dense(16, activation = "relu") |> layer_dense(1, activation = "sigmoid", name = "pheno")
  
  model <- keras_model(
    inputs = list(input_meta, input_metab),
    outputs = list(recon_meta = dec_meta, recon_metab = dec_metab, pheno = pheno_out, latent_meta = enc_meta, latent_metab = enc_metab)
  )
  
  optimizer <- optimizer_adam(learning_rate = 0.001)
  n_train <- length(train_idx); n_batches <- ceiling(n_train / 32)
  
  # Training Loop
  for (epoch in 1:40) {
    idx <- sample(n_train)
    for (b in 1:n_batches) {
      bi <- idx[((b-1)*32 + 1):min(b*32, n_train)]
      xm  <- tf$convert_to_tensor(X_tr_m[bi, , drop = FALSE], dtype = "float32")
      xmb <- tf$convert_to_tensor(X_tr_p[bi, , drop = FALSE], dtype = "float32")
      yb  <- tf$convert_to_tensor(matrix(y[train_idx][bi], ncol = 1), dtype = "float32")
      with(tf$GradientTape() %as% tape, {
        preds <- model(list(xm, xmb), training = TRUE)
        l_rec_m <- tf$reduce_mean(tf$square(xm - preds[[1]])); l_rec_p <- tf$reduce_mean(tf$square(xmb - preds[[2]]))
        l_pheno <- tf$reduce_mean(tf$keras$losses$binary_crossentropy(yb, preds[[3]]))
        l_align <- mmd_loss(preds[[4]], preds[[5]])
        total_loss <- l_rec_m + l_rec_p + l_pheno + 0.1 * l_align
      })
      grads <- tape$gradient(total_loss, model$trainable_weights)
      optimizer$apply_gradients(zip_lists(grads, model$trainable_weights))
    }
  }
  
  # Evaluation
  xm_te <- tf$convert_to_tensor(X_te_m, dtype = "float32"); xmb_te <- tf$convert_to_tensor(X_te_p, dtype = "float32")
  preds_te <- model(list(xm_te, xmb_te), training = FALSE)
  prob_te <- as.numeric(preds_te[[3]])
  
  oof_y_test <- c(oof_y_test, y[test_idx]); oof_prob_test <- c(oof_prob_test, prob_te)
  roc_cv <- roc(y[test_idx], prob_te, quiet = TRUE)
  auc_results[fold] <- as.numeric(auc(roc_cv))
}

cat("Mean AUC:", round(mean(auc_results), 4), "±", round(sd(auc_results), 4), "\n")

# (Plotting code for ROC, PCA, and Feature Importance follows as previously provided)