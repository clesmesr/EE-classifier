## =========================================================
##  MAKE ALL FIGURES (synthetic + real) — parallel & reproducible
##  - Generates panels: curves + EE_E + EE_H
##  - Optionally reruns outer-CV and saves results boxplots
##  - File names match the manuscript exactly
## =========================================================

## ========== SETTINGS ==========
MAKE_RESULTS   <- TRUE      # FALSE = só painéis (rápido); TRUE = também boxplots (mais demorado)
N_CORES        <- 32        # nº de cores para paralelização
SEED_GLOBAL    <- 45        # semente global para reprodutibilidade
INNER_FOLDS    <- 5         # CV interno (caret)
K_SYNTH        <- 100       # folds sintéticos
K_GROWTH       <- 20
K_TECATOR      <- 50
K_MCO          <- 30

## ========== PACKAGES ==========
req <- c("roahd","fda.usc","caret","ggplot2","lattice",
         "doParallel","foreach","parallel","data.table","xtable","kernlab","randomForest")
for (p in req) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
lapply(req, require, character.only = TRUE)

set.seed(SEED_GLOBAL)
if (!dir.exists("images")) dir.create("images", recursive = TRUE)
## ========== UTIL: fast Fraiman–Muniz depth (no extra pkgs) ==========
## FM(z; R) = average_t [ 1 - |0.5 - F_t(z(t))| ], F_t empirical CDF
fm_against_ref <- function(ref, newdata) {
  n_ref <- nrow(ref); P <- ncol(ref); n_new <- nrow(newdata)
  out <- numeric(n_new)
  ref <- as.matrix(ref); newdata <- as.matrix(newdata)
  for (i in 1:n_new) {
    x  <- newdata[i, ]
    u  <- colMeans(ref <= matrix(x, nrow = n_ref, ncol = P, byrow = TRUE))  # F_t(x)
    dt <- 1 - abs(0.5 - u)  # correct FM scaling (range [0.5, 1])
    out[i] <- mean(dt)
  }
  out
}



## ========== EE features (no leakage helpers) ==========
mhi_against_ref <- function(ref, newdata) {
  n_ref <- nrow(ref); P <- ncol(ref); n_new <- nrow(newdata)
  out <- numeric(n_new)
  for (i in 1:n_new) {
    x <- matrix(newdata[i, ], nrow = n_ref, ncol = P, byrow = TRUE)
    out[i] <- mean(ref <= x)
  }
  out
}
mei_against_ref <- function(ref, newdata) {
  n_ref <- nrow(ref); P <- ncol(ref); n_new <- nrow(newdata)
  out <- numeric(n_new)
  for (i in 1:n_new) {
    x <- matrix(newdata[i, ], nrow = n_ref, ncol = P, byrow = TRUE)
    out[i] <- mean(ref >= x)
  }
  out
}

## ========== Plot helpers
plot_curves_two_classes <- function(X, Y,
                                    xlab = "time index",
                                    ylab = "value",
                                    legend_labels = c("class 0", "class 1")) {
  
  P <- ncol(X)
  tgrid <- seq_len(P)
  
  to_df <- function(M, g_label) {
    data.frame(
      id = rep(seq_len(nrow(M)), each = ncol(M)),
      t  = rep(tgrid, times = nrow(M)),
      y  = as.vector(t(M)),
      g  = g_label
    )
  }
  
  df <- rbind(
    to_df(X, "0"),
    to_df(Y, "1")
  )
  
  ggplot(df, aes(
    x = t,
    y = y,
    group = interaction(g, id),
    color = g
  )) +
    geom_line(alpha = 0.6) +
    scale_color_manual(
      values = c("0" = "black", "1" = "grey50"),
      labels = legend_labels,
      name = "group"
    ) +
    labs(x = xlab, y = ylab) +
    theme_minimal() +
    theme(
      legend.position = "none" 
    )
}

plot_ee_full <- function(X, Y, index = c("epi","hipo")) {
  index <- match.arg(index)
  if (index == "epi") {
    xx <- mei_against_ref(X, rbind(X, Y)); yy <- mei_against_ref(Y, rbind(X, Y))
    xlab <- "mei (ref class 0)"; ylab <- "mei (ref class 1)"
  } else {
    xx <- mhi_against_ref(X, rbind(X, Y)); yy <- mhi_against_ref(Y, rbind(X, Y))
    xlab <- "mhi (ref class 0)"; ylab <- "mhi (ref class 1)"
  }
  g <- factor(c(rep(0, nrow(X)), rep(1, nrow(Y))))
  df <- data.frame(x = xx, y = yy, g = g)
  ggplot(df, aes(x = x, y = y, shape = g)) +
    geom_point(size = 2) +
    scale_shape_manual(values = c("0" = 4, "1" = 1), name = "group",
                       labels = c("class 0","class 1")) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    labs(x = xlab, y = ylab) +
    theme_minimal() +
    theme(legend.position = "none")
}
plot_results_box <- function(folds_df, outfile_pdf) {
  method_levels <- c("EEh_lda","EEh_qda","EEh_knn","EEh_svm","EEh_rf",
                     "EEe_lda","EEe_qda","EEe_knn","EEe_svm","EEe_rf",
                     "DD_lda","DD_qda","DD_knn")
  folds_df$model <- factor(folds_df$model, levels = method_levels)
  p <- ggplot(folds_df, aes(x = model, y = acc)) +
    geom_boxplot(outlier.size = 0.8) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "method", y = "outer-cv accuracy") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(outfile_pdf, p, width = 6, height = 5, units = "in")
}

## ========== Outer CV (parallel-aware) ==========
make_stratified_folds <- function(n_x, n_y, K, seed = 7) {
  set.seed(seed)
  f0 <- caret::createFolds(factor(rep(0, n_x)), k = K, list = TRUE)
  f1 <- caret::createFolds(factor(rep(1, n_y)), k = K, list = TRUE)
  list(f0 = f0, f1 = f1)
}

run_outer_cv <- function(X, Y, K = 30, seed = 45,
                         models_ee = c("lda","qda","knn","svmRadial","rf"),
                         models_dd = c("lda","qda","knn"),
                         inner_number = 5,
                         do_checks = TRUE,
                         allow_parallel = TRUE) {
  rownames(X) <- if (is.null(rownames(X))) paste0("X_", seq_len(nrow(X))) else rownames(X)
  rownames(Y) <- if (is.null(rownames(Y))) paste0("Y_", seq_len(nrow(Y))) else rownames(Y)
  
  folds <- make_stratified_folds(nrow(X), nrow(Y), K, seed)
  acc_list <- list()
  P <- ncol(X)
  tol <- 1/P + 1e-6
  
  inner_ctrl <- caret::trainControl(method = "cv", number = inner_number, allowParallel = allow_parallel)
  
  for (k in seq_len(K)) {
    idx_val_X <- folds$f0[[k]]
    idx_val_Y <- folds$f1[[k]]
    
    X_train <- X[-idx_val_X, , drop = FALSE]; X_val <- X[idx_val_X, , drop = FALSE]
    Y_train <- Y[-idx_val_Y, , drop = FALSE]; Y_val <- Y[idx_val_Y, , drop = FALSE]
    
    if (do_checks) {
      stopifnot(length(intersect(rownames(X_train), rownames(X_val))) == 0)
      stopifnot(length(intersect(rownames(Y_train), rownames(Y_val))) == 0)
    }
    
    eeh_train <- data.frame(
      x = mhi_against_ref(X_train, rbind(X_train, Y_train)),
      y = mhi_against_ref(Y_train, rbind(X_train, Y_train)),
      g = factor(c(rep(0, nrow(X_train)), rep(1, nrow(Y_train))))
    )
    eeh_val <- data.frame(
      x = mhi_against_ref(X_train, rbind(X_val, Y_val)),
      y = mhi_against_ref(Y_train, rbind(X_val, Y_val)),
      g = factor(c(rep(0, nrow(X_val)), rep(1, nrow(Y_val))))
    )
    if (do_checks) {
      mei_x_val <- mei_against_ref(X_train, rbind(X_val, Y_val))
      mei_y_val <- mei_against_ref(Y_train, rbind(X_val, Y_val))
      stopifnot(median(abs(eeh_val$x + mei_x_val - 1)) < tol)
      stopifnot(median(abs(eeh_val$y + mei_y_val - 1)) < tol)
    }
    
    eee_train <- data.frame(
      x = mei_against_ref(X_train, rbind(X_train, Y_train)),
      y = mei_against_ref(Y_train, rbind(X_train, Y_train)),
      g = eeh_train$g
    )
    eee_val <- data.frame(
      x = mei_against_ref(X_train, rbind(X_val, Y_val)),
      y = mei_against_ref(Y_train, rbind(X_val, Y_val)),
      g = eeh_val$g
    )
    
    dd_train <- data.frame(
      x = fm_against_ref(X_train, rbind(X_train, Y_train)),
      y = fm_against_ref(Y_train, rbind(X_train, Y_train)),
      g = eeh_train$g
    )
    dd_val <- data.frame(
      x = fm_against_ref(X_train, rbind(X_val, Y_val)),
      y = fm_against_ref(Y_train, rbind(X_val, Y_val)),
      g = eeh_val$g
    )
    
    ## --- fit/predict loops (these serão paralelizados via caret) ---
    for (m in models_ee) {
      set.seed(7)
      mname <- if (m == "svmRadial") "EEh_svm" else paste0("EEh_", m)
      fit <- caret::train(eeh_train[, 1:2], eeh_train$g, method = m,
                          metric = "Accuracy", trControl = inner_ctrl)
      pred <- predict(fit, eeh_val[, 1:2])
      acc  <- mean(pred == eeh_val$g)
      acc_list[[length(acc_list) + 1]] <- data.frame(fold = k, model = mname, acc = acc)
    }
    for (m in models_ee) {
      set.seed(7)
      mname <- if (m == "svmRadial") "EEe_svm" else paste0("EEe_", m)
      fit <- caret::train(eee_train[, 1:2], eee_train$g, method = m,
                          metric = "Accuracy", trControl = inner_ctrl)
      pred <- predict(fit, eee_val[, 1:2])
      acc  <- mean(pred == eee_val$g)
      acc_list[[length(acc_list) + 1]] <- data.frame(fold = k, model = mname, acc = acc)
    }
    for (m in models_dd) {
      set.seed(7)
      mname <- paste0("DD_", m)
      fit <- caret::train(dd_train[, 1:2], dd_train$g, method = m,
                          metric = "Accuracy", trControl = inner_ctrl)
      pred <- predict(fit, dd_val[, 1:2])
      acc  <- mean(pred == dd_val$g)
      acc_list[[length(acc_list) + 1]] <- data.frame(fold = k, model = mname, acc = acc)
    }
  }
  
  acc_tbl <- data.table::rbindlist(acc_list)
  summary_tbl <- do.call(rbind, lapply(split(acc_tbl$acc, acc_tbl$model), function(v) {
    c(`Min.` = min(v, na.rm = TRUE),
      `1st Qu.` = as.numeric(quantile(v, .25, na.rm = TRUE)),
      `Median` = median(v, na.rm = TRUE),
      `Mean` = mean(v, na.rm = TRUE),
      `3rd Qu.` = as.numeric(quantile(v, .75, na.rm = TRUE)),
      `Max` = max(v, na.rm = TRUE),
      `NA's` = sum(is.na(v)))
  }))
  list(summary = as.data.frame(summary_tbl), folds = as.data.frame(acc_tbl))
}

pretty_summary <- function(summary_tbl) {
  rn <- rownames(summary_tbl)
  keep_order <- c("EEh_lda","EEh_qda","EEh_knn","EEh_svm","EEh_rf",
                  "EEe_lda","EEe_qda","EEe_knn","EEe_svm","EEe_rf",
                  "DD_lda","DD_qda","DD_knn")
  keep <- intersect(keep_order, rn)
  out <- summary_tbl[keep, , drop = FALSE]
  out
}

## ========== Synthetic data generators + naming to match manuscript ==========
exp_base_name <- function(e) {
  switch(as.character(e),
         "1" = "2_sen_sen",
         "2" = "1_sen_cos",
         "3" = "3_sen_sen",
         "4" = "4_amp",
         "5" = "5_rui",
         "6" = "6_func")
}
generate_synthetic_experiment <- function(exp_id, N = 200, P = 100, seed = SEED_GLOBAL) {
  set.seed(seed)
  tgrid <- seq(-1, 1, length.out = P)
  if (exp_id == 1) { mu1 <- sin(pi * tgrid);     mu2 <- sin(pi * tgrid);     alpha <- 0.2; beta <- 0.3
  } else if (exp_id == 2) { mu1 <- sin(pi * tgrid);     mu2 <- cos(pi * tgrid);     alpha <- 0.2; beta <- 0.3
  } else if (exp_id == 3) { mu1 <- sin(pi * tgrid) - 2; mu2 <- sin(pi * tgrid);     alpha <- 0.2; beta <- 0.3
  } else if (exp_id == 4) { mu1 <- sin(2*pi * tgrid);   mu2 <- cos(2*pi * tgrid);   alpha <- 0.7; beta <- 0.3
  } else if (exp_id == 5) { mu1 <- sin(2*pi * tgrid);   mu2 <- cos(2*pi * tgrid);   alpha <- 0.2; beta <- 0.9
  } else if (exp_id == 6) { mu1 <- tgrid^2 - 1;         mu2 <- tgrid^2;            alpha <- 0.2; beta <- 0.9
  } else stop("exp_id must be 1..6")
  C <- roahd::exp_cov_function(tgrid, alpha = alpha, beta = beta)
  X <- roahd::generate_gauss_fdata(N, mu1, C)
  Y <- roahd::generate_gauss_fdata(N, mu2, C)
  rownames(X) <- paste0("X_", seq_len(nrow(X))); rownames(Y) <- paste0("Y_", seq_len(nrow(Y)))
  list(X = X, Y = Y)
}

## ========== Real datasets loader ==========
load_real <- function(choice) {
  if (choice == "growth") {
    data(growth, package = "fda.usc")
    X <- t(growth$hgtf)  # girls -> class 0
    Y <- t(growth$hgtm)  # boys  -> class 1
  } else if (choice == "tecator") {
    data("tecator", package = "fda.usc")
    X <- tecator$absorp.fdata$data[tecator$y$Fat <  20, ]
    Y <- tecator$absorp.fdata$data[tecator$y$Fat >= 20, ]
  } else if (choice == "mco") {
    data(MCO, package = "fda.usc")
    X <- MCO$intact$data[MCO$classintact == 1, ]  # control
    Y <- MCO$intact$data[MCO$classintact == 2, ]  # treatment
  } else stop("choice must be 'growth' | 'tecator' | 'mco'")
  rownames(X) <- paste0("X_", seq_len(nrow(X))); rownames(Y) <- paste0("Y_", seq_len(nrow(Y)))
  list(X = X, Y = Y)
}

## ========== Parallel backend ==========
cl <- parallel::makePSOCKcluster(N_CORES)
doParallel::registerDoParallel(cl)
on.exit({ try(stopCluster(cl), silent = TRUE) }, add = TRUE)

## =========================================================
##  A) SYNTHETIC — panels + (optional) results boxplots
## =========================================================
for (e in 1:6) {
  message(sprintf(">>> synthetic experiment %d", e))
  dat <- generate_synthetic_experiment(e, N = 200, P = 100, seed = SEED_GLOBAL)
  base <- exp_base_name(e)
  
  # panels (curves + EE)
  p_curves <- plot_curves_two_classes(dat$X, dat$Y,
                                      xlab = "time index", ylab = "value",
                                      legend_labels = c("class 0","class 1"))
  ggsave(file.path("images", paste0(base, "_sim.pdf")), p_curves, width = 6, height = 5, units = "in")
  ggsave(file.path("images", paste0(base, "_epi_plot.pdf")),  plot_ee_full(dat$X, dat$Y, "epi"),
         width = 5, height = 4, units = "in")
  ggsave(file.path("images", paste0(base, "_hipo_plot.pdf")), plot_ee_full(dat$X, dat$Y, "hipo"),
         width = 5, height = 4, units = "in")
  
  # results boxplot
  if (MAKE_RESULTS) {
    res <- run_outer_cv(dat$X, dat$Y,
                        K = K_SYNTH, seed = SEED_GLOBAL,
                        models_ee = c("lda","qda","knn","svmRadial","rf"),
                        models_dd = c("lda","qda","knn"),
                        inner_number = INNER_FOLDS,
                        do_checks = TRUE,
                        allow_parallel = TRUE)
    plot_results_box(res$folds, outfile_pdf = file.path("images", sprintf("%d_results.pdf", e)))
  }
}

## =========================================================
##  B) REAL — panels + (optional) results boxplots
## =========================================================

## Berkeley Growth
message(">>> real: growth")
dat_g <- load_real("growth")
ggsave("images/7_growth_sim.pdf",
       plot_curves_two_classes(dat_g$X, dat_g$Y, xlab = "age index", ylab = "height"),
       width = 6, height = 5, units = "in")
ggsave("images/7_growth_epi_plot.pdf", plot_ee_full(dat_g$X, dat_g$Y, "epi"),
       width = 5, height = 4, units = "in")
ggsave("images/7_growth_hipo_plot.pdf", plot_ee_full(dat_g$X, dat_g$Y, "hipo"),
       width = 5, height = 4, units = "in")
if (MAKE_RESULTS) {
  res_g <- run_outer_cv(dat_g$X, dat_g$Y,
                        K = K_GROWTH, seed = SEED_GLOBAL,
                        models_ee = c("lda","qda","knn","svmRadial","rf"),
                        models_dd = c("lda","qda","knn"),
                        inner_number = INNER_FOLDS, do_checks = TRUE, allow_parallel = TRUE)
  plot_results_box(res_g$folds, outfile_pdf = "images/7_results.pdf")
}

## Tecator
message(">>> real: tecator")
dat_t <- load_real("tecator")
ggsave("images/8_tecator_sim.pdf",
       plot_curves_two_classes(dat_t$X, dat_t$Y, xlab = "wavelength index", ylab = "absorbance"),
       width = 6, height = 5, units = "in")
ggsave("images/8_tecator_epi_plot.pdf", plot_ee_full(dat_t$X, dat_t$Y, "epi"),
       width = 5, height = 4, units = "in")
ggsave("images/8_tecator_hipo_plot.pdf", plot_ee_full(dat_t$X, dat_t$Y, "hipo"),
       width = 5, height = 4, units = "in")
if (MAKE_RESULTS) {
  res_t <- run_outer_cv(dat_t$X, dat_t$Y,
                        K = K_TECATOR, seed = SEED_GLOBAL,
                        models_ee = c("lda","qda","knn","svmRadial","rf"),
                        models_dd = c("lda","qda","knn"),
                        inner_number = INNER_FOLDS, do_checks = TRUE, allow_parallel = TRUE)
  plot_results_box(res_t$folds, outfile_pdf = "images/8_results.pdf")
}

## MCO
message(">>> real: mco")
dat_m <- load_real("mco")
ggsave("images/9_mco_sim.pdf",
       plot_curves_two_classes(dat_m$X, dat_m$Y, xlab = "time index (10 s)", ylab = "mco intensity"),
       width = 6, height = 5, units = "in")
ggsave("images/9_mco_epi_plot.pdf", plot_ee_full(dat_m$X, dat_m$Y, "epi"),
       width = 5, height = 4, units = "in")
ggsave("images/9_mco_hipo_plot.pdf", plot_ee_full(dat_m$X, dat_m$Y, "hipo"),
       width = 5, height = 4, units = "in")
if (MAKE_RESULTS) {
  res_m <- run_outer_cv(dat_m$X, dat_m$Y,
                        K = K_MCO, seed = SEED_GLOBAL,
                        models_ee = c("lda","qda","knn","svmRadial","rf"),
                        models_dd = c("lda","qda","knn"),
                        inner_number = INNER_FOLDS, do_checks = TRUE, allow_parallel = TRUE)
  plot_results_box(res_m$folds, outfile_pdf = "images/9_results.pdf")
}

message("Done. All figures saved under 'images/'.")