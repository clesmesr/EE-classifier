## ============================================================
##  ARTICLE REPRO — FULL PIPELINE (ALL DATASETS)
##  - Outer K-fold CV as in the paper
##  - Validation features computed ONLY vs training folds (no leakage)
##  - Inner CV (caret) with parallelization (up to 32 cores)
##  - FM depth implemented in-script (no DepthProc)
##  - All plots saved to images/
##  - All summary tables printed to console (and LaTeX printed to console)
##  - FIX: remove all.inside=TRUE to get correct CDF tails
## ============================================================

options(stringsAsFactors = FALSE)

## ---------------------------
##  packages (auto-install)
## ---------------------------
pkgs <- c(
  "roahd","caret","fda","fda.usc","ggplot2","lattice",
  "xtable","doParallel","kernlab","randomForest","MASS"
)
to_install <- pkgs[!suppressWarnings(sapply(pkgs, requireNamespace, quietly = TRUE))]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(pkgs, require, character.only = TRUE))

## ---------------------------
##  parallel (up to 32 cores)
## ---------------------------
cores <- min(32, parallel::detectCores())
cl <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)
on.exit({ try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)

## ---------------------------
##  EE (MHI / MEI)  — fast & no leakage
##  Using per-column sorting + findInterval (vectorized)
## ---------------------------
mhi_against_ref <- function(ref, newdata) {
  n_ref <- nrow(ref); P <- ncol(ref); n_new <- nrow(newdata)
  ref_sorted <- apply(ref, 2, sort)
  out <- numeric(n_new)
  for (i in 1:n_new) {
    ## fix: do NOT clamp outside range; keep exact tail counts
    cnt_le <- vapply(seq_len(P),
                     function(t) findInterval(newdata[i, t], ref_sorted[, t]),
                     numeric(1))
    out[i] <- mean(cnt_le / n_ref)
  }
  out
}
mei_against_ref <- function(ref, newdata) {
  ## ref >= x  <=>  -ref <= -x
  n_ref <- nrow(ref); P <- ncol(ref); n_new <- nrow(newdata)
  ref_sorted_rev <- apply(-ref, 2, sort)
  out <- numeric(n_new)
  for (i in 1:n_new) {
    ## fix: do NOT clamp outside range; keep exact tail counts
    cnt_ge <- vapply(seq_len(P),
                     function(t) findInterval(-newdata[i, t], ref_sorted_rev[, t]),
                     numeric(1))
    out[i] <- mean(cnt_ge / n_ref)
  }
  out
}

## ---------------------------
##  FM (Fraiman–Muniz) — implemented here
##  FM(x) = mean_t { 1 - | 0.5 - F_t(x_t) | }, with F_t empirical CDF of ref[,t]
## ---------------------------
fm_against_ref <- function(ref, newdata) {
  n_ref <- nrow(ref); P <- ncol(ref); n_new <- nrow(newdata)
  ref_sorted <- apply(ref, 2, sort)
  out <- numeric(n_new)
  for (i in 1:n_new) {
    ## fix: do NOT clamp outside range; keep exact tail CDF
    Fvals <- vapply(seq_len(P),
                    function(t) findInterval(newdata[i, t], ref_sorted[, t]) / n_ref,
                    numeric(1))
    out[i] <- mean(1 - abs(0.5 - Fvals))
  }
  out
}

## ---------------------------
##  synthetic data generator (Exp 1–6)
## ---------------------------
generate_synthetic_experiment <- function(exp_id, N = 200, P = 100, seed = 45) {
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
  ## keep original behavior; generate_gauss_fdata may return a matrix
  X <- roahd::generate_gauss_fdata(N, mu1, C)  # class 0
  Y <- roahd::generate_gauss_fdata(N, mu2, C)  # class 1
  rownames(X) <- paste0("X_", seq_len(nrow(X)))
  rownames(Y) <- paste0("Y_", seq_len(nrow(Y)))
  list(X = X, Y = Y)
}

## ---------------------------
##  real datasets (growth, tecator, mco)
## ---------------------------
load_real <- function(which = c("growth","tecator","mco")) {
  which <- match.arg(which)
  if (which == "growth") {
    data(growth) # fda.usc
    X <- t(growth$hgtf)  # girls (0)
    Y <- t(growth$hgtm)  # boys (1)
  } else if (which == "tecator") {
    data("tecator", package = "fda.usc")
    X <- tecator$absorp.fdata$data[tecator$y$Fat <  20, ]  # low fat (0)
    Y <- tecator$absorp.fdata$data[tecator$y$Fat >= 20, ]  # high fat (1)
  } else {
    data(MCO)
    X <- MCO$intact$data[MCO$classintact == 1, ]  # control (0)
    Y <- MCO$intact$data[MCO$classintact == 2, ]  # treatment (1)
  }
  rownames(X) <- paste0("X_", seq_len(nrow(X)))
  rownames(Y) <- paste0("Y_", seq_len(nrow(Y)))
  list(X = X, Y = Y)
}

## ---------------------------
##  stratified folds
## ---------------------------
make_stratified_folds <- function(n_x, n_y, K, seed = 7) {
  set.seed(seed)
  f0 <- caret::createFolds(factor(rep(0, n_x)), k = K, list = TRUE)
  f1 <- caret::createFolds(factor(rep(1, n_y)), k = K, list = TRUE)
  list(f0 = f0, f1 = f1)
}

## ---------------------------
##  inner training with explicit grids
## ---------------------------
train_with_grid <- function(method, x, y, inner_ctrl) {
  x <- as.data.frame(x); y <- factor(y)
  
  if (method == "svmRadial") {
    sig <- as.numeric(kernlab::sigest(as.matrix(x), frac = 1)[1])
    tg  <- expand.grid(sigma = sig, C = 2^seq(-2, 3))
    fit <- suppressMessages(suppressWarnings(
      caret::train(x, y, method = "svmRadial",
                   preProcess = c("center","scale"),
                   trControl = inner_ctrl, metric = "Accuracy",
                   tuneGrid = tg)
    ))
    return(fit)
  }
  
  if (method == "rf") {
    tg <- data.frame(mtry = sort(unique(pmax(1, pmin(ncol(x), c(1,2))))))
    fit <- suppressMessages(suppressWarnings(
      caret::train(x, y, method = "rf",
                   trControl = inner_ctrl, metric = "Accuracy",
                   tuneGrid = tg, ntree = 500)
    ))
    return(fit)
  }
  
  if (method == "knn") {
    max_k <- max(3, min(25, nrow(x) - 1))
    tg <- data.frame(k = unique(seq(1, max_k, by = 2)))
    fit <- suppressMessages(suppressWarnings(
      caret::train(x, y, method = "knn",
                   preProcess = c("center","scale"),
                   trControl = inner_ctrl, metric = "Accuracy",
                   tuneGrid = tg)
    ))
    return(fit)
  }
  
  ## lda / qda (no tuning)
  fit <- suppressMessages(suppressWarnings(
    caret::train(x, y, method = method,
                 preProcess = c("center","scale"),
                 trControl = inner_ctrl, metric = "Accuracy")
  ))
  fit
}

## ---------------------------
##  outer CV (core)
## ---------------------------
run_outer_cv <- function(X, Y, K, seed = 45,
                         models_ee = c("lda","qda","knn","svmRadial","rf"),
                         models_dd = c("lda","qda","knn"),
                         inner_number = 5,
                         do_checks = TRUE) {
  
  rownames(X) <- if (is.null(rownames(X))) paste0("X_", seq_len(nrow(X))) else rownames(X)
  rownames(Y) <- if (is.null(rownames(Y))) paste0("Y_", seq_len(nrow(Y))) else rownames(Y)
  
  folds <- make_stratified_folds(nrow(X), nrow(Y), K, seed)
  acc_list <- list()
  P <- ncol(X); tol <- 1/P + 1e-6
  
  for (k in seq_len(K)) {
    idx_val_X <- folds$f0[[k]]; idx_val_Y <- folds$f1[[k]]
    X_train <- X[-idx_val_X, , drop = FALSE]
    X_val   <- X[ idx_val_X, , drop = FALSE]
    Y_train <- Y[-idx_val_Y, , drop = FALSE]
    Y_val   <- Y[ idx_val_Y, , drop = FALSE]
    
    if (do_checks) {
      stopifnot(length(intersect(rownames(X_train), rownames(X_val))) == 0)
      stopifnot(length(intersect(rownames(Y_train), rownames(Y_val))) == 0)
    }
    
    ## EEh (MHI)
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
    
    ## EEe (MEI)
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
    
    ## DD (FM)
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
    
    inner_ctrl <- caret::trainControl(method = "cv", number = inner_number, allowParallel = TRUE)
    
    ## EEh
    for (m in models_ee) {
      set.seed(7)
      mname <- if (m == "svmRadial") "EEh_svm" else paste0("EEh_", m)
      fit <- train_with_grid(m, eeh_train[, 1:2], eeh_train$g, inner_ctrl)
      pred <- predict(fit, eeh_val[, 1:2])
      acc  <- mean(pred == eeh_val$g)
      acc_list[[length(acc_list) + 1]] <- data.frame(fold = k, model = mname, acc = acc)
    }
    ## EEe
    for (m in models_ee) {
      set.seed(7)
      mname <- if (m == "svmRadial") "EEe_svm" else paste0("EEe_", m)
      fit <- train_with_grid(m, eee_train[, 1:2], eee_train$g, inner_ctrl)
      pred <- predict(fit, eee_val[, 1:2])
      acc  <- mean(pred == eee_val$g)
      acc_list[[length(acc_list) + 1]] <- data.frame(fold = k, model = mname, acc = acc)
    }
    ## DD
    for (m in models_dd) {
      set.seed(7)
      mname <- paste0("DD_", m)
      fit <- train_with_grid(m, dd_train[, 1:2], dd_train$g, inner_ctrl)
      pred <- predict(fit, dd_val[, 1:2])
      acc  <- mean(pred == dd_val$g)
      acc_list[[length(acc_list) + 1]] <- data.frame(fold = k, model = mname, acc = acc)
    }
  }
  
  acc_tbl <- do.call(rbind, acc_list)
  
  summary_tbl <- do.call(rbind, lapply(split(acc_tbl$acc, acc_tbl$model), function(v) {
    c(Min = min(v, na.rm = TRUE),
      `1st Qu.` = as.numeric(quantile(v, .25, na.rm = TRUE)),
      Median = median(v, na.rm = TRUE),
      Mean = mean(v, na.rm = TRUE),
      `3rd Qu.` = as.numeric(quantile(v, .75, na.rm = TRUE)),
      Max = max(v, na.rm = TRUE),
      `NA's` = sum(is.na(v)))
  }))
  summary_tbl <- as.data.frame(summary_tbl)
  
  list(summary = summary_tbl, folds = acc_tbl)
}

## ---------------------------
##  pretty summary (console + latex)
## ---------------------------
row_order <- c("EEh_lda","EEh_qda","EEh_knn","EEh_svm","EEh_rf",
               "EEe_lda","EEe_qda","EEe_knn","EEe_svm","EEe_rf",
               "DD_lda","DD_qda","DD_knn")
col_order <- c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max","NA's")
canonicalize_cols <- function(df) {
  nm <- colnames(df)
  nm <- sub("^Min$", "Min.", nm)
  nm <- sub("^1st Qu\\..*$", "1st Qu.", nm)
  nm <- sub("^3rd Qu\\..*$", "3rd Qu.", nm)
  nm <- sub("^NA.?s$", "NA's", nm, ignore.case = TRUE)
  colnames(df) <- nm
  df
}
pretty_summary <- function(summary_tbl) {
  out <- canonicalize_cols(summary_tbl)
  keep_rows <- intersect(row_order, rownames(out))
  out <- out[keep_rows, , drop = FALSE]
  keep_cols <- intersect(col_order, colnames(out))
  out <- out[, keep_cols, drop = FALSE]
  out
}
print_latex_table <- function(summary_tbl, caption, label = NULL, digits = 4) {
  out <- pretty_summary(summary_tbl)
  df  <- data.frame(method = tolower(rownames(out)), out, row.names = NULL, check.names = FALSE)
  xt  <- xtable::xtable(df, caption = caption, label = label)
  ## dynamic alignment: 1 (overall) + ncol(df)
  xtable::align(xt) <- c("l", rep("c", ncol(df)))
  print(xt,
        include.rownames = FALSE,
        booktabs = TRUE,
        sanitize.text.function = identity,
        caption.placement = "top",
        digits = c(0, 0, rep(digits, ncol(df) - 1)))
}

## ---------------------------
##  plotting helpers (lower-case, no title)
## ---------------------------
if (!dir.exists("images")) dir.create("images", recursive = TRUE)

plot_curves_two_classes <- function(X, Y, xlab = "time index", ylab = "value",
                                    legend_labels = c("class 0","class 1")) {
  P <- ncol(X); tgrid <- seq_len(P)
  to_df <- function(M, g_label) data.frame(t = rep(tgrid, each = nrow(M)),
                                           y = as.vector(t(M)), g = g_label)
  df <- rbind(to_df(X, "0"), to_df(Y, "1"))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = t, y = y, group = interaction(g, t), color = g)) +
    ggplot2::geom_line(alpha = 0.6) +
    ggplot2::scale_color_manual(values = c("0" = "black", "1" = "grey50"),
                                labels = legend_labels, name = "group") +
    ggplot2::labs(x = xlab, y = ylab) +  # no title
    ggplot2::theme_minimal()
  p
}

plot_ee_full <- function(X, Y, index = c("epi", "hipo")) {
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
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, shape = g)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_shape_manual(values = c("0" = 4, "1" = 1), name = "group",
                                labels = c("class 0","class 1")) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = xlab, y = ylab) +  # no title
    ggplot2::theme_minimal()
  p
}

plot_results_box <- function(folds_df, outfile_pdf) {
  # lower-case method labels
  folds_df$model <- factor(tolower(as.character(folds_df$model)),
                           levels = tolower(c("EEh_lda","EEh_qda","EEh_knn","EEh_svm","EEh_rf",
                                              "EEe_lda","EEe_qda","EEe_knn","EEe_svm","EEe_rf",
                                              "DD_lda","DD_qda","DD_knn")))
  p <- ggplot2::ggplot(folds_df, ggplot2::aes(x = model, y = acc)) +
    ggplot2::geom_boxplot(outlier.size = 0.8) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "accuracy") +  # no title
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(filename = outfile_pdf, plot = p, width = 6, height = 5, units = "in")
  invisible(p)
}

## ---------------------------
##  RUN: Synthetic (1–6) and Real (growth, tecator, mco)
## ---------------------------
set.seed(45)

## ========== Synthetic Experiments 1–6 (K=100) ==========
for (e in 1:6) {
  cat(sprintf("\n=== synthetic experiment %d (100-fold) ===\n", e))
  dat_e <- generate_synthetic_experiment(e, N = 200, P = 100, seed = 45)
  res_e <- run_outer_cv(dat_e$X, dat_e$Y,
                        K = 100, seed = 45,
                        models_ee = c("lda","qda","knn","svmRadial","rf"),
                        models_dd = c("lda","qda","knn"),
                        inner_number = 5, do_checks = TRUE)
  # console table
  pretty <- pretty_summary(res_e$summary)
  print(round(pretty, 6))
  # LaTeX table (printed to console)
  print_latex_table(res_e$summary,
                    caption = sprintf("results experiment %d", e),
                    label   = sprintf("tab:exp%d", e),
                    digits  = 4)
  # boxplot to images/<e>_results.pdf
  plot_results_box(res_e$folds, outfile_pdf = file.path("images", sprintf("%d_results.pdf", e)))
}

## ========== Real: Berkeley Growth (K=20) ==========
cat("\n=== berkeley growth (20-fold) ===\n")
dat_g <- load_real("growth")
res_g <- run_outer_cv(dat_g$X, dat_g$Y,
                      K = 20, seed = 45,
                      models_ee = c("lda","qda","knn","svmRadial","rf"),
                      models_dd = c("lda","qda","knn"),
                      inner_number = 5, do_checks = TRUE)
print(round(pretty_summary(res_g$summary), 6))
print_latex_table(res_g$summary,
                  caption = "results berkeley growth (20-fold cv)",
                  label   = "tab:growth", digits = 4)
# figures (lower-case; no title)
p_curves <- plot_curves_two_classes(dat_g$X, dat_g$Y,
                                    xlab = "age index", ylab = "height",
                                    legend_labels = c("class 0","class 1"))
ggplot2::ggsave("images/7_growth_sim.pdf", p_curves, width = 6, height = 5, units = "in")
ggplot2::ggsave("images/7_growth_epi_plot.pdf", plot_ee_full(dat_g$X, dat_g$Y, "epi"),
                width = 5, height = 4, units = "in")
ggplot2::ggsave("images/7_growth_hipo_plot.pdf", plot_ee_full(dat_g$X, dat_g$Y, "hipo"),
                width = 5, height = 4, units = "in")
plot_results_box(res_g$folds, outfile_pdf = "images/7_results.pdf")

## ========== Real: Tecator (K=50) ==========
cat("\n=== tecator (50-fold) ===\n")
dat_t <- load_real("tecator")
res_t <- run_outer_cv(dat_t$X, dat_t$Y,
                      K = 50, seed = 45,
                      models_ee = c("lda","qda","knn","svmRadial","rf"),
                      models_dd = c("lda","qda","knn"),
                      inner_number = 5, do_checks = TRUE)
print(round(pretty_summary(res_t$summary), 6))
print_latex_table(res_t$summary,
                  caption = "results tecator (50-fold cv)",
                  label   = "tab:tecator", digits = 4)
# figures
p_curves <- plot_curves_two_classes(dat_t$X, dat_t$Y,
                                    xlab = "wavelength index", ylab = "absorbance",
                                    legend_labels = c("class 0","class 1"))
ggplot2::ggsave("images/8_tecator_sim.pdf", p_curves, width = 6, height = 5, units = "in")
ggplot2::ggsave("images/8_tecator_epi_plot.pdf", plot_ee_full(dat_t$X, dat_t$Y, "epi"),
                width = 5, height = 4, units = "in")
ggplot2::ggsave("images/8_tecator_hipo_plot.pdf", plot_ee_full(dat_t$X, dat_t$Y, "hipo"),
                width = 5, height = 4, units = "in")
plot_results_box(res_t$folds, outfile_pdf = "images/8_results.pdf")

## ========== Real: MCO (K=30) ==========
cat("\n=== mco (30-fold) ===\n")
dat_m <- load_real("mco")
res_m <- run_outer_cv(dat_m$X, dat_m$Y,
                      K = 30, seed = 45,
                      models_ee = c("lda","qda","knn","svmRadial","rf"),
                      models_dd = c("lda","qda","knn"),
                      inner_number = 5, do_checks = TRUE)
print(round(pretty_summary(res_m$summary), 6))
print_latex_table(res_m$summary,
                  caption = "results mco (30-fold cv)",
                  label   = "tab:mco", digits = 4)
# figures
p_curves <- plot_curves_two_classes(dat_m$X, dat_m$Y,
                                    xlab = "time index (10 s)", ylab = "mco intensity",
                                    legend_labels = c("class 0","class 1"))
ggplot2::ggsave("images/9_mco_sim.pdf", p_curves, width = 6, height = 5, units = "in")
ggplot2::ggsave("images/9_mco_epi_plot.pdf", plot_ee_full(dat_m$X, dat_m$Y, "epi"),
                width = 5, height = 4, units = "in")
ggplot2::ggsave("images/9_mco_hipo_plot.pdf", plot_ee_full(dat_m$X, dat_m$Y, "hipo"),
                width = 5, height = 4, units = "in")
plot_results_box(res_m$folds, outfile_pdf = "images/9_results.pdf")



