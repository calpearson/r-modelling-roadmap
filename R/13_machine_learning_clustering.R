# ============================================================
# explain_rf.R — plain-language interpretation of a
# randomForest::randomForest() model, organised around four questions:
#
#   1. How good is the model?                (OOB error, AUC, Brier, calibration)
#   2. Which variables are most important?   (permutation + Gini importance)
#   3. How well does it classify patients?   (confusion matrix, sens/spec/PPV/NPV)
#   4. How do predicted risks differ across patients?
#                                            (risk distribution, deciles, subgroups,
#                                             partial dependence)
#
# NOTE: this REPLACES the earlier explain_rf() in explain_trees.R (same
# name). If you source both files, whichever is sourced LAST wins — source
# this one after explain_trees.R, or delete the old explain_rf() from there.
#
# EVERY performance number here uses OUT-OF-BAG (OOB) predictions — each
# patient is predicted only by the trees that never saw them during
# training. This is an honest internal validation (same idea as
# predict(out_rf, type = "prob") in your own code), but it is still
# internal: it does not replace validation on new/external data.
#
# Text output (always printed): forest settings, OOB error, AUC (+CI),
# Brier score, calibration slope/intercept, top variables, confusion
# matrix + classification metrics, predicted-risk summary, caveats.
#
# Optional tables/figures (all FALSE/NULL by default):
#   importance_table       — full importance table (permutation + Gini, with CIs)
#   importance_plot        — importance: point-range (permutation) + bars (Gini)
#   confusion_plot         — confusion-matrix heatmap
#   roc_plot               — ROC curve with AUC, Youden point, default-threshold point
#   risk_distribution_plot — histogram of predicted risk by actual outcome
#   decile_table           — observed vs predicted event rate by risk decile
#   calibration_plot       — observed vs predicted risk by decile (with 95% CIs)
#   gains_plot             — cumulative gains: % of events captured vs % screened
#   error_plot             — OOB error vs number of trees (has the forest converged?)
#   pdp_vars               — partial dependence plot(s): a character vector of
#                            variable names, or "auto" for the top 3 by importance
#   pdp_ice                — TRUE adds individual (ICE) lines to each PDP
#   risk_by                — character vector of variable names: table + plot of
#                            mean predicted vs observed risk by subgroup
#   all_outputs            — TRUE switches on everything above (pdp "auto",
#                            risk_by = the factor predictors)
#
# `data` is the training data frame (needed only for pdp_vars / risk_by).
# If omitted, the function tries to find it from the model's call.
#
# REQUIRED: randomForest. OPTIONAL: ggplot2 (all figures), pROC (DeLong AUC
# CI — a Hanley-McNeil approximation is used if pROC isn't installed),
# crayon (bold text).
#   install.packages(c("randomForest", "ggplot2", "pROC"))
# ============================================================

library(randomForest)

# ---- Shared helpers ----
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) crayon::bold(x) else x
}
format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}

.rf_pct <- function(x, d = 1) {
  ifelse(is.finite(x), sprintf(paste0("%.", d, "f%%"), 100 * x), "NA")
}

.rf_need_ggplot <- function(what) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning(sprintf("%s requires the 'ggplot2' package. Run install.packages('ggplot2').", what), call. = FALSE)
    return(FALSE)
  }
  TRUE
}

.rf_print_plot <- function(p) {
  # explicit print: ggplot objects only auto-draw at top level, not inside functions
  if (!is.null(p)) print(p)
  invisible(p)
}

# ---- ROC curve + AUC computed directly (no pROC dependency) ----
# Handles tied predicted values correctly (ties become diagonal segments).
.rf_roc <- function(is_pos, prob) {
  ok <- !is.na(is_pos) & !is.na(prob)
  is_pos <- is_pos[ok]; prob <- prob[ok]
  n_pos <- sum(is_pos); n_neg <- sum(!is_pos)
  if (n_pos == 0 || n_neg == 0) return(NULL)
  o <- order(prob, decreasing = TRUE)
  p <- prob[o]; y <- is_pos[o]
  tps <- cumsum(y); fps <- cumsum(!y)
  last_of_tie <- c(diff(p) != 0, TRUE)
  tpr <- c(0, tps[last_of_tie] / n_pos)
  fpr <- c(0, fps[last_of_tie] / n_neg)
  thr <- c(Inf, p[last_of_tie])
  auc <- sum(diff(fpr) * (utils::head(tpr, -1) + utils::tail(tpr, -1)) / 2)
  list(fpr = fpr, tpr = tpr, thr = thr, auc = auc, n_pos = n_pos, n_neg = n_neg)
}

# Hanley-McNeil approximate CI for an AUC (fallback when pROC is unavailable)
.rf_auc_ci_hm <- function(auc, n_pos, n_neg, level = 0.95) {
  q1 <- auc / (2 - auc); q2 <- 2 * auc^2 / (1 + auc)
  se <- sqrt((auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) + (n_neg - 1) * (q2 - auc^2)) / (n_pos * n_neg))
  z <- stats::qnorm(1 - (1 - level) / 2)
  c(lower = max(0, auc - z * se), upper = min(1, auc + z * se))
}

# Wilson score interval for a proportion (vectorised)
.rf_wilson <- function(x, n, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  data.frame(lower = pmax(0, centre - half), upper = pmin(1, centre + half))
}

# Same bucket rule as dplyr::ntile() (ties broken by row order)
.rf_ntile <- function(x, n) {
  len <- length(x)
  as.integer(floor(n * (rank(x, ties.method = "first") - 1) / len + 1))
}

# Binary classification metrics from logical vectors
.rf_class_metrics <- function(is_pos, pred_pos) {
  tp <- sum(is_pos & pred_pos); fp <- sum(!is_pos & pred_pos)
  fn <- sum(is_pos & !pred_pos); tn <- sum(!is_pos & !pred_pos)
  n <- tp + fp + fn + tn
  sens <- tp / (tp + fn); spec <- tn / (tn + fp)
  ppv <- tp / (tp + fp);  npv <- tn / (tn + fn)
  acc <- (tp + tn) / n
  pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n^2
  kappa <- (acc - pe) / (1 - pe)
  list(tp = tp, fp = fp, fn = fn, tn = tn, n = n,
       sensitivity = sens, specificity = spec, ppv = ppv, npv = npv,
       accuracy = acc, balanced_accuracy = (sens + spec) / 2, kappa = kappa)
}

# Group a variable for subgroup tables (factors as-is; numerics into quartiles)
.rf_group <- function(x) {
  if (is.factor(x) || is.character(x) || is.logical(x)) return(factor(x))
  ux <- unique(x[!is.na(x)])
  if (length(ux) <= 6) return(factor(x))
  br <- unique(stats::quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE))
  if (length(br) < 3) return(factor(x))
  cut(x, breaks = br, include.lowest = TRUE, dig.lab = 4)
}

# Prediction for the "focus" outcome: P(positive class) for classification,
# the predicted value for regression
.rf_focus_pred <- function(model, nd, is_class, focus) {
  if (is_class) as.numeric(stats::predict(model, newdata = nd, type = "prob")[, focus])
  else as.numeric(stats::predict(model, newdata = nd))
}

# Partial dependence (and optional ICE) computed directly: set the variable
# to each grid value for every row, predict, average. Rows are subsampled
# deterministically (evenly spaced) if there are many — no RNG used.
.rf_pdp <- function(model, data, var, is_class, focus, ice = FALSE,
                    max_rows = 500, n_grid = 20, n_ice = 50) {
  x <- data[[var]]
  idx <- if (nrow(data) > max_rows) unique(round(seq(1, nrow(data), length.out = max_rows))) else seq_len(nrow(data))
  d <- data[idx, , drop = FALSE]
  
  is_cat <- is.factor(x) || is.character(x)
  if (is_cat) {
    lv <- if (is.factor(x)) levels(x) else sort(unique(x))
    grid <- lv
  } else {
    ux <- sort(unique(x[!is.na(x)]))
    grid <- if (length(ux) <= 10) ux else unique(unname(stats::quantile(x, probs = seq(0.05, 0.95, length.out = n_grid), na.rm = TRUE)))
  }
  
  ice_rows <- unique(round(seq(1, nrow(d), length.out = min(n_ice, nrow(d)))))
  mean_vals <- numeric(length(grid))
  ice_list <- vector("list", length(grid))
  
  for (k in seq_along(grid)) {
    nd <- d
    nd[[var]] <- if (is_cat) factor(grid[k], levels = lv) else grid[k]
    pr <- .rf_focus_pred(model, nd, is_class, focus)
    mean_vals[k] <- mean(pr, na.rm = TRUE)
    if (ice) ice_list[[k]] <- data.frame(x_value = grid[k], yhat = pr[ice_rows], row_id = ice_rows)
  }
  list(pdp = data.frame(x_value = if (is_cat) factor(grid, levels = lv) else grid, yhat = mean_vals),
       ice = if (ice) do.call(rbind, ice_list) else NULL,
       is_cat = is_cat, n_rows_used = nrow(d))
}


# ============================================================
# explain_rf()
# ============================================================
explain_rf <- function(model, data = NULL, positive_class = NULL, conf_level = 0.95,
                       importance_table = FALSE, importance_plot = FALSE,
                       confusion_plot = FALSE, roc_plot = FALSE,
                       risk_distribution_plot = FALSE, decile_table = FALSE,
                       calibration_plot = FALSE, gains_plot = FALSE,
                       error_plot = FALSE, pdp_vars = NULL, pdp_ice = FALSE,
                       risk_by = NULL, all_outputs = FALSE) {
  
  # ==========================================================
  # 0. VALIDATION + SETUP
  # ==========================================================
  if (!inherits(model, "randomForest")) {
    stop("explain_rf() expects a model fit with randomForest::randomForest().")
  }
  model_type <- tryCatch({
    t <- model$type
    if (is.null(t)) NA_character_ else t
  }, error = function(e) NA_character_)
  if (is.na(model_type) || !(model_type %in% c("classification", "regression"))) {
    stop("explain_rf() supports classification and regression forests only (not unsupervised).")
  }
  is_class <- model_type == "classification"
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  alpha_pct <- conf_level * 100
  
  y_obs <- model$y
  if (is.null(y_obs)) {
    stop("This forest does not contain its response values (model$y is NULL), so OOB performance cannot be evaluated.")
  }
  n_model <- length(y_obs)
  ntree <- model$ntree
  outcome <- tryCatch(deparse(model$terms[[2]]), error = function(e) "the outcome")
  
  # ---- Locate training data (only needed for pdp_vars / risk_by) ----
  data_note <- NULL
  if (is.null(data)) {
    data <- tryCatch(eval(model$call$data, envir = parent.frame()), error = function(e) NULL)
    if (is.data.frame(data)) data_note <- "(training data located automatically from the model's call)"
  }
  if (!is.null(data) && !is.data.frame(data)) data <- NULL
  if (!is.null(data) && nrow(data) != n_model) {
    data_note <- sprintf("(a data frame was found but has %d rows vs. %d in the model — ignored; pass data = <your training data> explicitly)",
                         nrow(data), n_model)
    data <- NULL
  }
  
  predictors <- tryCatch({
    pn <- rownames(model$importance)
    if (is.null(pn)) attr(model$terms, "term.labels") else pn
  }, error = function(e) character(0))
  
  # ---- all_outputs switches everything on ----
  if (isTRUE(all_outputs)) {
    importance_table <- importance_plot <- confusion_plot <- roc_plot <- TRUE
    risk_distribution_plot <- decile_table <- calibration_plot <- gains_plot <- error_plot <- TRUE
    if (is.null(pdp_vars)) pdp_vars <- "auto"
    if (is.null(risk_by) && !is.null(data)) {
      fac_preds <- predictors[predictors %in% names(data)]
      fac_preds <- fac_preds[vapply(fac_preds, function(v) is.factor(data[[v]]) || is.character(data[[v]]), logical(1))]
      if (length(fac_preds) > 0) risk_by <- fac_preds
    }
  }
  
  # ---- Classes / positive class ----
  classes <- NULL; pos <- NULL; neg <- NULL; is_binary <- FALSE
  if (is_class) {
    if (!is.factor(y_obs)) stop("Classification forest, but model$y is not a factor — unexpected object structure.")
    classes <- levels(y_obs)
    is_binary <- length(classes) == 2
    if (is_binary) {
      pos <- if (!is.null(positive_class)) positive_class else if ("Yes" %in% classes) "Yes" else classes[2]
      if (!(pos %in% classes)) stop(sprintf("positive_class = '%s' is not one of this model's outcome levels: %s", pos, paste(classes, collapse = ", ")))
      neg <- setdiff(classes, pos)
    }
  }
  focus <- if (is_class) { if (!is.null(positive_class) && positive_class %in% classes) positive_class else if (is_binary) pos else classes[1] } else NULL
  
  # ==========================================================
  # 1. OOB PROBABILITIES, TRUTH, DIAGNOSTICS
  # ==========================================================
  prob_pos <- NULL; is_pos <- NULL; ok <- NULL; y01 <- NULL
  if (is_class && is_binary) {
    votes <- tryCatch(as.matrix(unclass(model$votes)), error = function(e) NULL)
    if (is.null(votes) || !(pos %in% colnames(votes))) {
      stop("Could not read the OOB vote matrix (model$votes) for this forest — it is needed for AUC/calibration.")
    }
    rs <- rowSums(votes)
    prob_mat <- votes / rs
    prob_mat[rs == 0, ] <- NA   # patients never out-of-bag (very rare) have no OOB prediction
    prob_pos <- prob_mat[, pos]
    is_pos <- (y_obs == pos)
    ok <- !is.na(prob_pos) & !is.na(is_pos)
    y01 <- as.integer(is_pos)
  }
  
  diagnostics <- character(0)
  
  if (is_class && is_binary) {
    if (sum(!ok) > 0) {
      diagnostics <- c(diagnostics, sprintf(
        "has %d patient(s) with no OOB prediction (never left out of any tree) — excluded from the performance figures below.", sum(!ok)))
    }
    n_ev <- sum(is_pos[ok]); n_ne <- sum(!is_pos[ok])
    if (min(n_ev, n_ne) < 50) {
      diagnostics <- c(diagnostics, sprintf(
        "has only %d patients in its smaller outcome group. With so few, every performance figure below (AUC, sensitivity, calibration) is statistically noisy — read the confidence intervals, not just the point estimates.",
        min(n_ev, n_ne)))
    }
    if (max(n_ev, n_ne) / max(1, min(n_ev, n_ne)) > 9) {
      diagnostics <- c(diagnostics,
                       "has strongly IMBALANCED outcome classes (>9:1). Overall accuracy can look excellent just by predicting the majority class — rely on AUC, sensitivity/specificity and calibration instead, and consider randomForest(..., sampsize/strata/classwt) or a lower classification threshold (see the Youden threshold below).")
    }
  }
  
  err_curve <- tryCatch({
    if (is_class) model$err.rate[, "OOB"] else model$mse
  }, error = function(e) NULL)
  if (!is.null(err_curve) && length(err_curve) >= 20) {
    last_10pct <- utils::tail(err_curve, ceiling(length(err_curve) * 0.1))
    mid_ref <- err_curve[ceiling(length(err_curve) * 0.5)]
    if (is.finite(mid_ref) && mid_ref != 0) {
      rel_drop <- (mid_ref - mean(last_10pct)) / mid_ref
      if (is.finite(rel_drop) && rel_drop > 0.02) {
        diagnostics <- c(diagnostics, sprintf(
          "may benefit from MORE trees — error was still improving by ~%.1f%% between the midpoint and the end of the %d trees. Refit with a larger ntree and check the curve flattens (error_plot = TRUE).",
          rel_drop * 100, ntree))
      }
    }
  }
  
  has_importance <- !is.null(model$importance) && length(model$importance) > 0
  if (!has_importance) {
    diagnostics <- c(diagnostics, "contains no variable importance values. Refit with randomForest(..., importance = TRUE).")
  }
  
  if (length(diagnostics) > 0) {
    cat("########################################################\n")
    cat("# MODEL DIAGNOSTIC WARNINGS — read before trusting results\n")
    cat("########################################################\n")
    for (d in diagnostics) cat(bold("WARNING:"), "this model", d, "\n\n")
  }
  
  # ==========================================================
  # 2. HEADER
  # ==========================================================
  cat("========================================================\n")
  cat("RANDOM FOREST —", ifelse(is_class, "classification", "regression"), "\n")
  cat("MODEL:", tryCatch(deparse(model$call$formula), error = function(e) "(formula not available)"), "\n")
  cat("========================================================\n\n")
  
  cat("What a random forest does, in plain English:\n")
  cat("It grows hundreds of decision trees, each on a random resample of patients and\n")
  cat("a random subset of variables at every split, then lets them vote. The crowd is\n")
  cat("far more stable than any single tree, but you can no longer read off simple\n")
  cat("rules — you judge it by how well it PREDICTS, not by coefficients.\n\n")
  
  mtry <- tryCatch(model$mtry, error = function(e) NA)
  cat(sprintf("Forest: %s trees, %s variables tried at each split (mtry), %s patients.\n",
              bold(ntree), bold(ifelse(is.null(mtry), "?", mtry)), bold(format(n_model, big.mark = ","))))
  if (is_class && is_binary) {
    cat(sprintf("Outcome '%s': EVENT = '%s' (%s patients, %s of the cohort); non-event = '%s' (%s).\n",
                outcome, pos, format(sum(is_pos[ok]), big.mark = ","), .rf_pct(mean(is_pos[ok])),
                neg, format(sum(!is_pos[ok]), big.mark = ",")))
    cat("(Change which level counts as the event with positive_class = \"...\")\n")
  } else if (is_class) {
    cat(sprintf("Outcome '%s' has %d classes: %s.\n", outcome, length(classes), paste(classes, collapse = ", ")))
  }
  cat("\nEvery performance figure below uses OUT-OF-BAG (OOB) predictions: each patient\n")
  cat("is scored only by the ~37% of trees that never saw them. That is an honest\n")
  cat("internal check without a separate test set — but it is still the SAME dataset,\n")
  cat("so it can't replace validation on genuinely new patients.\n")
  if (!is.null(data_note)) cat(data_note, "\n")
  cat("\n")
  
  # Containers for the returned object
  out <- list()
  
  # ==========================================================
  # Q1. HOW GOOD IS THE MODEL?
  # ==========================================================
  cat("========================================================\n")
  cat("Q1. HOW GOOD IS THE MODEL?\n")
  cat("========================================================\n\n")
  
  roc_obj <- NULL; youden_thr <- NA_real_
  if (is_class && is_binary) {
    pp <- prob_pos[ok]; ip <- is_pos[ok]; yy <- y01[ok]
    n_ok <- length(pp)
    
    # ---- OOB error + no-information-rate test ----
    pred_default <- (model$predicted == pos)[ok]
    n_correct <- sum(pred_default == ip)
    acc <- n_correct / n_ok
    nir <- max(mean(ip), 1 - mean(ip))
    nir_class <- if (mean(ip) >= 0.5) pos else neg
    nir_test <- tryCatch(stats::binom.test(n_correct, n_ok, p = nir, alternative = "greater"), error = function(e) NULL)
    
    cat(sprintf("OOB accuracy: %s (%s of %s patients classified correctly at the default rule\n",
                bold(.rf_pct(acc)), format(n_correct, big.mark = ","), format(n_ok, big.mark = ",")))
    cat("'flag as event if at least half the trees vote event').\n")
    cat(sprintf("A model that always guessed '%s' would score %s (the no-information rate).\n", nir_class, .rf_pct(nir)))
    if (!is.null(nir_test)) {
      cat(sprintf("One-sided binomial test, accuracy > no-information rate: p %s.\n", format_pval(nir_test$p.value)))
      if (nir_test$p.value >= 0.05) {
        cat(bold("CAUTION:"), "accuracy is NOT clearly better than always guessing the majority class —\n")
        cat("accuracy alone is a misleading yardstick here; use the AUC and the sensitivity/\n")
        cat("specificity below.\n")
      }
    }
    cat("\n")
    
    # ---- ROC / AUC ----
    roc_obj <- .rf_roc(ip, pp)
    if (!is.null(roc_obj)) {
      auc <- roc_obj$auc
      auc_ci <- NULL; ci_method <- NULL
      if (requireNamespace("pROC", quietly = TRUE)) {
        auc_ci <- tryCatch({
          r <- pROC::roc(response = yy, predictor = pp, levels = c(0, 1), direction = "<", quiet = TRUE)
          ci <- as.numeric(pROC::ci.auc(r, conf.level = conf_level, method = "delong"))
          c(lower = ci[1], upper = ci[3])
        }, error = function(e) NULL)
        if (!is.null(auc_ci)) ci_method <- "DeLong"
      }
      if (is.null(auc_ci)) {
        auc_ci <- .rf_auc_ci_hm(auc, roc_obj$n_pos, roc_obj$n_neg, conf_level)
        ci_method <- "Hanley-McNeil approximation"
      }
      auc_p <- tryCatch(stats::wilcox.test(pp[ip], pp[!ip], alternative = "greater", exact = FALSE)$p.value,
                        error = function(e) NA)
      
      band <- if (auc < 0.6) "essentially no better than chance"
      else if (auc < 0.7) "poor discrimination"
      else if (auc < 0.8) "acceptable discrimination"
      else if (auc < 0.9) "excellent discrimination"
      else "outstanding discrimination (check for leakage if this is unexpected)"
      
      cat(sprintf("AUC (area under the ROC curve): %s (%.0f%% CI %.3f to %.3f, %s)\n",
                  bold(sprintf("%.3f", auc)), alpha_pct, auc_ci[["lower"]], auc_ci[["upper"]], ci_method))
      cat(sprintf("Test of AUC > 0.5 (Mann-Whitney): p %s.\n", format_pval(auc_p)))
      cat(sprintf("In plain English: if you pick one patient who had the event and one who didn't,\n"))
      cat(sprintf("the model gives the event patient the HIGHER predicted risk %.0f%% of the time.\n", auc * 100))
      cat(sprintf("0.5 = a coin flip, 1.0 = perfect ranking. %s is %s.\n\n", sprintf("%.2f", auc), band))
      
      j <- roc_obj$tpr - roc_obj$fpr
      youden_thr <- roc_obj$thr[which.max(j)]
      out$auc <- auc; out$auc_ci <- auc_ci
    }
    
    # ---- Brier score + calibration ----
    brier <- mean((pp - yy)^2)
    brier_ref <- mean(yy) * (1 - mean(yy))
    scaled_brier <- if (brier_ref > 0) 1 - brier / brier_ref else NA
    cat(sprintf("Brier score: %.4f (a model that predicted everyone's risk as the overall average\n", brier))
    cat(sprintf("would score %.4f; lower is better). Scaled Brier = %.1f%% improvement over that.\n\n",
                brier_ref, 100 * scaled_brier))
    
    p_clip <- pmin(pmax(pp, 0.01), 0.99)
    lp <- stats::qlogis(p_clip)
    slope_fit <- tryCatch(suppressWarnings(stats::glm(yy ~ lp, family = stats::binomial())), error = function(e) NULL)
    citl_fit  <- tryCatch(suppressWarnings(stats::glm(yy ~ 1 + offset(lp), family = stats::binomial())), error = function(e) NULL)
    cal_slope <- NA; cal_slope_ci <- c(NA, NA); cal_citl <- NA
    if (!is.null(slope_fit)) {
      cal_slope <- unname(stats::coef(slope_fit)[2])
      cal_slope_ci <- tryCatch(unname(stats::confint.default(slope_fit, level = conf_level)[2, ]), error = function(e) c(NA, NA))
    }
    if (!is.null(citl_fit)) cal_citl <- unname(stats::coef(citl_fit)[1])
    
    cat("Calibration (do predicted risks match observed risks?):\n")
    if (is.finite(cal_citl)) {
      cat(sprintf("  Calibration-in-the-large (intercept): %.2f — ", cal_citl))
      cat(if (abs(cal_citl) < 0.1) "average predicted risk is about right.\n"
          else if (cal_citl > 0) "risk is UNDER-predicted on average (observed > predicted).\n"
          else "risk is OVER-predicted on average (observed < predicted).\n")
    }
    if (is.finite(cal_slope)) {
      cat(sprintf("  Calibration slope: %.2f (%.0f%% CI %.2f to %.2f) — ", cal_slope, alpha_pct,
                  cal_slope_ci[1], cal_slope_ci[2]))
      cat(if (cal_slope > 1.1) "predictions are not extreme enough (compressed toward the middle), typical of random forests.\n"
          else if (cal_slope < 0.9) "predictions are too extreme (a sign of over-fitting).\n"
          else "spread of predictions is about right.\n")
      cat("  (Ideal: intercept 0, slope 1. Forest probabilities are coarse and often exactly 0 or 1,\n")
      cat("  so these two numbers are approximate — judge them alongside the calibration plot.)\n")
    }
    cat("\n")
    out$brier <- brier; out$calibration <- list(slope = cal_slope, slope_ci = cal_slope_ci, intercept = cal_citl)
    
  } else if (is_class) {
    cat(sprintf("Overall OOB error rate: %s (%s correct).\n\n", bold(.rf_pct(tail(model$err.rate[, "OOB"], 1))),
                .rf_pct(1 - tail(model$err.rate[, "OOB"], 1))))
    cat("(AUC, Brier score and calibration are only reported for binary outcomes.)\n\n")
  } else {
    mse_last <- tail(model$mse, 1); rsq_last <- tail(model$rsq, 1)
    cat(sprintf("OOB mean squared error: %.4f (RMSE %.4f, in the outcome's own units).\n", mse_last, sqrt(mse_last)))
    cat(sprintf("OOB percent variance explained: %s — a random-forest version of R-squared.\n\n", bold(.rf_pct(rsq_last))))
    out$oob_mse <- mse_last; out$oob_rsq <- rsq_last
  }
  
  # ==========================================================
  # Q2. WHICH VARIABLES ARE MOST IMPORTANT?
  # ==========================================================
  cat("========================================================\n")
  cat("Q2. WHICH VARIABLES ARE MOST IMPORTANT?\n")
  cat("========================================================\n\n")
  
  imp_df <- NULL
  if (has_importance) {
    imp_raw <- as.matrix(model$importance)
    perm_col <- intersect(c("MeanDecreaseAccuracy", "%IncMSE"), colnames(imp_raw))
    node_col <- intersect(c("MeanDecreaseGini", "IncNodePurity"), colnames(imp_raw))
    perm_col <- if (length(perm_col) > 0) perm_col[1] else NA_character_
    node_col <- if (length(node_col) > 0) node_col[1] else NA_character_
    
    vars <- rownames(imp_raw)
    perm <- if (!is.na(perm_col)) imp_raw[, perm_col] else rep(NA_real_, length(vars))
    node <- if (!is.na(node_col)) imp_raw[, node_col] else rep(NA_real_, length(vars))
    
    # importanceSD holds randomForest's own "standard errors" of the permutation
    # importance (it is the denominator importance(scale = TRUE) divides by)
    perm_se <- rep(NA_real_, length(vars))
    sd_raw <- tryCatch(model$importanceSD, error = function(e) NULL)
    if (!is.na(perm_col) && !is.null(sd_raw)) {
      if (is.null(dim(sd_raw))) {
        if (!is.null(names(sd_raw))) perm_se <- as.numeric(sd_raw[vars])
      } else if (perm_col %in% colnames(sd_raw)) {
        perm_se <- as.numeric(sd_raw[vars, perm_col])
      }
    }
    
    imp_df <- data.frame(variable = vars, permutation = as.numeric(perm), perm_se = perm_se,
                         perm_lower = as.numeric(perm) - z * perm_se,
                         perm_upper = as.numeric(perm) + z * perm_se,
                         gini = as.numeric(node), row.names = NULL, stringsAsFactors = FALSE)
    sort_key <- if (!is.na(perm_col)) imp_df$permutation else imp_df$gini
    imp_df <- imp_df[order(-sort_key), , drop = FALSE]
    rownames(imp_df) <- NULL
    imp_df$distinguishable_from_zero <- ifelse(is.finite(imp_df$perm_lower), imp_df$perm_lower > 0, NA)
    
    perm_name <- if (is_class) "MeanDecreaseAccuracy" else "%IncMSE"
    node_name <- if (is_class) "MeanDecreaseGini" else "IncNodePurity"
    
    if (!is.na(perm_col)) {
      cat(sprintf("%s (permutation importance): how much WORSE the model's predictions get if\n", perm_name))
      cat("this variable's values are randomly shuffled, breaking its link to the outcome.\n")
      cat("Bigger = the model leans on it more. This is the measure to trust most.\n\n")
    } else {
      cat(bold("NOTE:"), "no permutation importance was stored (refit with importance = TRUE);\n")
      cat("only the Gini measure below is available.\n\n")
    }
    
    cat(sprintf("%-22s %14s %26s %12s\n", "Variable", "Permutation", sprintf("%.0f%% CI", alpha_pct), "Gini"))
    for (i in seq_len(nrow(imp_df))) {
      ci_txt <- if (is.finite(imp_df$perm_lower[i])) sprintf("%.4f to %.4f", imp_df$perm_lower[i], imp_df$perm_upper[i]) else "NA"
      cat(sprintf("%-22s %14s %26s %12s\n", imp_df$variable[i],
                  ifelse(is.finite(imp_df$permutation[i]), sprintf("%.4f", imp_df$permutation[i]), "NA"),
                  ci_txt,
                  ifelse(is.finite(imp_df$gini[i]), sprintf("%.2f", imp_df$gini[i]), "NA")))
    }
    cat("\n")
    
    top <- imp_df$variable[1]
    cat(sprintf("In plain English: '%s' matters most to this model's predictions.\n", top))
    non_sig <- imp_df$variable[which(imp_df$distinguishable_from_zero %in% FALSE)]
    if (length(non_sig) > 0) {
      cat(sprintf("%s: the permutation-importance CI includes zero (shuffling them barely hurts the model).\n",
                  paste(non_sig, collapse = ", ")))
      cat("Treat those as contributing little — the interval is approximate and descriptive, not a formal test.\n")
    }
    pos_perm <- pmax(imp_df$permutation, 0)
    if (!is.na(perm_col) && sum(pos_perm, na.rm = TRUE) > 0 && max(pos_perm, na.rm = TRUE) / sum(pos_perm, na.rm = TRUE) > 0.5) {
      cat(bold("NOTE:"), sprintf("'%s' alone accounts for over half of total permutation importance — worth\n", top))
      cat("checking this is clinically sensible and not a leaked or proxy-for-outcome variable.\n")
    }
    cat("\n")
    cat(bold("How NOT to read this:"), "importance says WHICH variables the model uses, not\n")
    cat("whether they raise or lower risk, nor by how much — and it is NOT a causal effect.\n")
    cat("Correlated predictors share (dilute) each other's importance. Gini importance is\n")
    cat("biased toward continuous variables and those with many categories; prefer permutation\n")
    cat("importance when the two disagree. For direction, use the partial dependence plots (pdp_vars).\n\n")
    out$importance <- imp_df
  } else {
    cat("No variable importance available for this forest.\n\n")
  }
  
  # ==========================================================
  # Q3. HOW WELL DOES IT CLASSIFY PATIENTS?
  # ==========================================================
  cat("========================================================\n")
  cat("Q3. HOW WELL DOES IT CLASSIFY PATIENTS?\n")
  cat("========================================================\n\n")
  
  conf_tab <- NULL; metrics_default <- NULL; metrics_youden <- NULL
  if (is_class) {
    pred_class <- model$predicted
    keep <- if (is_binary) ok else !is.na(pred_class)
    conf_tab <- table(Actual = y_obs[keep], Predicted = pred_class[keep])
    cat("Confusion matrix (rows = what actually happened, columns = what the model predicted):\n\n")
    print(conf_tab)
    cat("\n")
    out$confusion <- conf_tab
    
    if (is_binary) {
      pp <- prob_pos[ok]; ip <- is_pos[ok]
      m1 <- .rf_class_metrics(ip, (model$predicted == pos)[ok])
      metrics_default <- m1
      
      cat("At the DEFAULT rule (flag as event if at least 50% of trees vote event):\n")
      cat(sprintf("  Of the %s patients who actually HAD the event, the model flagged %s:  sensitivity = %s\n",
                  format(m1$tp + m1$fn, big.mark = ","), format(m1$tp, big.mark = ","), bold(.rf_pct(m1$sensitivity))))
      cat(sprintf("  Of the %s patients who did NOT, the model correctly cleared %s:        specificity = %s\n",
                  format(m1$tn + m1$fp, big.mark = ","), format(m1$tn, big.mark = ","), bold(.rf_pct(m1$specificity))))
      cat(sprintf("  Of the %s patients FLAGGED, %s truly had the event:                   PPV (precision) = %s\n",
                  format(m1$tp + m1$fp, big.mark = ","), format(m1$tp, big.mark = ","), .rf_pct(m1$ppv)))
      cat(sprintf("  Of the %s patients CLEARED, %s truly did not:                         NPV = %s\n",
                  format(m1$tn + m1$fn, big.mark = ","), format(m1$tn, big.mark = ","), .rf_pct(m1$npv)))
      cat(sprintf("  Balanced accuracy = %s, Cohen's kappa = %.3f\n\n", .rf_pct(m1$balanced_accuracy), m1$kappa))
      
      if (is.finite(m1$sensitivity) && is.finite(m1$specificity) && m1$sensitivity < 0.5 && m1$specificity > 0.8) {
        cat(bold("NOTE:"), "the default 50% rule is missing most true events (low sensitivity) while rarely\n")
        cat("raising false alarms — common when the event is uncommon. A lower threshold trades\n")
        cat("false alarms for catching more events; see the Youden threshold below.\n\n")
      }
      
      if (is.finite(youden_thr)) {
        m2 <- .rf_class_metrics(ip, pp >= youden_thr)
        metrics_youden <- m2
        cat(sprintf("At the YOUDEN-optimal threshold (predicted risk >= %.3f, the cut-off that maximises\n", youden_thr))
        cat("sensitivity + specificity - 1):\n")
        cat(sprintf("  sensitivity = %s, specificity = %s, PPV = %s, NPV = %s, balanced accuracy = %s\n\n",
                    .rf_pct(m2$sensitivity), .rf_pct(m2$specificity), .rf_pct(m2$ppv), .rf_pct(m2$npv),
                    .rf_pct(m2$balanced_accuracy)))
        cat("The best threshold depends on the cost of a missed event vs. a false alarm — Youden\n")
        cat("weights them equally, which is rarely exactly right clinically. Also remember PPV/NPV\n")
        cat("depend on how common the event is: the same model has a lower PPV in a lower-risk population.\n\n")
        out$threshold_youden <- youden_thr
      }
      out$metrics_default <- m1; out$metrics_youden <- metrics_youden
    } else {
      cat("Per-class OOB error rate (class.error): \n")
      print(round(model$confusion[, ncol(model$confusion), drop = FALSE], 3))
      cat("\n(Sensitivity/specificity-style metrics are reported for binary outcomes only.)\n\n")
    }
  } else {
    cat("(Classification metrics apply to classification forests only.)\n\n")
  }
  
  # ==========================================================
  # Q4. HOW DO PREDICTED RISKS DIFFER ACROSS PATIENTS?
  # ==========================================================
  cat("========================================================\n")
  cat("Q4. HOW DO PREDICTED RISKS DIFFER ACROSS PATIENTS?\n")
  cat("========================================================\n\n")
  
  tab_dec <- NULL
  if (is_class && is_binary) {
    pp <- prob_pos[ok]; yy <- y01[ok]; ip <- is_pos[ok]
    qs <- stats::quantile(pp, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    cat("Distribution of each patient's OOB predicted risk of the event:\n")
    cat(sprintf("  minimum %s | 25th pct %s | median %s | 75th pct %s | maximum %s | mean %s\n",
                .rf_pct(qs[[1]]), .rf_pct(qs[[2]]), .rf_pct(qs[[3]]), .rf_pct(qs[[4]]), .rf_pct(qs[[5]]), .rf_pct(mean(pp))))
    cat(sprintf("  Mean predicted risk (%s) vs. observed event rate (%s) overall.\n",
                .rf_pct(mean(pp)), .rf_pct(mean(yy))))
    if (length(unique(pp)) < 50) {
      cat("  (Few distinct predicted values — forest probabilities are coarse, so many patients tie.\n")
      cat("  Deciles below split ties by row order, which adds some arbitrary noise.)\n")
    }
    cat("\n")
    
    # ---- Decile table (always computed; printed only if requested) ----
    dec <- .rf_ntile(pp, 10)
    dec_levels <- sort(unique(dec))
    tab_dec <- do.call(rbind, lapply(dec_levels, function(d) {
      idx <- dec == d
      data.frame(risk_decile = d, n = sum(idx), events = sum(yy[idx]),
                 observed_event_rate = mean(yy[idx]), mean_predicted_risk = mean(pp[idx]))
    }))
    ci_dec <- .rf_wilson(tab_dec$events, tab_dec$n, conf_level)
    tab_dec$obs_lower <- ci_dec$lower; tab_dec$obs_upper <- ci_dec$upper
    tab_dec <- tab_dec[order(-tab_dec$risk_decile), , drop = FALSE]
    overall_rate <- mean(yy)
    tab_dec$lift <- tab_dec$observed_event_rate / overall_rate
    tab_dec$cum_pct_patients <- cumsum(tab_dec$n) / sum(tab_dec$n)
    tab_dec$cum_pct_events_captured <- cumsum(tab_dec$events) / sum(tab_dec$events)
    rownames(tab_dec) <- NULL
    out$deciles <- tab_dec
    
    top_row <- tab_dec[1, ]
    cat(sprintf("Top-risk decile (highest 10%% by predicted risk): observed event rate %s, versus %s overall\n",
                bold(.rf_pct(top_row$observed_event_rate)), .rf_pct(overall_rate)))
    if (nrow(tab_dec) >= 3) {
      cat(sprintf("(%.1fx the average), and it contains %s of ALL events. Screening just the top 3 deciles\n",
                  top_row$lift, .rf_pct(top_row$cum_pct_events_captured)))
      cat(sprintf("would capture %s of all events.\n\n", .rf_pct(tab_dec$cum_pct_events_captured[3])))
    } else {
      cat(sprintf("(%.1fx the average), and it contains %s of ALL events.\n\n",
                  top_row$lift, .rf_pct(top_row$cum_pct_events_captured)))
    }
    
    if (isTRUE(decile_table)) {
      cat("Risk decile table (decile 10 = highest predicted risk):\n")
      show_tab <- tab_dec
      for (cn in c("observed_event_rate", "mean_predicted_risk", "obs_lower", "obs_upper",
                   "cum_pct_patients", "cum_pct_events_captured")) show_tab[[cn]] <- round(show_tab[[cn]], 3)
      show_tab$lift <- round(show_tab$lift, 2)
      print(show_tab, row.names = FALSE)
      cat("\n")
    }
  } else if (!is_class) {
    cat("(Risk deciles and calibration apply to binary classification forests.)\n\n")
  } else {
    cat("(Risk deciles and calibration are only reported for binary outcomes.)\n\n")
  }
  
  # ---- Subgroup risk table (risk_by) ----
  risk_by_results <- list()
  if (!is.null(risk_by)) {
    if (is.null(data)) {
      cat(bold("NOTE:"), "risk_by needs the training data — pass data = <your data frame>.\n\n")
    } else if (!(is_class && is_binary)) {
      cat(bold("NOTE:"), "risk_by is only available for binary classification forests.\n\n")
    } else {
      for (v in risk_by) {
        if (!(v %in% names(data))) {
          warning(sprintf("risk_by variable '%s' was not found in `data` — skipped.", v), call. = FALSE)
          next
        }
        g <- .rf_group(data[[v]])[ok]
        pp <- prob_pos[ok]; yy <- y01[ok]
        lv <- levels(droplevels(g))
        rb <- do.call(rbind, lapply(lv, function(l) {
          idx <- which(g == l)
          data.frame(variable = v, group = l, n = length(idx), events = sum(yy[idx]),
                     observed_rate = mean(yy[idx]), mean_predicted = mean(pp[idx]))
        }))
        ci_rb <- .rf_wilson(rb$events, rb$n, conf_level)
        rb$obs_lower <- ci_rb$lower; rb$obs_upper <- ci_rb$upper
        risk_by_results[[v]] <- rb
        
        cat(sprintf("Predicted vs. observed risk by '%s':\n", v))
        show_rb <- rb
        for (cn in c("observed_rate", "mean_predicted", "obs_lower", "obs_upper")) show_rb[[cn]] <- round(show_rb[[cn]], 3)
        print(show_rb[, setdiff(names(show_rb), "variable")], row.names = FALSE)
        cat("\n")
      }
      if (length(risk_by_results) > 0) {
        cat("In plain English: where 'mean_predicted' sits close to 'observed_rate' the model is\n")
        cat("well calibrated for that subgroup; a large gap means it systematically over- or\n")
        cat("under-predicts risk for those patients (a fairness/transportability concern).\n\n")
        out$risk_by <- risk_by_results
      }
    }
  }
  
  # ==========================================================
  # FIGURES
  # ==========================================================
  
  # ---- Importance plot ----
  if (isTRUE(importance_plot)) {
    if (is.null(imp_df)) {
      warning("importance_plot = TRUE, but this forest has no importance values.", call. = FALSE)
    } else if (.rf_need_ggplot("importance_plot")) {
      perm_label <- sprintf("Permutation importance (%s, with %.0f%% CI)", ifelse(is_class, "MeanDecreaseAccuracy", "%IncMSE"), alpha_pct)
      gini_label <- sprintf("Node-purity importance (%s)", ifelse(is_class, "MeanDecreaseGini", "IncNodePurity"))
      ord <- imp_df$variable
      d_list <- list()
      if (any(is.finite(imp_df$permutation))) {
        d_list[[1]] <- data.frame(variable = imp_df$variable, measure = perm_label, kind = "point",
                                  value = imp_df$permutation,
                                  lo = ifelse(is.finite(imp_df$perm_lower), imp_df$perm_lower, imp_df$permutation),
                                  hi = ifelse(is.finite(imp_df$perm_upper), imp_df$perm_upper, imp_df$permutation))
      }
      if (any(is.finite(imp_df$gini))) {
        d_list[[length(d_list) + 1]] <- data.frame(variable = imp_df$variable, measure = gini_label, kind = "bar",
                                                   value = imp_df$gini, lo = NA_real_, hi = NA_real_)
      }
      d <- do.call(rbind, d_list)
      d$variable <- factor(d$variable, levels = rev(ord))
      d$measure <- factor(d$measure, levels = unique(d$measure))
      p <- ggplot2::ggplot(d, ggplot2::aes(x = variable, y = value)) +
        ggplot2::geom_col(data = d[d$kind == "bar", , drop = FALSE], fill = "steelblue") +
        ggplot2::geom_pointrange(data = d[d$kind == "point", , drop = FALSE],
                                 ggplot2::aes(ymin = lo, ymax = hi), colour = "darkblue") +
        ggplot2::facet_wrap(~ measure, scales = "free_y") +
        ggplot2::coord_flip() +
        ggplot2::labs(title = sprintf("Random Forest variable importance for '%s'", outcome),
                      subtitle = "Higher = the model relies on the variable more. Says nothing about direction of effect. CI is approximate.",
                      x = NULL, y = "Importance") +
        ggplot2::theme_minimal(base_size = 13)
      if (any(d$kind == "point")) {
        p <- p + ggplot2::geom_hline(data = data.frame(measure = factor(perm_label, levels = levels(d$measure)), y0 = 0),
                                     ggplot2::aes(yintercept = y0), linetype = "dashed", colour = "red", inherit.aes = FALSE)
      }
      .rf_print_plot(p)
    }
  }
  
  # ---- Importance table (full) ----
  if (isTRUE(importance_table) && !is.null(imp_df)) {
    cat("\nFull variable-importance table:\n")
    show_imp <- imp_df
    show_imp[, c("permutation", "perm_se", "perm_lower", "perm_upper")] <- round(show_imp[, c("permutation", "perm_se", "perm_lower", "perm_upper")], 5)
    show_imp$gini <- round(show_imp$gini, 3)
    print(show_imp, row.names = FALSE)
    cat("\n")
  }
  
  # ---- Error curve ----
  if (isTRUE(error_plot)) {
    if (.rf_need_ggplot("error_plot")) {
      if (is_class && !is.null(model$err.rate)) {
        er <- as.data.frame(model$err.rate)
        er$trees <- seq_len(nrow(er))
        long <- do.call(rbind, lapply(setdiff(names(er), "trees"), function(nm) {
          data.frame(trees = er$trees, error = er[[nm]], series = ifelse(nm == "OOB", "Overall OOB error", paste("Error for class", nm)))
        }))
      } else {
        long <- data.frame(trees = seq_along(model$mse), error = model$mse, series = "OOB mean squared error")
      }
      p <- ggplot2::ggplot(long, ggplot2::aes(x = trees, y = error, colour = series)) +
        ggplot2::geom_line(linewidth = 0.8) +
        ggplot2::labs(title = "Random Forest: out-of-bag error as trees are added",
                      subtitle = "Curves should flatten — if still falling at the right-hand edge, grow more trees",
                      x = "Number of trees", y = ifelse(is_class, "Error rate", "Mean squared error"), colour = NULL) +
        ggplot2::theme_minimal(base_size = 13)
      .rf_print_plot(p)
    }
  }
  
  # ---- Confusion matrix heatmap ----
  if (isTRUE(confusion_plot)) {
    if (is.null(conf_tab)) {
      warning("confusion_plot = TRUE requires a classification forest.", call. = FALSE)
    } else if (.rf_need_ggplot("confusion_plot")) {
      cdf <- as.data.frame(conf_tab, stringsAsFactors = FALSE)
      names(cdf) <- c("Actual", "Predicted", "n")
      row_tot <- tapply(cdf$n, cdf$Actual, sum)
      cdf$pct <- cdf$n / row_tot[cdf$Actual]
      cdf$label <- sprintf("%d\n(%.1f%% of actual %s)", cdf$n, 100 * cdf$pct, cdf$Actual)
      p <- ggplot2::ggplot(cdf, ggplot2::aes(x = Predicted, y = Actual, fill = pct)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::geom_text(ggplot2::aes(label = label), size = 4) +
        ggplot2::scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1), name = "% of actual") +
        ggplot2::labs(title = "Confusion matrix: OOB predictions vs. actual outcome",
                      subtitle = "Each row sums to 100% — the diagonal shows how often each true outcome is classified correctly",
                      x = "Predicted outcome", y = "Actual outcome") +
        ggplot2::theme_minimal(base_size = 13)
      .rf_print_plot(p)
    }
  }
  
  # ---- ROC curve ----
  if (isTRUE(roc_plot)) {
    if (is.null(roc_obj)) {
      warning("roc_plot = TRUE requires a binary classification forest with both outcome groups present.", call. = FALSE)
    } else if (.rf_need_ggplot("roc_plot")) {
      rdf <- data.frame(fpr = roc_obj$fpr, tpr = roc_obj$tpr)
      pts <- data.frame(fpr = 1 - metrics_default$specificity, tpr = metrics_default$sensitivity, label = "Default rule (50% of trees)")
      if (!is.null(metrics_youden)) {
        pts <- rbind(pts, data.frame(fpr = 1 - metrics_youden$specificity, tpr = metrics_youden$sensitivity,
                                     label = sprintf("Youden threshold (>= %.3f)", youden_thr)))
      }
      p <- ggplot2::ggplot(rdf, ggplot2::aes(x = fpr, y = tpr)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::geom_path(linewidth = 1, colour = "steelblue") +
        ggplot2::geom_point(data = pts, ggplot2::aes(colour = label), size = 3) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(title = "ROC curve: Random Forest (out-of-bag predictions)",
                      subtitle = sprintf("AUC = %.3f (%.0f%% CI %.3f to %.3f). Dashed line = chance.",
                                         roc_obj$auc, alpha_pct, out$auc_ci[["lower"]], out$auc_ci[["upper"]]),
                      x = "1 - Specificity (false-positive rate)", y = "Sensitivity (true-positive rate)", colour = NULL) +
        ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(legend.position = "bottom")
      .rf_print_plot(p)
    }
  }
  
  # ---- Predicted-risk distribution ----
  if (isTRUE(risk_distribution_plot)) {
    if (!(is_class && is_binary)) {
      warning("risk_distribution_plot = TRUE requires a binary classification forest.", call. = FALSE)
    } else if (.rf_need_ggplot("risk_distribution_plot")) {
      ddf <- data.frame(risk = prob_pos[ok], actual = factor(ifelse(is_pos[ok], paste0("Actually ", pos), paste0("Actually ", neg)),
                                                             levels = c(paste0("Actually ", neg), paste0("Actually ", pos))))
      p <- ggplot2::ggplot(ddf, ggplot2::aes(x = risk, fill = actual)) +
        ggplot2::geom_histogram(bins = 30, colour = "white") +
        ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", colour = "red") +
        ggplot2::facet_wrap(~ actual, ncol = 1, scales = "free_y") +
        ggplot2::scale_x_continuous(limits = c(-0.02, 1.02)) +
        ggplot2::labs(title = sprintf("Predicted risk of '%s' by what actually happened", pos),
                      subtitle = "Good discrimination: the event group sits to the right. Dashed line = default 50% rule.",
                      x = "OOB predicted risk (share of trees voting event)", y = "Number of patients", fill = NULL) +
        ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(legend.position = "none")
      .rf_print_plot(p)
    }
  }
  
  # ---- Calibration plot ----
  if (isTRUE(calibration_plot)) {
    if (is.null(tab_dec)) {
      warning("calibration_plot = TRUE requires a binary classification forest.", call. = FALSE)
    } else if (.rf_need_ggplot("calibration_plot")) {
      m <- max(c(tab_dec$obs_upper, tab_dec$mean_predicted_risk), na.rm = TRUE)
      p <- ggplot2::ggplot(tab_dec, ggplot2::aes(x = mean_predicted_risk, y = observed_event_rate)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = obs_lower, ymax = obs_upper), width = 0, colour = "steelblue") +
        ggplot2::geom_line(colour = "steelblue") +
        ggplot2::geom_point(size = 3, colour = "darkblue") +
        ggplot2::coord_equal(xlim = c(0, m), ylim = c(0, m)) +
        ggplot2::labs(title = "Calibration plot: observed vs. predicted risk by decile",
                      subtitle = sprintf("Points on the dashed line = well calibrated. Error bars: %.0f%% Wilson CI of the observed rate.", alpha_pct),
                      x = "Mean predicted risk (OOB)", y = "Observed event rate") +
        ggplot2::theme_minimal(base_size = 13)
      .rf_print_plot(p)
    }
  }
  
  # ---- Gains plot ----
  if (isTRUE(gains_plot)) {
    if (is.null(tab_dec)) {
      warning("gains_plot = TRUE requires a binary classification forest.", call. = FALSE)
    } else if (.rf_need_ggplot("gains_plot")) {
      gdf <- data.frame(pct_patients = c(0, tab_dec$cum_pct_patients), pct_events = c(0, tab_dec$cum_pct_events_captured))
      p <- ggplot2::ggplot(gdf, ggplot2::aes(x = pct_patients, y = pct_events)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::geom_line(linewidth = 1, colour = "steelblue") +
        ggplot2::geom_point(size = 2.5, colour = "darkblue") +
        ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
        ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
        ggplot2::labs(title = "Cumulative gains: how many events are captured by screening the highest-risk patients first",
                      subtitle = "Further above the dashed line (random selection) = the model concentrates events in the top deciles",
                      x = "% of patients screened (highest predicted risk first)", y = "% of all events captured") +
        ggplot2::theme_minimal(base_size = 13)
      .rf_print_plot(p)
    }
  }
  
  # ---- Risk by subgroup plots ----
  if (length(risk_by_results) > 0 && .rf_need_ggplot("risk_by plots")) {
    for (v in names(risk_by_results)) {
      rb <- risk_by_results[[v]]
      long <- rbind(
        data.frame(group = rb$group, series = "Observed event rate", value = rb$observed_rate, lo = rb$obs_lower, hi = rb$obs_upper),
        data.frame(group = rb$group, series = "Mean predicted risk", value = rb$mean_predicted, lo = NA_real_, hi = NA_real_)
      )
      long$group <- factor(long$group, levels = rb$group)
      p <- ggplot2::ggplot(long, ggplot2::aes(x = group, y = value, colour = series, shape = series)) +
        ggplot2::geom_pointrange(ggplot2::aes(ymin = lo, ymax = hi), position = ggplot2::position_dodge(width = 0.4), size = 0.7) +
        ggplot2::labs(title = sprintf("Predicted vs. observed risk by %s", v),
                      subtitle = "Model is well calibrated for a subgroup when the two markers sit together",
                      x = v, y = sprintf("Risk of '%s'", pos), colour = NULL, shape = NULL) +
        ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(legend.position = "bottom")
      .rf_print_plot(p)
    }
  }
  
  # ---- Partial dependence (+ ICE) ----
  pdp_results <- list()
  if (!is.null(pdp_vars)) {
    if (is.null(data)) {
      cat(bold("NOTE:"), "pdp_vars needs the training data — pass data = <your data frame>.\n\n")
    } else {
      if (identical(pdp_vars, "auto")) {
        cand <- if (!is.null(imp_df)) imp_df$variable else predictors
        pdp_vars <- utils::head(cand[cand %in% names(data)], 3)
      }
      if (length(pdp_vars) > 0 && .rf_need_ggplot("pdp_vars")) {
        cat("========================================================\n")
        cat("PARTIAL DEPENDENCE — how predicted risk changes as one variable changes\n")
        cat("========================================================\n")
        cat("For each value of the chosen variable, EVERY patient is given that value (keeping\n")
        cat("their own values for everything else), the model predicts, and the predictions\n")
        cat(sprintf("are averaged. The curve shows the DIRECTION and SHAPE of the model's use of that\n"))
        cat("variable — which importance alone cannot.\n")
        cat("Caution: if predictors are strongly correlated, this averages over unrealistic\n")
        cat("combinations (e.g. very old patients with very young-patient values elsewhere).\n")
        if (isTRUE(pdp_ice)) cat("Grey lines are individual patients (ICE) — fanning lines mean the effect differs by patient.\n")
        cat("\n")
        focus_label <- if (is_class) sprintf("Average predicted probability of '%s'", focus) else sprintf("Average predicted %s", outcome)
        
        for (v in pdp_vars) {
          if (!(v %in% names(data))) {
            warning(sprintf("pdp variable '%s' was not found in `data` — skipped.", v), call. = FALSE)
            next
          }
          pd <- tryCatch(.rf_pdp(model, data, v, is_class, focus, ice = isTRUE(pdp_ice)), error = function(e) {
            warning(sprintf("Could not compute partial dependence for '%s': %s", v, conditionMessage(e)), call. = FALSE)
            NULL
          })
          if (is.null(pd)) next
          pdp_results[[v]] <- pd$pdp
          
          rng <- range(pd$pdp$yhat, na.rm = TRUE)
          cat(sprintf("> %s: average predicted %s ranges from %s to %s across its values",
                      v, ifelse(is_class, "risk", outcome),
                      ifelse(is_class, .rf_pct(rng[1]), sprintf("%.3f", rng[1])),
                      ifelse(is_class, .rf_pct(rng[2]), sprintf("%.3f", rng[2]))))
          if (!pd$is_cat) {
            ordr <- pd$pdp[order(pd$pdp$x_value), ]
            trend <- stats::cor(ordr$x_value, ordr$yhat, method = "spearman")
            cat(sprintf("; overall %s.\n",
                        if (!is.finite(trend)) "flat"
                        else if (trend > 0.7) "rising (higher values -> higher predicted risk)"
                        else if (trend < -0.7) "falling (higher values -> lower predicted risk)"
                        else "non-monotonic or fairly flat"))
          } else cat(".\n")
          
          p <- ggplot2::ggplot(pd$pdp, ggplot2::aes(x = x_value, y = yhat))
          if (!pd$is_cat) {
            if (isTRUE(pdp_ice) && !is.null(pd$ice)) {
              p <- p + ggplot2::geom_line(data = pd$ice, ggplot2::aes(x = x_value, y = yhat, group = row_id),
                                          colour = "grey70", alpha = 0.5, inherit.aes = FALSE)
            }
            p <- p + ggplot2::geom_line(linewidth = 1.3, colour = "darkblue") +
              ggplot2::geom_rug(data = data.frame(x_value = data[[v]]), ggplot2::aes(x = x_value), inherit.aes = FALSE, alpha = 0.2, sides = "b")
          } else {
            if (isTRUE(pdp_ice) && !is.null(pd$ice)) {
              pd$ice$x_value <- factor(pd$ice$x_value, levels = levels(pd$pdp$x_value))
              p <- p + ggplot2::geom_line(data = pd$ice, ggplot2::aes(x = x_value, y = yhat, group = row_id),
                                          colour = "grey70", alpha = 0.4, inherit.aes = FALSE)
            }
            p <- p + ggplot2::geom_col(fill = "steelblue", width = 0.6, alpha = 0.8)
          }
          p <- p + ggplot2::labs(title = sprintf("Partial dependence: %s", v),
                                 subtitle = sprintf("Based on %d patients%s. Shows how the model's prediction changes with %s, averaging over everyone else.",
                                                    pd$n_rows_used, ifelse(pd$n_rows_used < nrow(data), " (evenly spaced subsample)", ""), v),
                                 x = v, y = focus_label) +
            ggplot2::theme_minimal(base_size = 13)
          .rf_print_plot(p)
        }
        cat("\n")
        out$pdp <- pdp_results
      }
    }
  }
  
  # ==========================================================
  # CAVEATS / APPROPRIATENESS
  # ==========================================================
  cat("########################################################\n")
  cat("# IS A RANDOM FOREST THE RIGHT TOOL — AND WHAT IT CAN'T TELL YOU\n")
  cat("########################################################\n\n")
  cat(bold("Good for:"), "prediction with minimal tuning; capturing non-linearities and interactions\n")
  cat("automatically; a built-in honest error estimate (OOB); robustness to outliers.\n\n")
  cat(bold("Not for:"), "\n")
  cat(" - CAUSAL questions. This answers 'can we PREDICT the event?', not 'does X cause it?'.\n")
  cat("   Importance and partial dependence describe the model, not the biology, and are\n")
  cat("   confounded by anything left out. Use an adjusted regression (explain_glm()) for\n")
  cat("   effect estimates with confidence intervals.\n")
  cat(" - Formal inference. There are no p-values or coefficient CIs for forest predictors;\n")
  cat("   the CIs above are for performance metrics and an approximate importance interval.\n")
  cat(" - Final validation. OOB is internal validation: it is optimistic about performance in\n")
  cat("   a new setting. Before relying on this model, validate it on external or temporal data\n")
  cat("   (see the TRIPOD guidance on reporting prediction models) and re-check calibration there.\n")
  cat(" - Calibrated absolute risks. Forest probabilities are often compressed and coarse;\n")
  cat("   if you need trustworthy risk estimates (not just rankings), recalibrate (e.g. logistic\n")
  cat("   or isotonic recalibration) and report the calibration plot.\n\n")
  cat("A fair comparison: fit a logistic regression with the same predictors and compare\n")
  cat("cross-validated AUC and calibration. With few predictors and modest numbers of events,\n")
  cat("regression is often as accurate, far more interpretable, and better calibrated.\n")
  cat("========================================================\n")
  
  invisible(out)
}


# ============================================================
# USAGE EXAMPLES (using the cohort dataset from earlier)
# ============================================================
library(randomForest)

out_rf <- randomForest(
  event_f ~ age + sex + bmi + comorbidity_score + sbp,
  data = cohort, importance = TRUE
)

# Text answers to the four questions only:
explain_rf(out_rf)

# Add whichever tables/figures you want:
explain_rf(out_rf, data = cohort,
           importance_plot = TRUE, importance_table = TRUE,
           confusion_plot = TRUE, roc_plot = TRUE,
           risk_distribution_plot = TRUE, decile_table = TRUE,
           calibration_plot = TRUE, gains_plot = TRUE, error_plot = TRUE,
           pdp_vars = c("age", "bmi"), pdp_ice = TRUE,
           risk_by = c("sex"))

# ...or everything at once:
res <- explain_rf(out_rf, data = cohort, all_outputs = TRUE)
res$importance; res$deciles; res$metrics_youden; res$auc

# A different level as the "event":
explain_rf(out_rf, positive_class = "No")

# ============================================================
# EDGE CASES THIS FUNCTION HAS BEEN SPECIFICALLY CHECKED AGAINST
# ============================================================
#  1. Not a randomForest object / unsupervised forest        -> stops with a clear message
#  2. model$y missing                                        -> stops with a specific message
#  3. positive_class not among the outcome levels            -> stops, lists the valid levels
#  4. Default positive class                                 -> "Yes" if present, else the 2nd
#     level, and the choice is STATED in the output so it can't be silently inverted
#  5. Patients with no OOB prediction (votes row sums to 0)  -> set to NA, excluded, flagged
#  6. Small minority group (<50) / strong class imbalance    -> flagged before results
#  7. Error curve still falling at the last tree             -> flagged, suggests more trees
#  8. Fit without importance = TRUE                          -> flagged; importance sections degrade
#     gracefully (Gini only / skipped) instead of erroring
#  9. importanceSD stored as a matrix (classification) or a named vector (regression) -> both handled
# 10. pROC not installed or erroring                         -> Hanley-McNeil CI fallback, method named
# 11. Tied predicted probabilities (very common with forests)-> ROC/AUC handle ties as diagonal
#     segments; deciles break ties by row order and say so when predictions are coarse
# 12. Calibration slope with exact 0/1 forest votes          -> probabilities clipped to [0.01, 0.99]
#     before the logit, and the output says the slope/intercept are approximate
# 13. Multi-class outcome                                    -> confusion matrix + per-class error +
#     importance + error curve only; binary-only analyses are skipped with a message
# 14. Regression forest                                      -> OOB RMSE and % variance explained,
#     importance, error curve, partial dependence; classification-only outputs warn instead of error
# 15. data not supplied for pdp_vars/risk_by                 -> tries the model's own call; if that
#     fails or the row count differs, a specific message says to pass data = ...
# 16. A requested pdp/risk_by variable not in the data       -> warned and skipped, others continue
# 17. Numeric vs factor predictors in partial dependence     -> numeric: line over the 5th-95th
#     percentile grid (avoids extrapolating into the tails); factor: bars over its levels with
#     factor levels reconstructed exactly (predict() requires identical levels)
# 18. Large datasets in partial dependence                   -> evenly spaced subsample of up to 500
#     rows, deterministic (no RNG, so your random seed is never disturbed)
# 19. ggplot2 not installed                                  -> every figure warns with install
#     instructions; the text output and tables still work
# 20. Each figure's prerequisites missing (e.g. roc_plot on a regression forest) -> warns rather
#     than erroring
# 21. Numeric grouping for risk_by with few unique values or duplicate quantile breaks -> falls
#     back to treating the variable as categorical instead of crashing in cut()
# 22. A decile bucket with zero events or all events         -> Wilson intervals stay inside [0, 1]
#     (a plain Wald interval would not)
# ============================================================