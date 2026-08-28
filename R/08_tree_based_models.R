# 08_tree_based_models

# ============================================================
# 8. TREE-BASED METHODS — outcome: event_f
# ============================================================
library(rpart)
library(rpart.plot)

model_tree <- rpart(event_f ~ age + bmi + comorbidity_score + sex + treatment,
                    data = cohort, method = "class")
out_tree_printcp <- printcp(model_tree)
out_tree_printcp
# rpart.plot(model_tree)  # visual, not an object

library(randomForest)
model_rf <- randomForest(event_f ~ age + bmi + comorbidity_score + sex + treatment,
                         data = cohort, ntree = 500, importance = TRUE)
out_rf_print      <- model_rf
out_rf_importance <- importance(model_rf)
out_rf_print
out_rf_importance
# varImpPlot(model_rf)  # visual, not an object

library(xgboost)
x_xgb <- model.matrix(event_f ~ age + bmi + comorbidity_score + sex + treatment,
                      data = cohort)[, -1]
y_xgb <- cohort$event
dtrain <- xgb.DMatrix(data = x_xgb, label = y_xgb)
model_xgb <- xgb.train(data = dtrain, nrounds = 100,
                       params = list(objective = "binary:logistic", max_depth = 4, eta = 0.1),
                       verbose = 0)
out_xgb_importance <- xgb.importance(model = model_xgb)
out_xgb_importance








# ============================================================
# explain_trees.R — plain-language interpretation of tree-based
# models: rpart() (single decision tree), randomForest(), and
# xgboost() (gradient boosting).
#
# Like explain_glmnet.R, these methods do NOT produce p-values or
# confidence intervals — they're prediction-focused machine
# learning methods, not classical statistical inference. This
# file never fabricates a p-value/CI for any of them. Instead,
# each function reports what these methods DO give you honestly:
# prediction performance, variable importance (with important
# caveats about what "importance" does and doesn't mean), and
# model structure.
#
# Three functions:
#   explain_rpart(model, plot_tree = FALSE)
#   explain_rf(model, plot_importance = FALSE)
#   explain_xgb(model, plot_importance = FALSE)
# ============================================================

library(rpart)
library(rpart.plot)
library(randomForest)
library(xgboost)

# ---- Shared helper ----
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) crayon::bold(x) else x
}


# ============================================================
# 1. explain_rpart() — a single decision tree
# ============================================================
explain_rpart <- function(model, plot_tree = FALSE) {
  
  if (!inherits(model, "rpart")) {
    stop("explain_rpart() expects a model fit with rpart::rpart().")
  }
  
  method_used <- tryCatch({
    m <- model$method
    if (is.null(m)) NA else m
  }, error = function(e) NA)
  outcome <- tryCatch(deparse(model$terms[[2]]), error = function(e) "the outcome")
  
  cat("========================================================\n")
  cat("DECISION TREE (rpart)\n")
  cat("========================================================\n\n")
  
  cat("What a decision tree does, in plain English:\n")
  cat("It repeatedly asks the single most useful yes/no question it can find (e.g.\n")
  cat("'is age > 60?') to split patients into increasingly similar groups, then\n")
  cat("predicts the most common outcome (classification) or average outcome\n")
  cat("(regression) within each final group ('leaf'). The big advantage over other\n")
  cat("methods here is that you can literally read off the decision logic.\n\n")
  
  # ==========================================================
  # PRE-FLIGHT CHECKS
  # ==========================================================
  diagnostics <- character(0)
  
  n_splits <- tryCatch({
    sp <- model$splits
    if (is.null(sp)) 0L else nrow(sp)
  }, error = function(e) NA)
  
  if (is.na(n_splits) || n_splits == 0 || nrow(model$frame) == 1) {
    diagnostics <- c(diagnostics,
                     "made ZERO splits — the tree is just a single root node predicting the same value for every patient. This usually means none of your predictors were useful enough (given rpart's default complexity settings) to justify a split. Check your data, or try rpart(..., control = rpart.control(cp = 0)) to force at least some splits and inspect what's going on.")
  }
  
  if (is.na(method_used)) {
    diagnostics <- c(diagnostics, "has an undetermined fitting method (model$method was not readable) — some of the explanations below may not be quite right; double-check manually.")
  }
  
  if (length(diagnostics) > 0) {
    cat("########################################################\n")
    cat("# MODEL DIAGNOSTIC WARNINGS\n")
    cat("########################################################\n")
    for (d in diagnostics) cat(bold("WARNING:"), "this model", d, "\n\n")
    if (n_splits == 0 || is.na(n_splits)) {
      cat("========================================================\n")
      return(invisible(NULL))
    }
  }
  
  n_obs <- tryCatch(sum(model$frame$n[1]), error = function(e) NA)
  n_leaves <- tryCatch(sum(model$frame$var == "<leaf>"), error = function(e) NA)
  
  cat(sprintf("This tree was fit on %s patients, uses %s split(s), and ends in %s\n",
              bold(format(n_obs, big.mark = ",")), bold(n_splits), bold(n_leaves)))
  cat("final groups ('leaves') — each leaf is a distinct combination of yes/no\n")
  cat("answers that leads to its own prediction.\n\n")
  
  # ==========================================================
  # VARIABLES ACTUALLY USED
  # ==========================================================
  var_imp <- tryCatch(model$variable.importance, error = function(e) NULL)
  
  if (!is.null(var_imp) && length(var_imp) > 0) {
    cat("Variable importance (how much each variable contributed to reducing\n")
    cat("prediction error across all its splits and 'near-miss' alternative splits):\n")
    cat("--------------------------------------------------------\n")
    var_imp_sorted <- sort(var_imp, decreasing = TRUE)
    for (vn in names(var_imp_sorted)) {
      cat(sprintf("  %-20s %.2f\n", vn, var_imp_sorted[[vn]]))
    }
    cat("\nIn plain English: HIGHER numbers mean the variable was more useful for\n")
    cat("splitting patients into more accurate groups. This tells you WHICH variables\n")
    cat("mattered, but — unlike an odds ratio — NOT the direction or size of their\n")
    cat("effect. Read the tree diagram itself (plot_tree = TRUE) to see the actual\n")
    cat("direction of each split (e.g. whether older or younger patients go left).\n\n")
    
    unused <- setdiff(attr(model$terms, "term.labels"), names(var_imp))
    if (length(unused) > 0) {
      cat(sprintf("Variable(s) NEVER used in any split: %s\n", paste(unused, collapse = ", ")))
      cat("These were available to the tree but never chosen — not necessarily because\n")
      cat("they're irrelevant, but because other variables were more useful FIRST (once\n")
      cat("a variable like age splits the data well, a correlated variable often has\n")
      cat("little left to add).\n\n")
    }
  } else {
    cat("(No variable importance available — this can happen with a very simple tree.)\n\n")
  }
  
  # ==========================================================
  # CP TABLE — pruning guidance
  # ==========================================================
  cp_table <- tryCatch(model$cptable, error = function(e) NULL)
  
  if (!is.null(cp_table) && nrow(cp_table) > 0) {
    cat("Complexity / pruning table:\n")
    cat("--------------------------------------------------------\n")
    cat("Each row is a different 'how big should this tree be' option. CP (complexity\n")
    cat("parameter) controls how much each split must improve the fit to be worth\n")
    cat("keeping — smaller CP = bigger, more detailed (and more overfitting-prone) tree.\n")
    cat("xerror is the CROSS-VALIDATED error (a more honest estimate of how the tree\n")
    cat("will perform on NEW patients than the training error, which will look\n")
    cat("artificially good since the tree was built on the same data).\n\n")
    
    print(round(cp_table, 4))
    
    best_idx <- which.min(cp_table[, "xerror"])
    best_cp  <- cp_table[best_idx, "CP"]
    
    cat(sprintf("\nThe row with the LOWEST cross-validated error (xerror) uses CP = %.4f",
                best_cp))
    cat(sprintf(" (%d split(s)).\n", cp_table[best_idx, "nsplit"]))
    
    if (best_idx < nrow(cp_table)) {
      xerror_best <- cp_table[best_idx, "xerror"]
      xstd_best   <- cp_table[best_idx, "xstd"]
      within_1se  <- which(cp_table[, "xerror"] <= xerror_best + xstd_best)
      simplest_1se_idx <- min(within_1se)
      
      if (simplest_1se_idx < best_idx) {
        cat(sprintf("A SIMPLER tree (CP = %.4f, %d split(s)) is within 1 standard error of\n",
                    cp_table[simplest_1se_idx, "CP"], cp_table[simplest_1se_idx, "nsplit"]))
        cat("that best result — the common '1-SE rule' says this simpler tree is often\n")
        cat("preferable, since it's statistically about as good but less likely to be\n")
        cat("overfit to this particular dataset. To prune to it:\n")
        cat(sprintf("  model_pruned <- prune(model, cp = %.4f)\n", cp_table[simplest_1se_idx, "CP"]))
      }
    }
    cat("\n")
  }
  
  cat("--------------------------------------------------------\n")
  cat(bold("CAVEATS — is a single tree the right tool here?\n"))
  cat("A single decision tree is HIGHLY interpretable (you can trace the exact logic\n")
  cat("for any prediction) but also HIGH VARIANCE — a slightly different sample of\n")
  cat("patients can produce a meaningfully different tree. It's a great tool for\n")
  cat("exploring/communicating decision logic, but for the best possible PREDICTION\n")
  cat("accuracy, an ensemble method (Random Forest, XGBoost — see explain_rf()/\n")
  cat("explain_xgb()) will almost always outperform a single tree, at the cost of\n")
  cat("losing this direct interpretability.\n")
  cat("No p-values or confidence intervals are produced by this method — that's\n")
  cat("expected, not a gap in this function; trees are a prediction/exploration tool,\n")
  cat("not a classical inferential one.\n")
  cat("========================================================\n")
  
  if (isTRUE(plot_tree)) {
    tryCatch({
      rpart.plot::rpart.plot(model, shadow.col = "gray", extra = 101)
    }, error = function(e) {
      warning(sprintf("Could not draw the tree plot: %s", conditionMessage(e)))
    })
  }
  
  invisible(list(cp_table = cp_table, variable_importance = var_imp))
}


# ============================================================
# 2. explain_rf() — Random Forest
# ============================================================
explain_rf <- function(model, plot_importance = FALSE) {
  
  if (!inherits(model, "randomForest")) {
    stop("explain_rf() expects a model fit with randomForest::randomForest().")
  }
  
  model_type <- tryCatch({
    t <- model$type
    if (is.null(t)) NA else t
  }, error = function(e) NA)  # "classification"/"regression"/"unsupervised"
  is_classification <- !is.na(model_type) && model_type == "classification"
  
  cat("========================================================\n")
  cat("RANDOM FOREST\n")
  cat("========================================================\n\n")
  
  cat("What a random forest does, in plain English:\n")
  cat("It builds HUNDREDS of individual decision trees, each on a random subset of\n")
  cat("patients AND a random subset of variables at each split, then averages their\n")
  cat("predictions (votes, for classification). This 'wisdom of crowds' approach fixes\n")
  cat("a single tree's biggest weakness (high variance/instability) at the cost of\n")
  cat("losing the ability to read off simple decision logic.\n\n")
  
  # ==========================================================
  # PRE-FLIGHT CHECKS
  # ==========================================================
  diagnostics <- character(0)
  
  n_trees <- tryCatch(model$ntree, error = function(e) NA)
  err_rate <- tryCatch(model$err.rate, error = function(e) NULL)
  mse_path <- tryCatch(model$mse, error = function(e) NULL)
  
  # Check whether the OOB error curve has stabilized by the end — if it's still
  # dropping meaningfully in the last 10% of trees, more trees are likely to help
  error_curve <- if (is_classification) err_rate[, "OOB"] else mse_path
  if (!is.null(error_curve) && length(error_curve) >= 20) {
    last_10pct <- tail(error_curve, ceiling(length(error_curve) * 0.1))
    early_ref  <- error_curve[ceiling(length(error_curve) * 0.5)]
    if (!is.na(early_ref) && early_ref != 0) {
      relative_drop <- (early_ref - mean(last_10pct)) / early_ref
      if (relative_drop > 0.02) {
        diagnostics <- c(diagnostics, sprintf(
          "may benefit from MORE trees — the error was still improving by roughly %.1f%% between the midpoint and the end of the %d trees fitted. Consider refitting with a larger ntree (e.g. 1000+) and checking that the error curve has clearly flattened (plot(model) shows this visually).",
          relative_drop * 100, n_trees))
      }
    }
  }
  
  importance_available <- tryCatch(!is.null(model$importance) && ncol(model$importance) > 0,
                                   error = function(e) FALSE)
  if (!importance_available) {
    diagnostics <- c(diagnostics,
                     "was fit WITHOUT importance = TRUE, so no variable importance is available. Refit with randomForest(..., importance = TRUE) if you want it.")
  }
  
  if (is_classification) {
    confusion <- tryCatch(model$confusion, error = function(e) NULL)
    if (!is.null(confusion)) {
      class_counts <- rowSums(confusion[, colnames(confusion) != "class.error", drop = FALSE])
      if (max(class_counts) / min(class_counts) > 9) {
        diagnostics <- c(diagnostics, sprintf(
          "has strongly IMBALANCED outcome classes (%s). Raw OOB accuracy can look artificially good on imbalanced data just by favoring the majority class — check the PER-CLASS error rates below carefully, not just the overall error.",
          paste(sprintf("%s: %d", names(class_counts), class_counts), collapse = ", ")))
      }
    }
  }
  
  if (length(diagnostics) > 0) {
    cat("########################################################\n")
    cat("# MODEL DIAGNOSTIC WARNINGS\n")
    cat("########################################################\n")
    for (d in diagnostics) cat(bold("WARNING:"), "this model", d, "\n\n")
  }
  
  cat(sprintf("Forest size: %s trees | Model type: %s\n\n", bold(n_trees),
              bold(ifelse(is.na(model_type), "unknown", model_type))))
  
  # ==========================================================
  # OOB PERFORMANCE
  # ==========================================================
  cat("Out-of-Bag (OOB) performance estimate:\n")
  cat("--------------------------------------------------------\n")
  cat("Each tree in the forest is built using only ~63% of patients (randomly\n")
  cat("resampled with replacement); the other ~37% ('out-of-bag') are natural,\n")
  cat("automatic held-out test cases for THAT tree. Averaging each patient's\n")
  cat("prediction across only the trees that didn't see them gives an honest\n")
  cat("performance estimate WITHOUT needing to set aside a separate test set.\n\n")
  
  if (is_classification) {
    confusion <- tryCatch(model$confusion, error = function(e) NULL)
    if (!is.null(confusion)) {
      oob_error <- tail(err_rate[, "OOB"], 1)
      cat(sprintf("Overall OOB error rate: %.1f%% (i.e. %.1f%% of predictions were correct)\n\n",
                  oob_error * 100, 100 - oob_error * 100))
      
      cat("Confusion matrix (rows = TRUE outcome, columns = PREDICTED outcome):\n")
      print(confusion)
      cat("\nIn plain English: the 'class.error' column shows how often the model got\n")
      cat("THAT specific true outcome wrong. A model can have good overall accuracy while\n")
      cat("still being much worse at predicting one particular outcome than the other —\n")
      cat("always check both classes, not just the overall error rate.\n\n")
    }
  } else {
    pct_var_explained <- tryCatch(tail(model$rsq, 1) * 100, error = function(e) NA)
    if (!is.na(pct_var_explained)) {
      cat(sprintf("Percent variance explained (OOB): %.1f%%\n", pct_var_explained))
      cat("(This is a random-forest-flavoured version of R-squared — how much of the\n")
      cat("variation in the outcome the forest's predictions account for.)\n\n")
    }
  }
  
  # ==========================================================
  # VARIABLE IMPORTANCE
  # ==========================================================
  if (importance_available) {
    cat("Variable importance:\n")
    cat("--------------------------------------------------------\n")
    imp <- model$importance
    
    has_accuracy_col <- "MeanDecreaseAccuracy" %in% colnames(imp)
    has_gini_col <- "MeanDecreaseGini" %in% colnames(imp)
    
    if (has_accuracy_col) {
      cat("MeanDecreaseAccuracy — how much WORSE the model's predictions get if you\n")
      cat("randomly scramble this variable's values (breaking any real relationship it\n")
      cat("had with the outcome). Bigger = more important for PREDICTION ACCURACY.\n\n")
      sorted_acc <- imp[order(-imp[, "MeanDecreaseAccuracy"]), , drop = FALSE]
      for (vn in rownames(sorted_acc)) {
        cat(sprintf("  %-20s %.3f\n", vn, sorted_acc[vn, "MeanDecreaseAccuracy"]))
      }
      cat("\n")
    }
    
    if (has_gini_col) {
      cat("MeanDecreaseGini — how much this variable helps make the groups it splits\n")
      cat("PURER (more homogeneous) on average, summed across all its uses in the forest.\n")
      cat(sprintf("%s this measure is known to be BIASED toward continuous variables and\n", bold("CAVEAT:")))
      cat("those with many categories (they get more opportunities to find a good-looking\n")
      cat("split just by chance) — prefer MeanDecreaseAccuracy above when the two disagree.\n\n")
    }
    
    if (!has_accuracy_col && !has_gini_col) {
      print(imp)
      cat("\n")
    }
  }
  
  cat("--------------------------------------------------------\n")
  cat(bold("CAVEATS — is a random forest the right tool here?\n"))
  cat("Good for: strong PREDICTION accuracy with minimal tuning, automatically handles\n")
  cat("non-linear relationships and variable interactions without you specifying them,\n")
  cat("robust to outliers, gives you a built-in honest performance estimate (OOB) with\n")
  cat("no separate validation set needed.\n")
  cat("Weaker for: precise, directional effect-size interpretation (importance tells you\n")
  cat("WHICH variables matter, not HOW — not a substitute for a regression's odds\n")
  cat("ratios); correlated predictors DILUTE each other's apparent importance (splitting\n")
  cat("credit between them, making each look less important than it really is);\n")
  cat("extrapolation beyond the range of the training data is unreliable, since trees\n")
  cat("can only ever predict values they've seen.\n")
  cat("No p-values or confidence intervals are produced by this method.\n")
  cat("========================================================\n")
  
  if (isTRUE(plot_importance)) {
    if (!importance_available) {
      warning("plot_importance = TRUE was requested, but this model has no importance to plot (refit with importance = TRUE).")
    } else {
      tryCatch({
        randomForest::varImpPlot(model)
      }, error = function(e) {
        warning(sprintf("Could not draw the importance plot: %s", conditionMessage(e)))
      })
    }
  }
  
  invisible(list(
    oob_error = if (is_classification) tail(err_rate[, "OOB"], 1) else NA,
    importance = if (importance_available) model$importance else NULL
  ))
}


# ============================================================
# 3. explain_xgb() — XGBoost (gradient boosting)
# ============================================================
explain_xgb <- function(model, feature_names = NULL, plot_importance = FALSE) {
  
  if (!inherits(model, "xgb.Booster")) {
    stop("explain_xgb() expects a model fit with xgboost::xgboost() or xgb.train().")
  }
  
  cat("========================================================\n")
  cat("XGBOOST (Gradient Boosted Trees)\n")
  cat("========================================================\n\n")
  
  cat("What gradient boosting does, in plain English:\n")
  cat("Instead of building many INDEPENDENT trees and averaging them (like a random\n")
  cat("forest), boosting builds trees ONE AT A TIME, where each new tree focuses\n")
  cat("specifically on correcting the mistakes the trees before it made. This\n")
  cat("typically gives even better prediction accuracy than a random forest, but is\n")
  cat("MORE prone to overfitting if not carefully tuned/validated (each tree is\n")
  cat("chasing the previous tree's errors, including any that were just noise).\n\n")
  
  
  # ==========================================================
  # MODEL SETTINGS RECAP
  # ==========================================================
  params <- tryCatch(model$params, error = function(e) list())
  n_rounds <- tryCatch({
    nr <- model$niter
    if (is.null(nr)) NA else nr
  }, error = function(e) NA)
  
  cat("Model settings used:\n")
  cat("--------------------------------------------------------\n")
  cat(sprintf("  Number of boosting rounds (trees): %s\n", ifelse(is.na(n_rounds), "unknown", n_rounds)))
  cat(sprintf("  max_depth (how deep each individual tree can go): %s\n",
              ifelse(is.null(params$max_depth), "default", params$max_depth)))
  cat("    -> deeper trees can capture more complex patterns, but are more prone to\n")
  cat("       overfitting and memorizing noise in the training data.\n")
  cat(sprintf("  eta / learning rate (how much each tree corrects the last): %s\n",
              ifelse(is.null(params$eta), "default", params$eta)))
  cat("    -> smaller = slower, more cautious learning (usually needs MORE rounds to\n")
  cat("       compensate, but is less prone to overfitting than a large eta).\n")
  cat(sprintf("  Objective: %s\n\n", ifelse(is.null(params$objective), "unknown", params$objective)))
  
  
  # ==========================================================
  # THE SINGLE MOST IMPORTANT CAVEAT FOR THIS SPECIFIC WORKFLOW
  # ==========================================================
  cat("########################################################\n")
  cat("# CRITICAL CAVEAT — read this before trusting this model\n")
  cat("########################################################\n")
  cat(bold("This model was fit with xgboost() directly on the FULL training data, with a\n"))
  cat(bold("FIXED number of rounds/depth/learning rate chosen in advance — with no\n"))
  cat(bold("cross-validation or held-out test set to check whether it's overfitting.\n\n"))
  cat("XGBoost is powerful specifically because it can fit training data extremely\n")
  cat("well — which means WITHOUT a validation step, you have NO way to tell whether\n")
  cat("this model will generalize to new patients, or has simply memorized this\n")
  cat("particular dataset (including its noise). Before trusting this model's\n")
  cat("predictions on new data, strongly consider:\n")
  cat("  - xgb.cv(...) to cross-validate and see whether performance on held-out\n")
  cat("    folds keeps improving or has started getting WORSE (a sign of overfitting)\n")
  cat("  - early_stopping_rounds with a watchlist, to automatically stop adding trees\n")
  cat("    once held-out performance stops improving\n")
  cat("  - comparing this model's held-out performance against the simpler models in\n")
  cat("    this workflow (logistic regression, random forest) — the extra complexity\n")
  cat("    of boosting is only worth it if it demonstrably predicts better on data it\n")
  cat("    hasn't seen\n")
  cat("========================================================\n\n")
  
  
  # ==========================================================
  # FEATURE IMPORTANCE
  # ==========================================================
  fnames <- feature_names
  if (is.null(fnames)) {
    fnames <- tryCatch(model$feature_names, error = function(e) NULL)
  }
  
  importance_dt <- tryCatch(
    if (is.null(fnames)) xgb.importance(model = model) else xgb.importance(feature_names = fnames, model = model),
    error = function(e) NULL
  )
  
  if (is.null(importance_dt) || nrow(importance_dt) == 0) {
    cat(bold("WARNING:"), "could not compute feature importance for this model (or it came back\n")
    cat("empty). If this model's booster doesn't have feature names stored internally,\n")
    cat("pass them explicitly: explain_xgb(model, feature_names = colnames(x_xgb)).\n")
    cat("========================================================\n")
    return(invisible(NULL))
  }
  
  cat("Feature importance:\n")
  cat("--------------------------------------------------------\n")
  cat("Gain    — how much each feature's splits improved prediction accuracy, summed\n")
  cat("          across every tree. This is the MAIN measure of importance to look at.\n")
  cat("Cover   — roughly, how many patients passed through splits on this feature.\n")
  cat("Frequency — how many times this feature was used across all trees/splits.\n\n")
  
  print_cols <- intersect(c("Feature", "Gain", "Cover", "Frequency"), colnames(importance_dt))
  print(as.data.frame(importance_dt[, ..print_cols]))
  
  cat("\nIn plain English: features with high Gain contributed the most to this model's\n")
  cat("predictions. As with random forest importance, this tells you WHICH features\n")
  cat("mattered for prediction — NOT the direction of their effect (does higher age\n")
  cat("increase or decrease risk?) or a directly interpretable effect size like an odds\n")
  cat("ratio. For direction, consider SHAP values (SHAPforxgboost package) or partial\n")
  cat("dependence plots, which this function does not compute.\n\n")
  
  top_gain_pct <- if ("Gain" %in% colnames(importance_dt)) round(100 * importance_dt$Gain[1], 1) else NA
  if (!is.na(top_gain_pct) && top_gain_pct > 50) {
    cat(sprintf("%s a single feature ('%s') accounts for over half of the model's total\n",
                bold("NOTE:"), importance_dt$Feature[1]))
    cat("Gain. This model is leaning heavily on one predictor — worth checking whether\n")
    cat("that's expected/clinically sensible, or a sign something else went wrong (e.g.\n")
    cat("a leaked/proxy variable).\n\n")
  }
  
  cat("--------------------------------------------------------\n")
  cat(bold("CAVEATS — is XGBoost the right tool here?\n"))
  cat("Good for: squeezing out maximum prediction accuracy, especially with larger\n")
  cat("datasets and complex/non-linear relationships; handles missing data natively;\n")
  cat("widely used and well-validated across many prediction competitions/industries.\n")
  cat("Weaker for: interpretability (even more of a 'black box' than random forest —\n")
  cat("there's no single tree to trace); OVERFITTING risk is real and requires active\n")
  cat("management (see the critical caveat above) — unlike random forest, boosting does\n")
  cat("NOT have a built-in safeguard against fitting training-data noise; requires more\n")
  cat("careful hyperparameter tuning (depth, learning rate, rounds, regularization) to\n")
  cat("get right, and results can be sensitive to these choices.\n")
  cat("No p-values or confidence intervals are produced by this method.\n")
  cat("========================================================\n")
  
  if (isTRUE(plot_importance)) {
    tryCatch({
      xgb.plot.importance(importance_dt, top_n = min(15, nrow(importance_dt)))
    }, error = function(e) {
      warning(sprintf("Could not draw the importance plot: %s", conditionMessage(e)))
    })
  }
  
  invisible(importance_dt)
}


# ============================================================
# USAGE EXAMPLES (using the cohort dataset from earlier)
# ============================================================
library(rpart); library(rpart.plot)
model_tree <- rpart(event_f ~ age + bmi + comorbidity_score + sex + treatment,
                     data = cohort, method = "class")
tree_results <- explain_rpart(model_tree)
tree_results <- explain_rpart(model_tree, plot_tree = TRUE)

library(randomForest)
model_rf <- randomForest(event_f ~ age + bmi + comorbidity_score + sex + treatment,
                          data = cohort, ntree = 500, importance = TRUE)
rf_results <- explain_rf(model_rf)
rf_results <- explain_rf(model_rf, plot_importance = TRUE)

library(xgboost)
x_xgb <- model.matrix(event_f ~ age + bmi + comorbidity_score + sex + treatment,
                       data = cohort)[, -1]
y_xgb <- cohort$event
dtrain <- xgb.DMatrix(data = x_xgb, label = y_xgb)
model_xgb <- xgb.train(data = dtrain, nrounds = 100,
                       params = list(objective = "binary:logistic", max_depth = 4, eta = 0.1),
                       verbose = 0)
xgb_results <- explain_xgb(model_xgb)
xgb_results <- explain_xgb(model_xgb, plot_importance = TRUE)

# ============================================================
# EDGE CASES THIS FILE HAS BEEN SPECIFICALLY CHECKED AGAINST
# ============================================================
# explain_rpart():
#  1. A tree with zero splits (root node only)                  -> flagged, function stops cleanly
#     rather than erroring on downstream table access
#  2. model$method not readable                                 -> flagged, function still proceeds
#  3. No variable.importance available (very simple tree)       -> handled, section skipped gracefully
#  4. Variables available but never used in any split           -> explicitly listed and explained
#  5. 1-SE-rule simplification landing on the SAME row as the
#     minimum-xerror row                                        -> only shown if strictly simpler,
#     avoids a redundant/confusing suggestion
#  6. plot_tree = TRUE failing (e.g. graphics device issue)      -> caught, warns instead of crashing
#
# explain_rf():
#  7. Regression forest (numeric outcome) vs classification forest -> detected via model$type,
#     completely different performance section shown for each
#  8. Fit without importance = TRUE                              -> flagged, importance section
#     skipped safely instead of erroring on a NULL access
#  9. Severely imbalanced classification outcome                 -> flagged, points to per-class
#     error rates instead of misleading overall accuracy
# 10. Error curve check on very small ntree (<20 trees)           -> skipped (the check requires
#     enough trees to compare a meaningful early/late window)
# 11. Only one of MeanDecreaseAccuracy/MeanDecreaseGini present
#     (importance computed with a subset of options)             -> each handled independently,
#     doesn't assume both exist
# 12. plot_importance = TRUE without importance available         -> warns instead of erroring
#
# explain_xgb():
# 13. No feature names stored on the booster AND none supplied    -> clear, actionable warning
#     telling you exactly what argument to pass, instead of a cryptic xgb.importance() error
# 14. xgb.importance() returning zero rows                        -> flagged, function stops cleanly
# 15. Missing params (e.g. max_depth/eta not explicitly set,
#     defaults used instead)                                     -> defensive NULL checks, shown as
#     "default" rather than crashing on a missing list element
# 16. A single feature dominating total Gain (>50%)               -> explicitly flagged as worth a
#     sanity check (possible data leakage or proxy variable)
# 17. plot_importance = TRUE failing                              -> caught, warns instead of crashing
# 18. Regression/ranking/survival objectives (not just binary
#     classification)                                             -> objective is read and displayed
#     directly from the model rather than assumed to be binary
# 19. model$niter, model$method, or model$type being NULL rather than erroring when accessed
#     (R's `$` returns NULL for a missing list element instead of throwing) -> every such access is
#     explicitly checked for NULL and converted to NA before being used in an is.na()/if() test,
#     since R's is.na(NULL) returns a zero-length result that crashes a subsequent if()/||
#     condition rather than behaving like a normal NA — this was caught and fixed during review
#     in three separate places in this file
# ============================================================