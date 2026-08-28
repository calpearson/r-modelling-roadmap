# 09_tidymodels

# ============================================================
# 9. TIDYMODELS — outcome: event_f (same predictors as section 8)
# ============================================================
library(tidymodels)
library(tidyverse)
library(ranger)

split <- initial_split(cohort, prop = 0.8, strata = event_f)
train <- training(split)
test  <- testing(split)

rf_spec <- rand_forest(trees = 500) %>%
  set_engine("ranger") %>%
  set_mode("classification")

rf_wf <- workflow() %>%
  add_formula(event_f ~ age + bmi + comorbidity_score + sex + treatment) %>%
  add_model(rf_spec)

model_tidy_rf <- rf_wf %>% fit(data = train)

out_tidy_predictions <- predict(model_tidy_rf, test) %>% bind_cols(test)
out_tidy_metrics     <- out_tidy_predictions %>% metrics(truth = event_f, estimate = .pred_class)
out_tidy_predictions
out_tidy_metrics

folds <- vfold_cv(train, v = 5)
out_tidy_cv_metrics <- fit_resamples(rf_wf, folds) %>% collect_metrics()
out_tidy_cv_metrics






# ============================================================
# explain_tidymodels.R — plain-language interpretation of a
# tidymodels workflow: initial_split(), a fitted workflow,
# test-set predictions/metrics, and cross-validation results.
#
# Four functions, one per object you create along the way:
#   explain_split(split)                         — initial_split() object
#   explain_rf_workflow(model_fit)                — the fitted workflow
#   explain_test_performance(predictions, metrics, model_fit = NULL, test_data = NULL)
#   explain_cv_performance(cv_metrics, test_metrics = NULL)
#
# Beyond plain-English explanation, this file runs ACTUAL
# STATISTICAL TESTS of whether the analysis approach is well
# suited to your data (stratification balance check, a
# no-information-rate test, cross-validation stability check) —
# not just descriptive commentary — and produces ggplot2 figures
# with full titles for each key output, since ggplot2 is already
# a dependency of this specific workflow (via tidyverse).
# ============================================================

library(tidymodels)
library(tidyverse)

# ---- Shared helper ----
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) crayon::bold(x) else x
}


# ============================================================
# 1. explain_split() — the train/test split (initial_split)
# ============================================================
explain_split <- function(split) {
  
  if (!inherits(split, "rsplit")) {
    stop("explain_split() expects an object from rsample::initial_split().")
  }
  
  train_data <- tryCatch(training(split), error = function(e) NULL)
  test_data  <- tryCatch(testing(split), error = function(e) NULL)
  
  if (is.null(train_data) || is.null(test_data)) {
    stop("Could not extract training()/testing() data from this split object.")
  }
  
  cat("========================================================\n")
  cat("TRAIN / TEST SPLIT\n")
  cat("========================================================\n\n")
  
  cat("What this does, in plain English:\n")
  cat("The data is randomly divided into a TRAINING set (used to fit the model) and a\n")
  cat("TEST set (held back, completely unseen by the model during fitting, used only\n")
  cat("to check how well it performs on 'new' patients). This guards against fooling\n")
  cat("yourself by checking a model's accuracy on the same data it learned from.\n\n")
  
  n_train <- nrow(train_data); n_test <- nrow(test_data); n_total <- n_train + n_test
  cat(sprintf("Training set: %s patients (%.0f%%)\n", format(n_train, big.mark = ","), 100 * n_train / n_total))
  cat(sprintf("Test set:     %s patients (%.0f%%)\n\n", format(n_test, big.mark = ","), 100 * n_test / n_total))
  
  # ==========================================================
  # ATTRIBUTABLE TEST 1: did stratification actually balance the outcome
  # between train and test? (a real chi-squared test, not just eyeballing)
  # ==========================================================
  strata_var <- tryCatch(split$strata, error = function(e) NULL)
  
  cat("Did stratified splitting actually balance the outcome between the two sets?\n")
  cat("--------------------------------------------------------\n")
  
  outcome_candidates <- intersect(names(train_data), names(test_data))
  factor_candidates <- outcome_candidates[sapply(train_data[outcome_candidates], function(x) is.factor(x) && nlevels(x) <= 6)]
  
  if (length(factor_candidates) == 0) {
    cat("Could not identify a categorical outcome column to check balance on.\n\n")
  } else {
    check_var <- factor_candidates[1]
    if (!is.null(strata_var) && strata_var %in% factor_candidates) check_var <- strata_var
    
    tab <- rbind(
      Train = table(train_data[[check_var]]),
      Test  = table(test_data[[check_var]])
    )
    cat(sprintf("Class balance for '%s':\n\n", check_var))
    print(tab)
    prop_tab <- prop.table(tab, margin = 1)
    cat("\n")
    print(round(prop_tab, 3))
    
    test_result <- tryCatch(suppressWarnings(chisq.test(tab)), error = function(e) NULL)
    if (!is.null(test_result)) {
      cat(sprintf("\nChi-squared test for a difference in '%s' distribution between train and\n", check_var))
      cat(sprintf("test sets: chi-squared = %.2f, p = %.3f\n", test_result$statistic, test_result$p.value))
      if (test_result$p.value < 0.05) {
        cat(bold("\nCAUTION:"), "the outcome distribution differs significantly between your\n")
        cat("train and test sets (p < 0.05). If you intended to stratify on this variable,\n")
        cat("double check strata = ... in initial_split() names the right column — this\n")
        cat("result suggests the split may not be as balanced as intended.\n")
      } else {
        cat("\nNo significant difference detected (p >= 0.05) — train and test sets have a\n")
        cat("similar outcome distribution, consistent with successful stratification.\n")
      }
    }
  }
  
  # ==========================================================
  # ATTRIBUTABLE TEST 2: is the TEST set big enough for stable performance
  # estimates? (rule-of-thumb minimum event count check)
  # ==========================================================
  cat("\nIs the test set large enough for a reliable performance estimate?\n")
  cat("--------------------------------------------------------\n")
  
  if (length(factor_candidates) > 0) {
    check_var <- factor_candidates[1]
    min_class_n <- min(table(test_data[[check_var]]))
    cat(sprintf("Smallest outcome class in the TEST set: %d patients.\n", min_class_n))
    if (min_class_n < 30) {
      cat(bold("CAUTION:"), sprintf("only %d patients in the smallest test-set class. A single\n", min_class_n))
      cat("test-set performance metric (accuracy, kappa, etc.) computed from this few\n")
      cat("events can swing substantially just from which patients happened to land in\n")
      cat("the test set — treat any single test-set metric with real caution, and lean\n")
      cat("more heavily on the cross-validated results (see explain_cv_performance()),\n")
      cat("which average over many different splits rather than relying on just one.\n")
    } else {
      cat("This looks like a reasonable number of events for a moderately stable\n")
      cat("test-set performance estimate, though cross-validation (below) still gives\n")
      cat("you a more robust picture by using ALL the data across multiple splits.\n")
    }
  }
  
  cat("\n--------------------------------------------------------\n")
  cat(bold("Is a single train/test split the right approach here?\n"))
  cat("A single split is fast and simple, and fine for a first look or for VERY large\n")
  cat("datasets. But its performance estimate depends on which patients happened to\n")
  cat("land in the test set by chance — a different random split could give a\n")
  cat("meaningfully different number. K-FOLD CROSS-VALIDATION (see\n")
  cat("explain_cv_performance() below) fixes this by repeating the process across\n")
  cat("multiple splits and averaging — generally preferred unless your dataset is so\n")
  cat("large that a single split is already very stable, or your workflow needs the\n")
  cat("speed of fitting only once.\n")
  cat("========================================================\n")
  
  # ---- Figure: class balance bar chart, train vs test ----
  if (length(factor_candidates) > 0) {
    check_var <- factor_candidates[1]
    plot_df <- bind_rows(
      train_data %>% count(.data[[check_var]]) %>% mutate(set = "Train"),
      test_data  %>% count(.data[[check_var]]) %>% mutate(set = "Test")
    )
    p <- ggplot(plot_df, aes(x = .data[[check_var]], y = n, fill = set)) +
      geom_col(position = "dodge") +
      labs(title = "Outcome class balance: Training set vs. Test set",
           subtitle = sprintf("Checking that stratified splitting kept '%s' balanced across both sets", check_var),
           x = check_var, y = "Number of patients", fill = "Data set") +
      theme_minimal(base_size = 13)
    print(p)
  }
  
  invisible(list(n_train = n_train, n_test = n_test))
}


# ============================================================
# 2. explain_rf_workflow() — the fitted tidymodels workflow
# ============================================================
explain_rf_workflow <- function(model_fit) {
  
  if (!inherits(model_fit, "workflow")) {
    stop("explain_rf_workflow() expects a fitted tidymodels workflow object (from workflow() %>% fit()).")
  }
  if (!tryCatch(model_fit$trained, error = function(e) FALSE)) {
    stop("This workflow does not appear to be fitted yet — call fit(workflow, data = ...) first.")
  }
  
  cat("========================================================\n")
  cat("FITTED MODEL (tidymodels workflow)\n")
  cat("========================================================\n\n")
  
  spec <- tryCatch(extract_spec_parsnip(model_fit), error = function(e) NULL)
  engine_fit <- tryCatch(extract_fit_engine(model_fit), error = function(e) NULL)
  
  if (!is.null(spec)) {
    cat(sprintf("Model type: %s | Engine: %s | Mode: %s\n\n",
                bold(class(spec)[1]), bold(spec$engine), bold(spec$mode)))
  }
  
  cat("A tidymodels 'workflow' bundles together the FORMULA (which variables predict\n")
  cat("the outcome) and the MODEL SPECIFICATION (which algorithm, and its settings)\n")
  cat("into one object — this makes it easy to swap in a different algorithm later\n")
  cat("while keeping the rest of your pipeline (splitting, resampling, metrics)\n")
  cat("unchanged, which is the main advantage of this framework over calling\n")
  cat("randomForest()/glm()/etc. directly.\n\n")
  
  if (!is.null(engine_fit) && inherits(engine_fit, "ranger")) {
    cat("Underlying engine details (ranger — a fast Random Forest implementation):\n")
    cat("--------------------------------------------------------\n")
    cat(sprintf("  Number of trees: %s\n", engine_fit$num.trees))
    cat(sprintf("  Variables tried at each split (mtry): %s\n", engine_fit$mtry))
    cat(sprintf("  Out-of-Bag prediction error: %.1f%%\n\n",
                ifelse(is.null(engine_fit$prediction.error), NA, engine_fit$prediction.error * 100)))
    
    has_importance <- !is.null(engine_fit$variable.importance) && length(engine_fit$variable.importance) > 0
    
    if (has_importance) {
      cat("Variable importance:\n")
      imp <- sort(engine_fit$variable.importance, decreasing = TRUE)
      imp_df <- tibble(variable = names(imp), importance = as.numeric(imp))
      print(imp_df)
      cat("\nIn plain English: higher values mean the variable contributed more to this\n")
      cat("model's predictions. As with any random forest, this tells you WHICH\n")
      cat("variables mattered, not the DIRECTION of their effect.\n\n")
      
      p <- ggplot(imp_df, aes(x = reorder(variable, importance), y = importance)) +
        geom_col(fill = "steelblue") +
        coord_flip() +
        labs(title = "Variable importance — Random Forest (ranger engine)",
             subtitle = "Higher = more useful for prediction; does NOT indicate direction of effect",
             x = NULL, y = "Importance") +
        theme_minimal(base_size = 13)
      print(p)
    } else {
      cat(bold("SUGGESTED IMPROVEMENT:"), "this model was fit without variable importance\n")
      cat("enabled, so no importance figure can be shown. To enable it, refit with:\n\n")
      cat("  rf_spec <- rand_forest(trees = 500) %>%\n")
      cat("    set_engine(\"ranger\", importance = \"permutation\") %>%\n")
      cat("    set_mode(\"classification\")\n\n")
      cat("('permutation' importance is generally preferred over ranger's alternative,\n")
      cat("'impurity', which is biased toward continuous/high-cardinality variables —\n")
      cat("same caveat as randomForest's Gini importance.)\n\n")
    }
  } else if (!is.null(engine_fit)) {
    cat(sprintf("Engine object is of class '%s' — this function's detailed engine-specific\n",
                class(engine_fit)[1]))
    cat("summary only covers the 'ranger' engine; showing the raw fitted object instead:\n\n")
    print(engine_fit)
  }
  
  cat("--------------------------------------------------------\n")
  cat(bold("No p-values or confidence intervals are produced by this method — expected\n"))
  cat(bold("for a random forest, not a gap in this function; see explain_test_performance()\n"))
  cat(bold("and explain_cv_performance() for how to assess this model's actual PERFORMANCE.\n"))
  cat("========================================================\n")
  
  invisible(list(spec = spec, engine_fit = engine_fit))
}


# ============================================================
# 3. explain_test_performance() — held-out test set predictions + metrics
# ============================================================
explain_test_performance <- function(predictions, metrics, model_fit = NULL, test_data = NULL, truth_col = NULL) {
  
  if (!is.data.frame(predictions) || !(".pred_class" %in% names(predictions))) {
    stop("explain_test_performance() expects `predictions` to be a data frame containing a `.pred_class` column, as produced by predict(workflow, new_data) %>% bind_cols(new_data).")
  }
  
  if (is.null(truth_col)) {
    candidates <- names(predictions)[sapply(predictions, is.factor)]
    candidates <- setdiff(candidates, ".pred_class")
    if (length(candidates) == 1) {
      truth_col <- candidates
    } else if (length(candidates) > 1) {
      stop(sprintf("Could not automatically determine the truth column — multiple factor columns found (%s). Specify explicitly: explain_test_performance(predictions, metrics, truth_col = \"your_outcome_column\")",
                   paste(candidates, collapse = ", ")))
    } else {
      stop("Could not find a factor 'truth' column in `predictions` to compare against .pred_class. Specify truth_col explicitly.")
    }
  }
  
  cat("========================================================\n")
  cat("TEST-SET PERFORMANCE\n")
  cat("========================================================\n\n")
  
  n_test <- nrow(predictions)
  cat(sprintf("Evaluated on %s patients the model did NOT see during training.\n\n",
              bold(format(n_test, big.mark = ","))))
  
  cat("Performance metrics:\n")
  cat("--------------------------------------------------------\n")
  print(as.data.frame(metrics))
  cat("\n")
  
  metric_explanations <- list(
    accuracy = "the proportion of patients whose outcome the model predicted correctly, overall.",
    kap = "Cohen's Kappa — agreement between predicted and actual outcomes, ADJUSTED for the agreement you'd expect from chance alone. 0 = no better than chance; 1 = perfect agreement. Rough guide: <0.2 slight, 0.2-0.4 fair, 0.4-0.6 moderate, 0.6-0.8 substantial, >0.8 almost perfect agreement.",
    sens = "sensitivity (recall) — of patients who ACTUALLY had the event, what proportion did the model correctly flag?",
    spec = "specificity — of patients who did NOT have the event, what proportion did the model correctly identify as not having it?",
    precision = "of patients the model PREDICTED would have the event, what proportion actually did?",
    f_meas = "F1 score — a single number balancing precision and recall/sensitivity together.",
    roc_auc = "Area under the ROC curve — how well the model RANKS patients by risk, across every possible decision threshold, not just the default 0.5 cutoff. 0.5 = no better than random; 1.0 = perfect ranking."
  )
  
  for (i in seq_len(nrow(metrics))) {
    mname <- as.character(metrics$.metric[i])
    mval  <- metrics$.estimate[i]
    if (mname %in% names(metric_explanations)) {
      cat(sprintf("> %s = %.3f\n  %s\n\n", mname, mval, metric_explanations[[mname]]))
    }
  }
  
  cat("Confusion matrix:\n")
  cat("--------------------------------------------------------\n")
  cm <- tryCatch(conf_mat(predictions, truth = !!sym(truth_col), estimate = .pred_class),
                 error = function(e) NULL)
  if (!is.null(cm)) {
    print(cm)
    cat("\nRows = predicted outcome, columns = actual outcome (check your yardstick\n")
    cat("version's orientation if unsure — run `print(cm)` yourself to confirm).\n\n")
    
    p_cm <- autoplot(cm, type = "heatmap") +
      labs(title = "Confusion Matrix - Test Set Predictions vs. Actual Outcome",
           subtitle = sprintf("n = %d patients held out from model training", n_test)) +
      theme_minimal(base_size = 13)
    print(p_cm)
  }
  
  cat("--------------------------------------------------------\n")
  cat(bold("Is this model actually better than a trivial guess?\n"))
  cat("--------------------------------------------------------\n")
  
  truth_vec <- predictions[[truth_col]]
  class_props <- prop.table(table(truth_vec))
  no_info_rate <- max(class_props)
  majority_class <- names(class_props)[which.max(class_props)]
  
  n_correct <- sum(predictions$.pred_class == truth_vec)
  observed_accuracy <- n_correct / n_test
  
  cat(sprintf("'No-information rate': if you predicted every single patient was '%s'\n", majority_class))
  cat(sprintf("(the majority class) with no model at all, you'd be right %.1f%% of the time.\n\n", no_info_rate * 100))
  cat(sprintf("This model's actual accuracy: %.1f%%\n\n", observed_accuracy * 100))
  
  test_vs_baseline <- tryCatch(
    binom.test(n_correct, n_test, p = no_info_rate, alternative = "greater"),
    error = function(e) NULL
  )
  
  if (!is.null(test_vs_baseline)) {
    cat(sprintf("One-sided binomial test (is accuracy > no-information rate?): p = %.4f\n",
                test_vs_baseline$p.value))
    if (test_vs_baseline$p.value < 0.05) {
      cat(bold("\nGOOD SIGN:"), "this model performs SIGNIFICANTLY better than the trivial\n")
      cat("'always guess the majority class' baseline — it's learning a real signal, not\n")
      cat("just exploiting class imbalance.\n")
    } else {
      cat(bold("\nCAUTION:"), "this model is NOT significantly better than just always guessing\n")
      cat(sprintf("'%s' for every patient (p >= 0.05). ", majority_class))
      cat("With imbalanced classes, accuracy can look\n")
      cat("deceptively high while adding little real predictive value — check\n")
      cat("sensitivity/specificity for the minority class specifically, and consider\n")
      cat("whether this model (or this metric) is actually the right choice here.\n")
    }
  }
  
  p_baseline <- ggplot(
    tibble(type = c("No-information\n(always guess majority)", "This model"),
           acc = c(no_info_rate, observed_accuracy)),
    aes(x = type, y = acc, fill = type)
  ) +
    geom_col(width = 0.5, show.legend = FALSE) +
    geom_hline(yintercept = no_info_rate, linetype = "dashed", color = "red") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = "Model Accuracy vs. Trivial 'Always Guess Majority Class' Baseline",
         subtitle = "Dashed line = no-information rate. Model should clear this line by a meaningful margin.",
         x = NULL, y = "Accuracy") +
    theme_minimal(base_size = 13)
  print(p_baseline)
  
  if (!is.null(model_fit) && !is.null(test_data)) {
    cat("\n--------------------------------------------------------\n")
    cat("ROC curve (using the model + test data you supplied):\n")
    cat("--------------------------------------------------------\n")
    
    prob_preds <- tryCatch(predict(model_fit, test_data, type = "prob"), error = function(e) NULL)
    if (!is.null(prob_preds)) {
      prob_col <- names(prob_preds)[1]
      roc_data <- bind_cols(prob_preds, test_data)
      
      auc_val <- tryCatch(
        roc_auc(roc_data, truth = !!sym(truth_col), !!sym(prob_col))$.estimate,
        error = function(e) NA
      )
      cat(sprintf("AUC = %.3f\n\n", auc_val))
      
      p_roc <- tryCatch({
        roc_curve(roc_data, truth = !!sym(truth_col), !!sym(prob_col)) %>%
          autoplot() +
          labs(title = "ROC Curve - Test Set",
               subtitle = sprintf("AUC = %.3f (0.5 = no better than chance, 1.0 = perfect)", auc_val)) +
          theme_minimal(base_size = 13)
      }, error = function(e) NULL)
      if (!is.null(p_roc)) print(p_roc)
    } else {
      cat("Could not generate class probabilities from the supplied model/data —\n")
      cat("skipping the ROC curve.\n")
    }
  } else {
    cat("\n")
    cat(bold("SUGGESTED IMPROVEMENT:"), "your original predict() call only generated hard\n")
    cat("class labels (.pred_class), not probabilities — this means no ROC curve/AUC can\n")
    cat("be shown, and you're implicitly using a fixed 0.5 probability threshold for\n")
    cat("classification, which may not be optimal for your specific use case. To get\n")
    cat("probabilities too:\n\n")
    cat("  test_probs <- predict(model_tidy_rf, test, type = \"prob\") %>% bind_cols(test)\n\n")
    cat("Then pass `model_fit = model_tidy_rf, test_data = test` to this function for an\n")
    cat("automatic ROC curve and AUC.\n")
  }
  
  cat("========================================================\n")
  
  invisible(list(observed_accuracy = observed_accuracy, no_info_rate = no_info_rate,
                 beats_baseline_p = if (!is.null(test_vs_baseline)) test_vs_baseline$p.value else NA))
}


# ============================================================
# 4. explain_cv_performance() — k-fold cross-validation results
# ============================================================
explain_cv_performance <- function(cv_metrics, test_metrics = NULL) {
  
  if (!is.data.frame(cv_metrics) || !all(c(".metric", "mean") %in% names(cv_metrics))) {
    stop("explain_cv_performance() expects the data frame produced by fit_resamples(workflow, folds) %>% collect_metrics().")
  }
  
  cat("========================================================\n")
  cat("CROSS-VALIDATION PERFORMANCE\n")
  cat("========================================================\n\n")
  
  cat("What this does, in plain English:\n")
  cat("Instead of relying on ONE random train/test split (which could happen to be\n")
  cat("lucky or unlucky), the training data is divided into several 'folds'. The model\n")
  cat("is fit multiple times, each time holding out a DIFFERENT fold to test on and\n")
  cat("training on the rest — every patient in the training data gets used for both\n")
  cat("training AND testing at some point, just never in the same round. This gives\n")
  cat("you a DISTRIBUTION of performance estimates rather than a single number,\n")
  cat("showing you how much performance actually varies depending on which patients\n")
  cat("happen to be held out.\n\n")
  
  n_folds <- tryCatch(unique(cv_metrics$n), error = function(e) NA)
  if (length(n_folds) == 1 && !is.na(n_folds)) {
    cat(sprintf("Number of folds used: %s\n\n", n_folds))
  }
  
  cat("Results (mean across folds, with standard error):\n")
  cat("--------------------------------------------------------\n")
  print(as.data.frame(cv_metrics[, intersect(c(".metric", "mean", "n", "std_err"), names(cv_metrics))]))
  cat("\n")
  
  cat("How stable are these estimates?\n")
  cat("--------------------------------------------------------\n")
  
  unstable_metrics <- character(0)
  
  for (i in seq_len(nrow(cv_metrics))) {
    mname <- as.character(cv_metrics$.metric[i])
    mean_val <- cv_metrics$mean[i]
    se_val   <- if ("std_err" %in% names(cv_metrics)) cv_metrics$std_err[i] else NA
    
    if (!is.na(se_val) && mean_val != 0) {
      cv_pct <- 100 * se_val / abs(mean_val)
      cat(sprintf("> %s: mean = %.3f, SE = %.3f (SE is %.1f%% of the mean)\n", mname, mean_val, se_val, cv_pct))
      if (cv_pct > 10) {
        unstable_metrics <- c(unstable_metrics, mname)
      }
    }
  }
  
  if (length(unstable_metrics) > 0) {
    cat(sprintf("\n%s the fold-to-fold variability for %s is fairly large relative to its\n",
                bold("CAUTION:"), paste(unstable_metrics, collapse = ", ")))
    cat("mean value — this metric's estimate isn't very PRECISE yet. Consider using\n")
    cat("more folds (e.g. v = 10), or REPEATED cross-validation for more stable\n")
    cat("estimates:\n\n")
    cat("  folds_repeated <- vfold_cv(train, v = 5, repeats = 5)\n\n")
  } else {
    cat("\nThese estimates look reasonably stable (SE small relative to the mean).\n\n")
  }
  
  if (!is.null(test_metrics)) {
    cat("--------------------------------------------------------\n")
    cat("Comparing to your single train/test-split result:\n")
    cat("--------------------------------------------------------\n")
    for (i in seq_len(nrow(cv_metrics))) {
      mname <- as.character(cv_metrics$.metric[i])
      cv_mean <- cv_metrics$mean[i]
      cv_se   <- if ("std_err" %in% names(cv_metrics)) cv_metrics$std_err[i] else NA
      
      single_row <- test_metrics[as.character(test_metrics$.metric) == mname, ]
      if (nrow(single_row) == 1) {
        single_val <- single_row$.estimate[1]
        cat(sprintf("> %s: single test-set = %.3f | cross-validated = %.3f (SE %.3f)\n",
                    mname, single_val, cv_mean, ifelse(is.na(cv_se), NA, cv_se)))
        if (!is.na(cv_se) && abs(single_val - cv_mean) > 2 * cv_se) {
          cat("  This single test-set value is notably DIFFERENT from the cross-validated\n")
          cat("  average (more than ~2 SEs away) — a sign that your one train/test split may\n")
          cat("  have been somewhat lucky or unlucky. Trust the cross-validated number more.\n")
        }
      }
    }
    cat("\n")
  }
  
  if ("std_err" %in% names(cv_metrics)) {
    plot_df <- cv_metrics %>%
      mutate(lower = mean - std_err, upper = mean + std_err)
    
    p <- ggplot(plot_df, aes(x = .metric, y = mean)) +
      geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.8, color = "darkblue") +
      coord_flip() +
      labs(title = "Cross-Validated Performance Estimates (mean +/- 1 SE)",
           subtitle = sprintf("Averaged across %s resampling folds - error bars show fold-to-fold variability",
                              ifelse(length(n_folds) == 1 && !is.na(n_folds), n_folds, "multiple")),
           x = NULL, y = "Metric value") +
      theme_minimal(base_size = 13)
    print(p)
  }
  
  cat("--------------------------------------------------------\n")
  cat(bold("Is cross-validation the right approach here - and is this CV setup adequate?\n"))
  cat("Cross-validation is almost always preferable to a single train/test split when\n")
  cat("you can afford to fit the model multiple times (it's more computationally\n")
  cat("expensive, but generally worth it for a more honest, stable performance\n")
  cat("estimate). Watch out for: (1) if your outcome is rare, make sure folds are\n")
  cat("STRATIFIED (vfold_cv(..., strata = outcome)) so rare-outcome patients aren't\n")
  cat("concentrated in just a few folds; (2) if any wide SEs were flagged above,\n")
  cat("consider more folds or repeats for a more precise estimate; (3) cross-validation\n")
  cat("alone does NOT replace a genuinely held-out final test set if you're also using\n")
  cat("the CV results to TUNE hyperparameters - in that case you'd want a further,\n")
  cat("untouched test set for a final, unbiased performance check.\n")
  cat("========================================================\n")
  
  invisible(list(unstable_metrics = unstable_metrics))
}


# ============================================================
# USAGE EXAMPLES (using the cohort dataset and workflow from earlier)
# ============================================================
library(tidymodels); library(tidyverse); library(ranger)

split_results <- explain_split(split)

workflow_results <- explain_rf_workflow(model_tidy_rf)

test_results <- explain_test_performance(out_tidy_predictions, out_tidy_metrics)
# ...or with the ROC curve, if you generate probability predictions too:
# test_probs <- predict(model_tidy_rf, test, type = "prob") %>% bind_cols(test)
test_results <- explain_test_performance(out_tidy_predictions, out_tidy_metrics,
                                          model_fit = model_tidy_rf, test_data = test)

cv_results <- explain_cv_performance(out_tidy_cv_metrics, test_metrics = out_tidy_metrics)

# ============================================================
# EDGE CASES THIS FILE HAS BEEN SPECIFICALLY CHECKED AGAINST
# ============================================================
# explain_split():
#  1. No categorical/factor column found to check balance on        -> section skipped gracefully
#  2. split$strata missing/NULL (unstratified split)                -> falls back to the first
#     detected factor column instead of assuming a specific name
#  3. chisq.test() failing (e.g. expected cell counts too low)      -> caught, section skipped
#     rather than erroring the whole function
#
# explain_rf_workflow():
#  4. Workflow not yet fitted (fit() not called)                    -> stops with a clear, specific
#     message rather than erroring deep inside extract_fit_engine()
#  5. Engine other than 'ranger' (e.g. glmnet, xgboost via parsnip)  -> generic fallback: prints the
#     raw engine object instead of assuming ranger-specific fields exist
#  6. ranger fit without importance enabled                          -> explicitly detected, given a
#     specific "how to fix" suggestion rather than silently omitting the figure
#
# explain_test_performance():
#  7. `predictions` missing a `.pred_class` column                   -> stops immediately with a
#     specific, actionable message
#  8. Truth column not supplied and ambiguous (multiple factor columns) -> stops, asks you to specify
#     truth_col explicitly rather than guessing wrong
#  9. Truth column not supplied and NONE found                        -> stops with a clear message
# 10. conf_mat() failing for any reason                               -> caught, whole function
#     continues past the confusion-matrix section rather than halting entirely
# 11. binom.test() failing (e.g. n_test = 0)                          -> caught, that section skipped
# 12. `model_fit`/`test_data` not supplied for ROC curve              -> skipped gracefully, with a
#     specific suggestion for how to enable it, instead of erroring on missing args
# 13. predict(..., type = "prob") failing on the supplied model/data  -> caught, ROC section skipped
#     with an explanation instead of crashing
#
# explain_cv_performance():
# 14. `cv_metrics` missing expected columns (.metric/mean)            -> stops immediately with a
#     specific, actionable message
# 15. `std_err` column absent (e.g. older tidymodels versions/some resampling setups) -> stability
#     check and point-range figure both skipped gracefully rather than erroring
# 16. `test_metrics` not supplied for comparison                      -> comparison section skipped
#     cleanly, rest of the function still runs
# 17. A CV metric name with no matching row in `test_metrics`         -> that specific comparison
#     line is skipped, doesn't halt the loop for other metrics
# 18. mean value of exactly 0 for a metric (division-by-zero risk in % SE calc) -> guarded explicitly
# ============================================================
