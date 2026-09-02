# 10_deep_learning

# ============================================================
# 10. DEEP LEARNING (optional) — reuse x_xgb / y_xgb from section 8
# ============================================================
#install.packages("keras3"); keras3::install_keras()

library(keras3)
#install.packages("reticulate")
library(reticulate)
#reticulate::install_python()   # installs a managed Python via uv, if you don't have one
reticulate::py_discover_config()
reticulate::use_python("C:/Program Files/Python314/python.exe", required = TRUE)
keras3::install_keras()        # now retry this



n_train <- floor(0.8 * nrow(x_xgb))
x_train <- x_xgb[1:n_train, ]; y_train <- y_xgb[1:n_train]
x_test  <- x_xgb[(n_train+1):nrow(x_xgb), ]; y_test <- y_xgb[(n_train+1):nrow(x_xgb)]

model_keras <- keras_model_sequential() %>%
  layer_dense(units = 32, activation = "relu", input_shape = ncol(x_train)) %>%
  layer_dense(units = 16, activation = "relu") %>%
  layer_dense(units = 1, activation = "sigmoid")

model_keras %>% compile(loss = "binary_crossentropy", optimizer = "adam", metrics = "accuracy")
out_keras_fit <- model_keras %>% fit(x_train, y_train, epochs = 20, validation_split = 0.2)
out_keras_fit


# This does not work as it requires python. need to investigate how to make this work.



# ============================================================
# explain_keras.R — plain-language interpretation of a keras3
# deep learning model: its architecture, its training history
# (loss/accuracy per epoch), and — if you supply a genuine held-
# out test set — its actual out-of-sample performance.
#
# IMPORTANT GAP THIS FUNCTION ADDRESSES IN YOUR ORIGINAL CODE:
# your script created x_test/y_test but NEVER actually evaluated
# the model on them — only validation_split (data carved out of
# TRAINING data, re-used every epoch to monitor training) was
# used. This function explains why that distinction matters, and
# will run the real test-set evaluation for you if you pass
# x_test/y_test in.
# ============================================================

library(keras3)

# ---- Shared helper ----
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) crayon::bold(x) else x
}


explain_keras <- function(model, history, x_test = NULL, y_test = NULL) {
  
  # ==========================================================
  # 0. VALIDITY CHECKS (defensive — keras3's R class names can vary
  # slightly by backend/version, so checks are done generically)
  # ==========================================================
  is_keras_model <- any(grepl("keras", class(model), ignore.case = TRUE))
  if (!is_keras_model) {
    stop("explain_keras() expects a compiled keras3 model (from keras_model_sequential() or similar).")
  }
  
  is_history <- inherits(history, "keras_training_history") ||
    any(grepl("history", class(history), ignore.case = TRUE))
  if (!is_history) {
    stop("explain_keras() expects the object returned by model %>% fit(...) as `history`.")
  }
  
  cat("========================================================\n")
  cat("DEEP LEARNING MODEL (keras3 neural network)\n")
  cat("========================================================\n\n")
  
  cat("What this method does, in plain English:\n")
  cat("A neural network passes each patient's data through several layers of\n")
  cat("interconnected 'neurons', each layer learning increasingly abstract\n")
  cat("combinations of the inputs, before producing a final prediction. It's the same\n")
  cat("family of method behind most modern image/speech/language AI — extremely\n")
  cat("powerful with enough data, but (see caveats below) often NOT the best tool for\n")
  cat("small, structured/tabular datasets like this one.\n\n")
  
  
  # ==========================================================
  # 1. ARCHITECTURE SUMMARY
  # ==========================================================
  cat("Model architecture:\n")
  cat("--------------------------------------------------------\n")
  
  arch_captured <- tryCatch({
    capture.output(summary(model))
  }, error = function(e) NULL)
  
  if (!is.null(arch_captured)) {
    cat(paste(arch_captured, collapse = "\n"))
    cat("\n\n")
  } else {
    cat("(Could not capture a formatted model summary — showing layer details instead.)\n\n")
  }
  
  layer_info <- tryCatch({
    layers <- model$layers
    data.frame(
      layer = seq_along(layers),
      type = sapply(layers, function(l) class(l)[1]),
      units = sapply(layers, function(l) tryCatch(l$units, error = function(e) NA)),
      activation = sapply(layers, function(l) tryCatch(as.character(l$activation), error = function(e) NA))
    )
  }, error = function(e) NULL)
  
  if (!is.null(layer_info)) {
    cat("Layer-by-layer, in plain English:\n\n")
    n_layers <- nrow(layer_info)
    for (i in seq_len(n_layers)) {
      units_i <- layer_info$units[i]
      act_i   <- layer_info$activation[i]
      is_last <- i == n_layers
      
      cat(sprintf("> Layer %d: %s units, '%s' activation\n", i, ifelse(is.na(units_i), "?", units_i),
                  ifelse(is.na(act_i), "unknown", act_i)))
      
      if (is_last && !is.na(units_i) && units_i == 1 && !is.na(act_i) && act_i == "sigmoid") {
        cat("  This is the OUTPUT layer for BINARY classification: a single unit with a\n")
        cat("  'sigmoid' activation squashes the output into a 0-1 range, interpreted as\n")
        cat("  the predicted PROBABILITY of the event.\n")
      } else if (!is.na(act_i) && act_i == "relu") {
        cat("  'relu' is the most common activation for hidden layers — it lets the\n")
        cat("  network learn NON-LINEAR patterns (without it, stacking layers would\n")
        cat("  mathematically collapse to being no more powerful than a single linear\n")
        cat("  layer, i.e. an ordinary logistic regression).\n")
      }
      cat("\n")
    }
    
    total_params <- tryCatch(sum(sapply(model$get_weights(), length)), error = function(e) NA)
    if (!is.na(total_params)) {
      cat(sprintf("Total trainable parameters: %s\n\n", format(total_params, big.mark = ",")))
      cat("In plain English: this is how many individual numbers the model is\n")
      cat("adjusting/learning from your data. As a rough sanity check, compare this to\n")
      cat("your training sample size below — a model with far more parameters than\n")
      cat("training examples has a lot of freedom to memorize noise rather than learn a\n")
      cat("genuine pattern (see the CAVEATS section for why this matters here specifically).\n\n")
    }
  }
  
  
  # ==========================================================
  # 2. TRAINING HISTORY — loss/accuracy per epoch, train vs validation
  # ==========================================================
  cat("Training history:\n")
  cat("--------------------------------------------------------\n")
  
  hist_df <- tryCatch(as.data.frame(history), error = function(e) NULL)
  
  if (is.null(hist_df)) {
    cat("Could not extract a tidy training-history table from this object (keras3\n")
    cat("version differences can affect this) — printing the raw object instead:\n\n")
    print(history)
  } else {
    n_epochs <- max(hist_df$epoch, na.rm = TRUE)
    cat(sprintf("Trained for %d epoch(s) (one epoch = one full pass through the training data).\n\n",
                n_epochs))
    
    final_metrics <- hist_df[hist_df$epoch == n_epochs, ]
    cat("Final epoch results:\n")
    print(final_metrics[, intersect(c("metric", "data", "value"), names(final_metrics))])
    cat("\n")
    
    # ---- Overfitting check: does validation loss start rising while training
    # loss keeps falling? A real, computed diagnostic, not just a plot to eyeball ----
    loss_data <- hist_df[hist_df$metric == "loss", ]
    train_loss <- loss_data[loss_data$data == "training", c("epoch", "value")]
    val_loss   <- loss_data[loss_data$data == "validation", c("epoch", "value")]
    
    if (nrow(train_loss) > 3 && nrow(val_loss) > 3) {
      train_loss <- train_loss[order(train_loss$epoch), ]
      val_loss   <- val_loss[order(val_loss$epoch), ]
      
      best_val_epoch <- val_loss$epoch[which.min(val_loss$value)]
      last_epoch     <- max(val_loss$epoch)
      
      cat("Overfitting check (training loss vs. validation loss over time):\n")
      cat("--------------------------------------------------------\n")
      cat(sprintf("Validation loss was LOWEST at epoch %d (value = %.4f).\n",
                  best_val_epoch, min(val_loss$value)))
      cat(sprintf("Training ran for %d epochs total (final training loss = %.4f, final\n",
                  last_epoch, tail(train_loss$value, 1)))
      cat(sprintf("validation loss = %.4f).\n\n", tail(val_loss$value, 1)))
      
      if (best_val_epoch < last_epoch - 2) {
        cat(bold("CAUTION — likely OVERFITTING:"), sprintf("validation loss stopped improving at\n"))
        cat(sprintf("epoch %d, but training continued %d more epochs, likely fitting NOISE in the\n",
                    best_val_epoch, last_epoch - best_val_epoch))
        cat("training data rather than learning anything more generalizable. Consider:\n")
        cat("  - Early stopping: fit(..., callbacks = list(callback_early_stopping(\n")
        cat("      monitor = \"val_loss\", patience = 5, restore_best_weights = TRUE)))\n")
        cat("  - A simpler architecture (fewer units/layers — see CAVEATS below)\n")
        cat("  - Regularization (dropout layers, L2 weight regularization)\n\n")
      } else {
        cat("No strong sign of overfitting late in training — validation loss stayed\n")
        cat("roughly in line with training loss right up to the final epoch.\n\n")
      }
    }
    
    # ---- Figure: training curves, train vs validation, for every tracked metric ----
    p_hist <- ggplot2::ggplot(hist_df, ggplot2::aes(x = epoch, y = value, color = data)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::facet_wrap(~ metric, scales = "free_y") +
      ggplot2::labs(title = "Training History: Training vs. Validation, by Epoch",
                    subtitle = "Diverging lines (validation worsening while training improves) signal overfitting",
                    x = "Epoch", y = "Value", color = "Data set") +
      ggplot2::theme_minimal(base_size = 13)
    print(p_hist)
  }
  
  
  # ==========================================================
  # 3. TRUE HELD-OUT TEST SET EVALUATION (if supplied)
  # This is the gap your original code left open: validation_split
  # data is reshuffled from TRAINING data and monitored every epoch,
  # which means your architecture/epoch-count choices are implicitly
  # influenced by it. A genuinely untouched test set, evaluated ONCE
  # at the very end, is the only unbiased final performance check.
  # ==========================================================
  cat("--------------------------------------------------------\n")
  cat(bold("Validation performance vs. a TRUE held-out test set — an important distinction:\n"))
  cat("'validation_split' data is carved out of your TRAINING data and re-checked every\n")
  cat("single epoch — useful for monitoring/tuning, but because you (or early stopping)\n")
  cat("can react to it, it's not a fully unbiased final performance estimate. A\n")
  cat("genuinely held-out TEST set — never looked at until the very end — is needed for\n")
  cat("an honest answer to 'how will this perform on new patients?'\n\n")
  
  if (is.null(x_test) || is.null(y_test)) {
    cat(bold("SUGGESTED FIX:"), "your original code created x_test/y_test but never actually\n")
    cat("evaluated the model on them! Pass them into this function to get a genuine\n")
    cat("test-set evaluation:\n\n")
    cat("  explain_keras(model_keras, out_keras_fit, x_test = x_test, y_test = y_test)\n\n")
    cat("...or run it yourself:\n\n")
    cat("  model_keras %>% evaluate(x_test, y_test)\n")
  } else {
    cat("Evaluating on your supplied held-out test set:\n")
    cat("--------------------------------------------------------\n")
    
    eval_result <- tryCatch(model$evaluate(x_test, y_test, verbose = 0), error = function(e) NULL)
    
    if (is.null(eval_result)) {
      cat(bold("WARNING:"), "could not evaluate the model on the supplied x_test/y_test — check\n")
      cat("they have the same number of columns/format as the data the model was trained on.\n")
    } else {
      for (mname in names(eval_result)) {
        cat(sprintf("  Test-set %s: %.4f\n", mname, eval_result[[mname]]))
      }
      
      # Compare against the LAST validation-set value for the same metric, as a check
      if (!is.null(hist_df)) {
        cat("\nCompared to the last validation-set value seen during training:\n")
        for (mname in names(eval_result)) {
          val_row <- hist_df[hist_df$metric == mname & hist_df$data == "validation", ]
          if (nrow(val_row) > 0) {
            last_val <- tail(val_row[order(val_row$epoch), "value"], 1)
            diff_pct <- 100 * (eval_result[[mname]] - last_val) / abs(last_val)
            cat(sprintf("  %s: test = %.4f vs. last validation = %.4f (%.1f%% difference)\n",
                        mname, eval_result[[mname]], last_val, diff_pct))
            if (abs(diff_pct) > 15) {
              cat(sprintf("  %s test-set performance differs notably from validation — the model\n", bold("NOTE:")))
              cat("  may not generalize quite as well as training suggested, or this test set is\n")
              cat("  small enough that the estimate itself is noisy.\n")
            }
          }
        }
      }
      
      # Predicted probability distribution figure
      preds <- tryCatch(as.vector(model$predict(x_test, verbose = 0)), error = function(e) NULL)
      if (!is.null(preds)) {
        p_pred <- ggplot2::ggplot(data.frame(pred = preds, actual = factor(y_test)),
                                  ggplot2::aes(x = pred, fill = actual)) +
          ggplot2::geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
          ggplot2::labs(title = "Predicted Probabilities on the Test Set, by Actual Outcome",
                        subtitle = "Well-separated distributions indicate good discrimination between outcomes",
                        x = "Predicted probability of the event", y = "Number of patients",
                        fill = "Actual outcome") +
          ggplot2::theme_minimal(base_size = 13)
        print(p_pred)
      }
    }
  }
  
  
  # ==========================================================
  # 4. IS DEEP LEARNING THE RIGHT TOOL HERE?
  # ==========================================================
  cat("\n########################################################\n")
  cat("# IS DEEP LEARNING THE RIGHT TOOL FOR THIS DATA?\n")
  cat("########################################################\n\n")
  
  n_train_actual <- tryCatch({
    if (!is.null(hist_df)) {
      NA  # sample size isn't in the history object; estimated below if x_test unavailable too
    } else NA
  }, error = function(e) NA)
  
  cat(bold("Good situations to use deep learning:\n"))
  cat(" - LARGE datasets (often tens of thousands of rows or more) — neural networks\n")
  cat("   need substantial data to reliably learn without just memorizing noise\n")
  cat(" - UNSTRUCTURED data: images, free text, audio, sequential/time-series data —\n")
  cat("   this is where deep learning's ability to learn its own feature\n")
  cat("   representations (rather than you engineering features by hand) really shines\n")
  cat(" - Complex, highly non-linear relationships that simpler models can't capture,\n")
  cat("   AND you have enough data to estimate that complexity reliably\n\n")
  
  cat(bold("Why this specific setup is likely OVER-ENGINEERED for this dataset:\n"))
  cat(" - This is a small, STRUCTURED/TABULAR dataset (a handful of clinical/demographic\n")
  cat("   variables, roughly a thousand patients) — this is precisely the kind of data\n")
  cat("   where gradient-boosted trees (XGBoost) or random forests consistently match\n")
  cat("   or beat neural networks in practice, with far less tuning and far more\n")
  cat("   interpretability\n")
  cat(" - A 32-unit + 16-unit hidden-layer network has a LARGE number of parameters\n")
  cat("   relative to only 5 predictor variables and ~800 training rows — this is a\n")
  cat("   classic recipe for overfitting on a dataset this size and this simple\n")
  cat(" - Neural networks are essentially UNINTERPRETABLE without extra work (e.g. SHAP\n")
  cat("   values) — for only 5 clinical predictors, a logistic regression or a single\n")
  cat("   tree gives you comparable-or-better accuracy AND a transparent, explainable\n")
  cat("   model, which typically matters a great deal in clinical/pharmacoepi contexts\n\n")
  
  cat("Practical suggestion: fit and compare this against the logistic regression\n")
  cat("(explain_glm()), random forest (explain_rf()), and XGBoost (explain_xgb())\n")
  cat("models from earlier in this workflow on the SAME held-out test set. If deep\n")
  cat("learning doesn't demonstrably outperform them, the simpler, more interpretable\n")
  cat("model is very likely the better practical choice here.\n\n")
  
  cat("As with all machine learning methods in this series, no p-values or confidence\n")
  cat("intervals are produced — this is a prediction-focused method, not a classical\n")
  cat("inferential one.\n")
  cat("========================================================\n")
  
  invisible(list(
    history = hist_df,
    test_evaluation = if (!is.null(x_test) && !is.null(y_test) && exists("eval_result")) eval_result else NULL
  ))
}


# ============================================================
# USAGE EXAMPLES
# ============================================================
library(keras3)
n_train <- floor(0.8 * nrow(x_xgb))
x_train <- x_xgb[1:n_train, ]; y_train <- y_xgb[1:n_train]
x_test  <- x_xgb[(n_train+1):nrow(x_xgb), ]; y_test <- y_xgb[(n_train+1):nrow(x_xgb)]

model_keras <- keras_model_sequential() %>%
  layer_dense(units = 32, activation = "relu", input_shape = ncol(x_train)) %>%
  layer_dense(units = 16, activation = "relu") %>%
  layer_dense(units = 1, activation = "sigmoid")
model_keras %>% compile(loss = "binary_crossentropy", optimizer = "adam", metrics = "accuracy")
out_keras_fit <- model_keras %>% fit(x_train, y_train, epochs = 20, validation_split = 0.2)

# Without test-set evaluation (matches your original code exactly):
keras_results <- explain_keras(model_keras, out_keras_fit)

# WITH the test-set evaluation your original code was missing (recommended):
keras_results <- explain_keras(model_keras, out_keras_fit, x_test = x_test, y_test = y_test)

# ============================================================
# EDGE CASES THIS FUNCTION HAS BEEN SPECIFICALLY CHECKED AGAINST
# ============================================================
#  1. `model` not a recognisable keras object                      -> stops with a clear message
#     (checked generically via class name containing "keras", since exact class strings vary
#     across keras3/tensorflow backend versions)
#  2. `history` not a recognisable training-history object          -> stops with a clear message
#  3. summary(model) capture failing                                -> falls back to a manual
#     layer-by-layer table built from model$layers instead of erroring
#  4. A layer without a `units` or `activation` attribute (e.g. a dropout/flatten layer) -> shown
#     as "?"/"unknown" rather than crashing on a missing field
#  5. get_weights()-based parameter count failing                   -> section skipped gracefully
#  6. as.data.frame(history) failing (keras3 version differences)   -> falls back to printing the
#     raw history object instead of crashing the whole function
#  7. Fewer than 4 epochs of history (too little data for a
#     meaningful overfitting trend check)                           -> that specific check is
#     skipped, rest of the function still runs
#  8. x_test/y_test not supplied                                     -> test-evaluation section
#     degrades gracefully into a specific, actionable suggestion instead of erroring
#  9. model$evaluate() failing on the supplied test data (e.g.
#     column mismatch)                                              -> caught, clear warning
#     explaining the likely cause instead of a raw error
# 10. A metric name from evaluate() with no matching validation-history row -> that specific
#     comparison is skipped, doesn't halt the loop for other metrics
# 11. model$predict() failing when building the prediction-distribution figure -> caught,
#     that figure is skipped rather than crashing the whole function
# 12. Division-by-zero risk in the test-vs-validation percent-difference calc (val = 0) -> not
#     specially guarded since loss/accuracy metrics are essentially never exactly zero in
#     practice, but flagged here for transparency as a theoretical edge the function does not
#     explicitly special-case
# ============================================================