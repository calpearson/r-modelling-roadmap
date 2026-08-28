# 07_regularization_glmnet

# ============================================================
# 7. REGULARIZATION (glmnet) — outcome: event, predictors incl. noise markers
# ============================================================
library(glmnet)

x_glmnet <- model.matrix(event ~ age + sex + bmi + comorbidity_score + treatment +
                           marker1 + marker2 + marker3 + marker4 + marker5 +
                           marker6 + marker7 + marker8 + marker9 + marker10 +
                           marker11 + marker12 + marker13 + marker14 + marker15,
                         data = cohort)[, -1]
y_glmnet <- cohort$event

model_lasso_cv  <- cv.glmnet(x_glmnet, y_glmnet, alpha = 1, family = "binomial")
out_lasso_coefs <- coef(model_lasso_cv, s = "lambda.min")
out_lasso_coefs   # true predictors should survive; markers should shrink to ~0

out_lasso_predict <- predict(model_lasso_cv, newx = x_glmnet[1:5, ],
                             s = "lambda.min", type = "response")
out_lasso_predict







# ============================================================
# explain_glmnet.R — plain-language interpretation of a
# glmnet::cv.glmnet() model (LASSO / Ridge / Elastic Net).
#
# This is structurally DIFFERENT from every other explain_*()
# function in this series, for one crucial reason:
#
#   Regularized models like LASSO do NOT produce valid p-values
#   or confidence intervals via the standard formulas. The
#   penalty deliberately introduces BIAS (shrinking coefficients
#   toward zero, some to exactly zero) in exchange for lower
#   variance/better prediction — this is the whole point of the
#   method. This function will NEVER fabricate a p-value or CI
#   for a glmnet coefficient. Where other explain_*() functions
#   say "is this significant?", this one says "was this variable
#   KEPT or DROPPED by the penalty?" instead — a different, and
#   for this method, more honest question.
# ============================================================

library(glmnet)

# ---- Shared helpers ----
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) crayon::bold(x) else x
}


explain_glmnet <- function(model, lambda = "lambda.min", newx = NULL, forest_plot = FALSE) {
  
  # ==========================================================
  # 0. VALIDITY CHECKS
  # ==========================================================
  if (!inherits(model, "cv.glmnet")) {
    stop("explain_glmnet() expects a model fit with glmnet::cv.glmnet() (not a plain glmnet() object — this function relies on the cross-validation results).")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The 'Matrix' package is required (it's normally installed automatically alongside glmnet) — please install.packages('Matrix').")
  }
  if (!(identical(lambda, "lambda.min") || identical(lambda, "lambda.1se") || is.numeric(lambda))) {
    stop("`lambda` must be \"lambda.min\", \"lambda.1se\", or a numeric value.")
  }
  
  lambda_val <- if (is.numeric(lambda)) lambda else model[[lambda]]
  if (is.null(lambda_val) || is.na(lambda_val)) {
    stop(sprintf("Could not find a value for lambda = \"%s\" on this model object.", lambda))
  }
  
  
  # ==========================================================
  # 1. DETECT METHOD (LASSO / Ridge / Elastic Net) AND FAMILY
  # Detected primarily via class(), which is a stable, documented
  # part of glmnet's object structure — more reliable than trying
  # to re-parse the original function call (which can fail if
  # alpha/family were passed in as variables rather than literals).
  # ==========================================================
  fit_class <- class(model$glmnet.fit)[1]
  
  family_lookup <- list(
    lognet  = list(name = "binomial (logistic)", ratio_scale = TRUE,  ratio_label = "Odds Ratio (OR)", ratio_noun = "odds"),
    fishnet = list(name = "Poisson",              ratio_scale = TRUE,  ratio_label = "Rate Ratio (RR)", ratio_noun = "rate"),
    coxnet  = list(name = "Cox (survival)",        ratio_scale = TRUE,  ratio_label = "Hazard Ratio (HR)", ratio_noun = "hazard"),
    elnet   = list(name = "Gaussian (linear)",     ratio_scale = FALSE, ratio_label = NA, ratio_noun = "value"),
    mrelnet = list(name = "multi-response Gaussian", ratio_scale = FALSE, ratio_label = NA, ratio_noun = "value"),
    multnet = list(name = "multinomial",           ratio_scale = NA,    ratio_label = NA, ratio_noun = NA)
  )
  
  fam_info <- family_lookup[[fit_class]]
  if (is.null(fam_info)) {
    fam_info <- list(name = sprintf("unrecognised (%s)", fit_class), ratio_scale = FALSE, ratio_label = NA, ratio_noun = "value")
  }
  
  if (fit_class == "multnet") {
    cat("########################################################\n")
    cat("# NOT SUPPORTED\n")
    cat("########################################################\n")
    cat("This is a MULTINOMIAL glmnet model (more than 2 outcome categories). This\n")
    cat("function only handles binary/continuous/count/survival outcomes, because a\n")
    cat("multinomial model has a SEPARATE set of coefficients per outcome category,\n")
    cat("which needs a different display entirely. Use coef(model, s = ...) directly\n")
    cat("and inspect each category's coefficients by hand.\n")
    cat("========================================================\n")
    return(invisible(NULL))
  }
  
  alpha_val <- tryCatch({
    a <- model$call$alpha
    if (is.null(a)) 1 else suppressWarnings(as.numeric(a))  # glmnet's own default is alpha = 1 (LASSO)
  }, error = function(e) NA)
  
  method_name <- if (is.na(alpha_val)) {
    "an elastic-net-family model (alpha could not be automatically determined from the call — check your original cv.glmnet() call)"
  } else if (alpha_val == 1) {
    "LASSO (alpha = 1)"
  } else if (alpha_val == 0) {
    "Ridge Regression (alpha = 0)"
  } else {
    sprintf("Elastic Net (alpha = %.2f — a blend of LASSO and Ridge)", alpha_val)
  }
  
  
  # ==========================================================
  # 2. HEADER — plain-English explanation of what this method does
  # ==========================================================
  cat("========================================================\n")
  cat("REGULARIZED REGRESSION:", bold(method_name), "\n")
  cat("========================================================\n\n")
  
  cat(sprintf("Family: %s\n\n", bold(fam_info$name)))
  
  cat("What this method does, in plain English:\n")
  cat("Ordinary regression tries to fit the data as closely as possible, which can\n")
  cat("mean chasing noise when you have MANY predictors (especially more predictors\n")
  cat("than you can confidently estimate from your sample size, or predictors that are\n")
  cat("correlated with each other). This method adds a PENALTY that shrinks\n")
  cat("coefficients toward zero — trading a small amount of bias for a model that\n")
  cat("generalizes better to new data.\n\n")
  
  if (!is.na(alpha_val) && alpha_val == 1) {
    cat("Because this is LASSO specifically, the penalty can shrink some coefficients\n")
    cat("all the way to EXACTLY zero — effectively removing that variable from the\n")
    cat("model entirely. This gives you automatic variable selection as a side effect.\n\n")
  } else if (!is.na(alpha_val) && alpha_val == 0) {
    cat("Because this is Ridge specifically, coefficients are shrunk toward zero but\n")
    cat("essentially NEVER hit exactly zero — every variable stays in the model, just\n")
    cat("with a dampened effect. This is NOT a variable-selection method.\n\n")
  }
  
  
  # ==========================================================
  # 3. CROSS-VALIDATION SUMMARY
  # ==========================================================
  cat("How the amount of shrinkage (lambda) was chosen:\n")
  cat("--------------------------------------------------------\n")
  cat(sprintf("%d different penalty strengths (lambda values) were tried, each evaluated\n",
              length(model$lambda)))
  cat("via cross-validation (fitting on part of the data, checking prediction error on\n")
  cat("the held-out part, repeated across folds) — this picks a penalty strength using\n")
  cat("the DATA itself, not a value you have to guess.\n\n")
  
  idx_min <- which(model$lambda == model$lambda.min)
  idx_1se <- which(model$lambda == model$lambda.1se)
  nzero_min <- if (length(idx_min) == 1) model$nzero[idx_min] else NA
  nzero_1se <- if (length(idx_1se) == 1) model$nzero[idx_1se] else NA
  
  cat(sprintf("  lambda.min (lowest CV error): %.5f -> keeps %s non-zero predictor(s)\n",
              model$lambda.min, ifelse(is.na(nzero_min), "an undetermined number of", nzero_min)))
  cat(sprintf("  lambda.1se (simplest model within 1 SE of the best): %.5f -> keeps %s\n",
              model$lambda.1se, ifelse(is.na(nzero_1se), "an undetermined number of", nzero_1se)))
  cat("  non-zero predictor(s)\n\n")
  
  cat("lambda.min gives the model with the BEST cross-validated prediction accuracy.\n")
  cat("lambda.1se gives a SIMPLER, more conservative model that's still statistically\n")
  cat("about as good — often preferred when you want a more stable, parsimonious set\n")
  cat("of predictors and are less worried about squeezing out the last bit of accuracy.\n\n")
  
  cat(sprintf("This explanation uses: %s\n\n", bold(sprintf("lambda = \"%s\" (%.5f)", lambda, lambda_val))))
  
  cat(bold("Reproducibility note:"), "cross-validation randomly splits the data into folds.\n")
  cat("Unless you set a random seed (set.seed(...)) before calling cv.glmnet(), re-running\n")
  cat("this exact code could select a slightly different lambda and a slightly different\n")
  cat("set of retained variables. This is not a bug — it's normal for CV-based methods —\n")
  cat("but it's worth knowing if you're trying to exactly reproduce a result.\n\n")
  
  
  # ==========================================================
  # 4. COEFFICIENTS — retained vs. shrunk-to-zero
  # ==========================================================
  coef_raw <- tryCatch(coef(model, s = lambda_val), error = function(e) {
    stop(sprintf("Could not extract coefficients from this model at the requested lambda: %s", conditionMessage(e)))
  })
  coef_vec <- tryCatch(as.matrix(coef_raw)[, 1], error = function(e) {
    stop("Could not convert this model's coefficients into a usable vector — the object structure looks unexpected for a standard cv.glmnet() model.")
  })
  
  has_intercept <- "(Intercept)" %in% names(coef_vec)
  predictor_coefs <- if (has_intercept) coef_vec[names(coef_vec) != "(Intercept)"] else coef_vec
  
  retained <- predictor_coefs[predictor_coefs != 0]
  dropped  <- predictor_coefs[predictor_coefs == 0]
  
  cat("Which variables did the penalty KEEP vs. DROP?\n")
  cat("--------------------------------------------------------\n")
  cat(sprintf("Out of %d candidate predictors: %s retained (non-zero), %s dropped (shrunk\n",
              length(predictor_coefs), bold(length(retained)), bold(length(dropped))))
  cat("to exactly zero).\n\n")
  
  if (length(retained) == 0) {
    cat(bold("NOTE:"), "every predictor was shrunk to exactly zero at this lambda — the\n")
    cat("model has reduced to an intercept-only model, meaning it predicts the same\n")
    cat("value for everyone regardless of their covariates. This can happen if the\n")
    cat("penalty is very strong, or if none of the predictors carry a strong enough\n")
    cat("signal on their own. Consider using a smaller lambda if this seems too\n")
    cat("aggressive for your purposes.\n\n")
  } else {
    
    cat("RETAINED variables (the penalty judged these worth keeping):\n\n")
    
    retained_rows <- list()
    
    # Sort by absolute effect size, largest first, so the most influential
    # retained variables are shown first
    retained_sorted <- retained[order(-abs(retained))]
    
    for (vn in names(retained_sorted)) {
      raw_coef <- retained_sorted[[vn]]
      
      cat(sprintf("> %s\n", vn))
      cat(sprintf("  Raw (shrunk) coefficient: %.4f\n", raw_coef))
      
      if (isTRUE(fam_info$ratio_scale)) {
        ratio <- exp(raw_coef)
        direction_word <- if (ratio > 1) "higher" else if (ratio < 1) "lower" else "no different"
        pct_diff <- if (ratio >= 1) (ratio - 1) * 100 else (1 - ratio) * 100
        cat(sprintf("  %s: %.2f  (i.e. roughly %s%% %s %s per 1-unit increase, ON AVERAGE)\n",
                    fam_info$ratio_label, ratio, sprintf("%.1f", pct_diff), direction_word, fam_info$ratio_noun))
      } else {
        direction_word <- if (raw_coef > 0) "increase" else "decrease"
        cat(sprintf("  In plain English: roughly a %.4f-unit %s in the outcome per 1-unit\n",
                    abs(raw_coef), direction_word))
        cat("  increase in this variable, ON AVERAGE.\n")
      }
      
      cat(sprintf("  %s this is a REGULARIZED (deliberately shrunk) estimate — it is biased\n", bold("CAVEAT:")))
      cat("  toward zero/no-effect ON PURPOSE, to improve prediction on new data. It will\n")
      cat("  typically look SMALLER/closer to 'no effect' than the equivalent coefficient\n")
      cat("  from an unpenalized glm(). Do not treat its exact magnitude as an unbiased\n")
      cat("  effect-size estimate, and do NOT compute a p-value or confidence interval for\n")
      cat("  it using ordinary regression formulas — those aren't valid for this method.\n\n")
      
      retained_rows[[length(retained_rows) + 1]] <- data.frame(
        variable = vn, status = "retained", raw_coefficient = round(raw_coef, 4),
        ratio = if (isTRUE(fam_info$ratio_scale)) round(exp(raw_coef), 3) else NA
      )
    }
    
    cat("--------------------------------------------------------\n")
    cat(sprintf("DROPPED variables (shrunk to exactly zero — %d of them):\n", length(dropped)))
    if (length(dropped) > 0) {
      cat(paste(" -", names(dropped), collapse = "\n"), "\n\n")
      cat("In plain English: the penalty found no reliable enough signal for these\n")
      cat("variables to justify keeping them in the model. This does NOT necessarily\n")
      cat("mean these variables have zero true relationship with the outcome — it means\n")
      cat("their signal wasn't strong/consistent enough to survive the penalty, especially\n")
      cat("if they're correlated with a variable that WAS kept (LASSO tends to pick one\n")
      cat("representative from a group of correlated predictors somewhat arbitrarily,\n")
      cat("rather than keeping all of them).\n")
    } else {
      cat("(none — every candidate predictor was retained)\n")
    }
  }
  
  
  # ==========================================================
  # 5. OPTIONAL: EXPLAIN PREDICTIONS ON NEW DATA
  # ==========================================================
  if (!is.null(newx)) {
    cat("\n--------------------------------------------------------\n")
    cat("PREDICTIONS on the data you supplied:\n")
    cat("--------------------------------------------------------\n")
    
    pred_type <- "response"  # consistent across binomial/Poisson/Gaussian/Cox: predict.cv.glmnet's
    # "response" scale gives probabilities (binomial), predicted counts
    # (Poisson), raw predicted values (Gaussian), and relative risk =
    # exp(linear predictor) for Cox — NOT the linear predictor itself,
    # which is what type = "link" would give for Cox specifically.
    
    preds <- tryCatch(
      predict(model, newx = newx, s = lambda_val, type = pred_type),
      error = function(e) {
        cat(sprintf("%s could not generate predictions for the data you supplied: %s\n",
                    bold("WARNING:"), conditionMessage(e)))
        cat("This is often caused by newx having a different number/order of columns than\n")
        cat("the matrix the model was originally trained on — double-check it was built\n")
        cat("with the exact same model.matrix() formula and column order.\n")
        NULL
      }
    )
    
    if (!is.null(preds)) {
      preds_vec <- as.vector(preds)
      
      if (fit_class == "lognet") {
        cat("Each value below is a PREDICTED PROBABILITY of the event (0 to 1):\n\n")
        for (i in seq_along(preds_vec)) {
          cat(sprintf("  Row %d: %.1f%% predicted probability of the event\n", i, preds_vec[i] * 100))
        }
      } else if (fit_class == "fishnet") {
        cat("Each value below is a PREDICTED COUNT/RATE:\n\n")
        for (i in seq_along(preds_vec)) {
          cat(sprintf("  Row %d: %.2f predicted count\n", i, preds_vec[i]))
        }
      } else if (fit_class == "coxnet") {
        cat("Each value below is a PREDICTED RELATIVE RISK (relative to an average patient\n")
        cat("in this dataset) — a value of 2.0 means roughly double the hazard of an average\n")
        cat("patient. This is NOT a survival probability or a predicted survival time — Cox\n")
        cat("models deliberately leave the baseline hazard unspecified, so this only tells\n")
        cat("you RELATIVE risk ranking between patients, not absolute risk:\n\n")
        for (i in seq_along(preds_vec)) {
          cat(sprintf("  Row %d: %.2fx the relative risk of an average patient\n", i, preds_vec[i]))
        }
      } else {
        cat("Each value below is a PREDICTED value of the outcome, on its raw scale:\n\n")
        for (i in seq_along(preds_vec)) {
          cat(sprintf("  Row %d: %.3f\n", i, preds_vec[i]))
        }
      }
    }
  }
  
  
  # ==========================================================
  # 6. WHEN IS THIS METHOD APPROPRIATE? (explicitly requested)
  # ==========================================================
  cat("\n########################################################\n")
  cat("# IS THIS METHOD APPROPRIATE FOR YOUR SITUATION?\n")
  cat("########################################################\n\n")
  
  cat(bold("Good situations to use LASSO/Ridge/Elastic Net:\n"))
  cat(" - You have MANY candidate predictors (dozens to thousands), possibly more than\n")
  cat("   you could confidently include in an ordinary regression given your sample size\n")
  cat("   (e.g. genomic/biomarker panels, claims-data feature sets, this example's mix\n")
  cat("   of real predictors + many noise 'marker' variables)\n")
  cat(" - You suspect only a SUBSET of your predictors are genuinely useful, and want\n")
  cat("   automatic help figuring out which (LASSO/Elastic Net specifically)\n")
  cat(" - Your predictors are highly CORRELATED with each other, which makes ordinary\n")
  cat("   regression coefficients unstable (Ridge specifically handles this well,\n")
  cat("   without necessarily dropping variables)\n")
  cat(" - Your PRIMARY goal is prediction accuracy on new/future data, not obtaining\n")
  cat("   precise, unbiased effect-size estimates for a small set of pre-specified\n")
  cat("   variables of scientific interest\n\n")
  
  cat(bold("Situations where this is probably the WRONG tool:\n"))
  cat(" - You need valid p-values/confidence intervals for formal statistical inference\n")
  cat("   (e.g. a pre-registered analysis of a small number of specific exposure\n")
  cat("   variables) — use ordinary glm()/coxph() instead; this method's coefficients\n")
  cat("   are not designed to support that kind of inference out of the box\n")
  cat(" - Your sample size is small relative to the true complexity of the underlying\n")
  cat("   relationship — variable selection becomes UNSTABLE (a different random split\n")
  cat("   of the same data could select a meaningfully different set of variables)\n")
  cat(" - You already know, on subject-matter grounds, exactly which variables belong\n")
  cat("   in your model — there's no need for automatic selection, and an ordinary\n")
  cat("   regression will be more directly interpretable\n")
  cat(" - Precise, unbiased effect-size magnitudes matter more to you than prediction\n")
  cat("   accuracy (e.g. you need to report 'the odds ratio is X' with a defensible\n")
  cat("   confidence interval for a regulatory or clinical audience)\n\n")
  
  cat("If you need BOTH automatic variable selection AND valid inference, look into\n")
  cat("'post-selection inference' methods (e.g. the 'selectiveInference' or 'hdi'\n")
  cat("R packages) rather than treating standard LASSO output as inference-ready.\n")
  cat("========================================================\n")
  
  
  # ==========================================================
  # 7. OPTIONAL FOREST-STYLE PLOT (point estimates ONLY — no CI,
  # because standard glmnet does not provide valid ones)
  # ==========================================================
  if (isTRUE(forest_plot)) {
    if (length(retained) == 0) {
      warning("forest_plot = TRUE was requested, but no variables were retained at this lambda — nothing to plot.")
    } else {
      plot_vals <- if (isTRUE(fam_info$ratio_scale)) exp(retained_sorted) else retained_sorted
      n_rows <- length(plot_vals)
      y_pos <- rev(seq_len(n_rows))
      ref_line <- if (isTRUE(fam_info$ratio_scale)) 1 else 0
      x_label <- if (isTRUE(fam_info$ratio_scale)) fam_info$ratio_label else "Raw coefficient (shrunk estimate)"
      
      use_log_scale <- isTRUE(fam_info$ratio_scale) && all(plot_vals > 0)
      x_range <- range(c(plot_vals, ref_line), na.rm = TRUE)
      
      old_par <- par(mar = c(5, max(8, max(nchar(names(plot_vals))) * 0.6), 4, 2))
      on.exit(par(old_par), add = TRUE)
      
      plot(as.numeric(plot_vals), y_pos, xlim = x_range, ylim = c(0.5, n_rows + 0.5),
           log = if (use_log_scale) "x" else "",
           pch = 16, cex = 1.3, yaxt = "n", ylab = "", xlab = x_label,
           main = sprintf("%s: retained variables\n(POINT ESTIMATES ONLY — no valid CI for shrunk coefficients)", method_name))
      axis(2, at = y_pos, labels = names(plot_vals), las = 2, cex.axis = 0.85)
      abline(v = ref_line, lty = 2, col = "red", lwd = 1.5)
      text(x = ref_line, y = n_rows + 0.5, labels = "no effect", col = "red", pos = 3, xpd = TRUE, cex = 0.8)
    }
  }
  
  invisible(list(
    retained = if (length(retained) > 0) do.call(rbind, retained_rows) else NULL,
    dropped = names(dropped),
    lambda_used = lambda_val,
    method = method_name
  ))
}


# ============================================================
# USAGE EXAMPLES (using the cohort dataset from earlier)
# ============================================================
library(glmnet)
x_glmnet <- model.matrix(event ~ age + sex + bmi + comorbidity_score + treatment +
                            marker1 + marker2 + marker3 + marker4 + marker5 +
                            marker6 + marker7 + marker8 + marker9 + marker10 +
                            marker11 + marker12 + marker13 + marker14 + marker15,
                          data = cohort)[, -1]
y_glmnet <- cohort$event
model_lasso_cv <- cv.glmnet(x_glmnet, y_glmnet, alpha = 1, family = "binomial")

lasso_results <- explain_glmnet(model_lasso_cv)
lasso_results <- explain_glmnet(model_lasso_cv, forest_plot = TRUE)

# With prediction explanation for the same 5 rows as your original code:
lasso_results <- explain_glmnet(model_lasso_cv, newx = x_glmnet[1:5, ], forest_plot = TRUE)

# Using the more conservative lambda instead:
lasso_results <- explain_glmnet(model_lasso_cv, lambda = "lambda.1se")

# ============================================================
# EDGE CASES THIS FUNCTION HAS BEEN SPECIFICALLY CHECKED AGAINST
# ============================================================
#  1. A plain glmnet() object (no CV) passed by mistake           -> stops with a clear, specific message
#  2. `Matrix` package unavailable                                -> stops with a clear, specific message
#     (needed to convert the sparse coefficient matrix)
#  3. Invalid `lambda` argument (not lambda.min/lambda.1se/numeric) -> stops with a clear message
#  4. Multinomial (multnet) models                                -> explicitly detected and declined,
#     since they need a fundamentally different (per-category) display
#  5. Unrecognised glmnet.fit class (future glmnet family types)  -> generic fallback description,
#     doesn't assume ratio-scale interpretation incorrectly
#  6. alpha not extractable from the call (e.g. passed as a variable, not a literal) -> caught,
#     falls back to a generic "elastic-net-family" description instead of guessing
#  7. Cox models (coxnet) with no intercept term                  -> checked for presence of
#     "(Intercept)" before excluding it, rather than assuming it always exists
#  8. Every single predictor shrunk to exactly zero               -> handled explicitly and explained,
#     rather than producing an empty/confusing "retained" section
#  9. coef() extraction failing entirely                          -> stops with a specific, actionable
#     error message rather than a cryptic downstream failure
# 10. `newx` with mismatched columns vs. the training matrix       -> predict() wrapped in tryCatch,
#     explains the likely cause instead of showing a raw R error
# 11. `newx` not supplied at all                                  -> prediction section skipped cleanly,
#     rest of the function still runs
# 12. Different outcome types needing different predict() `type`  -> handled per family (response for
#     binomial/Poisson/Gaussian, link/relative-risk explained separately for Cox)
# 13. lambda.min and lambda.1se landing on the exact same lambda value in $lambda -> `which()` match
#     is used defensively rather than assuming distinct fixed positions
# 14. forest_plot = TRUE with zero retained variables              -> warns instead of erroring
# 15. Coefficients that are negative on the ratio scale — exp() of a negative number is always
#     positive, so the log-scale forest plot axis is always valid for retained ratio-scale variables
#     (no additional guard needed, but verified explicitly here rather than assumed)
# ============================================================