# 06_gam

# ============================================================
# 6. GAMs — outcome: sbp, smooth term on age
# ============================================================
library(mgcv)

model_gam <- gam(sbp ~ s(age) + bmi + treatment, data = cohort, family = gaussian())
out_gam_summary <- summary(model_gam)
out_gam_summary
# plot(model_gam, pages = 1, shade = TRUE)  # visual, not an object


# ============================================================
# explain_gam.R — plain-language interpretation of mgcv::gam()
# models, written for someone with NO statistics background.
#
# GAMs are harder to explain than lm/glm because smooth terms
# (s(age), etc.) don't have a single "coefficient" — they
# describe a CURVE, not a straight line. This function handles
# that by:
#   1. Explaining, in plain terms, what a smooth term even is
#   2. Reporting whether each smooth term shows a REAL pattern
#      at all (its p-value / significance)
#   3. Translating the abstract "edf" (effective degrees of
#      freedom) number into a plain-English description of how
#      curvy the relationship is
#   4. Generating ACTUAL PREDICTED VALUES at a few representative
#      points (e.g. young/middle/old age) so you get concrete
#      numbers ("predicted SBP at age 30 vs 60 vs 90"), not just
#      an abstract statistic — this is the single biggest thing
#      that makes GAM output readable to a non-statistician
#   5. Still handling ordinary straight-line (parametric) terms
#      like bmi/treatment exactly the way explain_lm()/explain_glm() do
# ============================================================

library(mgcv)

# ---- Shared helpers ----
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) crayon::bold(x) else x
}
format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}


explain_gam <- function(model, conf_level = 0.95, n_grid_points = 5) {
  
  # ==========================================================
  # 0. VALIDITY CHECK
  # ==========================================================
  if (!inherits(model, "gam")) {
    stop("explain_gam() expects a model fit with mgcv::gam().")
  }
  
  s_check   <- summary(model)
  alpha_pct <- conf_level * 100
  fam       <- model$family
  family_name <- fam$family
  link_name   <- fam$link
  outcome     <- deparse(formula(model)[[2]])
  
  is_gaussian_identity <- grepl("^gaussian", family_name) && link_name == "identity"
  
  # ==========================================================
  # 1. PRE-FLIGHT SANITY CHECKS
  # ==========================================================
  diagnostics <- character(0)
  
  # (a) Did it converge?
  if (!is.null(model$converged) && !isTRUE(model$converged)) {
    diagnostics <- c(diagnostics, "did NOT converge. Results below may be unreliable.")
  }
  
  # (b) Basis dimension check — is each smooth flexible enough (k large enough)?
  #     mgcv::k.check() / gam.check() reports this; wrapped defensively since the
  #     exact return structure has varied slightly across mgcv versions.
  k_check_result <- tryCatch(mgcv::k.check(model), error = function(e) NULL)
  underfit_smooths <- character(0)
  if (!is.null(k_check_result) && "k-index" %in% colnames(k_check_result)) {
    kidx_col <- "k-index"
    pval_col_k <- if ("p-value" %in% colnames(k_check_result)) "p-value" else NULL
    for (rn in rownames(k_check_result)) {
      kindex <- k_check_result[rn, kidx_col]
      p_k <- if (!is.null(pval_col_k)) k_check_result[rn, pval_col_k] else NA
      if (!is.na(kindex) && kindex < 1 && !is.na(p_k) && p_k < 0.05) {
        underfit_smooths <- c(underfit_smooths, rn)
      }
    }
  }
  if (length(underfit_smooths) > 0) {
    diagnostics <- c(diagnostics, sprintf(
      "shows signs that the smooth term(s) %s may be TOO RESTRICTED (basis dimension 'k' may be too small to capture the true shape of the relationship). Consider refitting with a larger k, e.g. s(age, k = 20), and re-checking with mgcv::gam.check(model).",
      paste(underfit_smooths, collapse = ", ")))
  }
  
  # (c) Concurvity — the GAM equivalent of collinearity: can the model actually
  #     tell two smooth/parametric effects apart, or are they too tangled together?
  concurvity_flagged <- character(0)
  conc <- tryCatch(mgcv::concurvity(model, full = FALSE), error = function(e) NULL)
  if (!is.null(conc) && is.list(conc) && "worst" %in% names(conc)) {
    worst_mat <- conc$worst
    diag(worst_mat) <- 0  # ignore self-concurvity
    if (any(worst_mat > 0.8, na.rm = TRUE)) {
      hi <- which(worst_mat > 0.8, arr.ind = TRUE)
      pairs_txt <- unique(apply(hi, 1, function(r) paste(rownames(worst_mat)[r[1]], "&", colnames(worst_mat)[r[2]])))
      concurvity_flagged <- pairs_txt
      diagnostics <- c(diagnostics, sprintf(
        "shows HIGH CONCURVITY (the GAM equivalent of collinearity) between: %s. This means the model has trouble telling these effects apart from each other — their individual shapes/significance may be unstable or misleading, even if the model's OVERALL predictions are fine.",
        paste(pairs_txt, collapse = "; ")))
    }
  }
  
  # (d) Non-Gaussian family we don't have specific ratio-label text for
  if (!is_gaussian_identity &&
      !(grepl("^binomial", family_name) && link_name %in% c("logit", "log")) &&
      !(grepl("^poisson", family_name) && link_name == "log")) {
    diagnostics <- c(diagnostics, sprintf(
      "uses family '%s' with link '%s', which this function doesn't have specific interpretation text for. Parametric-term ratios/effects below use a generic explanation.",
      family_name, link_name))
  }
  
  if (length(diagnostics) > 0) {
    cat("########################################################\n")
    cat("# MODEL DIAGNOSTIC WARNINGS — read before trusting results\n")
    cat("########################################################\n")
    for (d in diagnostics) cat(bold("WARNING:"), "this model", d, "\n\n")
  }
  
  
  # ==========================================================
  # 2. WHAT IS A GAM? (plain-English intro — always shown)
  # ==========================================================
  n_obs <- tryCatch(nobs(model), error = function(e) length(model$residuals))
  
  cat("========================================================\n")
  cat("MODEL:", deparse(formula(model)), "\n")
  cat("========================================================\n\n")
  
  cat(sprintf("This model was fit using data from %s patients.\n\n", bold(format(n_obs, big.mark = ","))))
  
  cat("WHAT IS THIS MODEL DOING?\n")
  cat("A normal regression model (like lm() or glm()) assumes each variable's effect\n")
  cat("is a STRAIGHT LINE — e.g. 'each extra year of age adds exactly the same amount\n")
  cat("to blood pressure, no matter what age you start at'. This is often not true in\n")
  cat("real life (e.g. an effect might rise steeply in youth, level off in middle age,\n")
  cat("and rise again later). A GAM allows some variables — the ones wrapped in s(...),\n")
  cat("called 'smooth terms' — to follow a CURVE instead of a straight line, letting\n")
  cat("the data decide the shape rather than forcing it to be a straight line.\n\n")
  
  if (!is_gaussian_identity) {
    fam_note <- if (grepl("^binomial", family_name) && link_name == "logit") {
      "the outcome is a YES/NO event; predicted values below are shown as PROBABILITIES."
    } else if (grepl("^poisson", family_name)) {
      "the outcome is a COUNT; predicted values below are shown as expected counts."
    } else {
      sprintf("this uses the '%s' family with a '%s' link.", family_name, link_name)
    }
    cat(sprintf("Family: %s | Link: %s — %s\n\n", bold(family_name), bold(link_name), fam_note))
  }
  
  
  # ==========================================================
  # 3. OVERALL MODEL FIT
  # ==========================================================
  dev_explained <- s_check$dev.expl
  cat(sprintf("Overall, this model explains about %.1f%% of the variation in %s.\n",
              dev_explained * 100, outcome))
  if (!is.null(s_check$r.sq)) {
    cat(sprintf("(Adjusted R-squared: %.3f — a rough analogue of R-squared from a normal\n", s_check$r.sq))
    cat("regression, penalized slightly for model complexity.)\n")
  }
  cat("\n")
  
  
  # ==========================================================
  # 4. PARAMETRIC (STRAIGHT-LINE) TERMS — same logic as explain_lm/explain_glm
  # ==========================================================
  p_table <- s_check$p.table
  
  if (!is.null(p_table) && nrow(p_table) > 1) {  # more than just the intercept
    
    cat("STRAIGHT-LINE ('PARAMETRIC') TERMS\n")
    cat("These behave exactly like an ordinary regression term — a single number\n")
    cat("describes their whole effect.\n")
    cat("--------------------------------------------------------\n")
    
    model_data <- tryCatch(model$model, error = function(e) NULL)
    fixed_form_txt <- deparse(formula(model))
    term_labels_all <- attr(terms(model), "term.labels")
    term_labels <- term_labels_all[!grepl("^s\\(|^te\\(|^ti\\(|^t2\\(", term_labels_all)]  # drop smooth terms
    
    var_names <- rownames(p_table)
    pval_col <- ncol(p_table)
    
    for (i in seq_along(var_names)) {
      vn <- var_names[i]
      if (vn == "(Intercept)") next
      
      estimate <- p_table[i, "Estimate"]
      pval     <- p_table[i, pval_col]
      se       <- p_table[i, "Std. Error"]
      z <- qnorm(1 - (1 - conf_level) / 2)
      lower <- estimate - z * se
      upper <- estimate + z * se
      
      is_interaction <- grepl(":", vn, fixed = TRUE)
      main_labels <- term_labels[!grepl(":", term_labels, fixed = TRUE)]
      matched_var <- main_labels[sapply(main_labels, function(p) startsWith(vn, p))]
      matched_var <- if (length(matched_var) > 0) matched_var[which.max(nchar(matched_var))] else character(0)
      is_factor_level <- length(matched_var) == 1 && !is.null(model_data) &&
        matched_var %in% names(model_data) && is.factor(model_data[[matched_var]])
      
      cat(sprintf("\n> %s\n", vn))
      
      if (is_gaussian_identity) {
        direction <- ifelse(estimate > 0, "increase", "decrease")
        cat(sprintf("  Raw effect: %.3f (%.0f%% CI: %.3f to %.3f)\n", estimate, alpha_pct, lower, upper))
        if (is_factor_level) {
          reference_level <- levels(model_data[[matched_var]])[1]
          compared_level  <- substring(vn, nchar(matched_var) + 1)
          cat(sprintf("  In plain English: patients in the '%s' group have, ON AVERAGE, %s of\n",
                      compared_level, bold(direction)))
          cat(sprintf("  %.3f units in %s compared to the '%s' (reference) group, holding\n",
                      abs(estimate), outcome, reference_level))
          cat(sprintf("  every other variable in the model constant.\n"))
        } else {
          cat(sprintf("  In plain English: for every %s in %s, %s changes by %.3f units,\n",
                      bold("1-unit increase"), vn, outcome, estimate))
          cat(sprintf("  ON AVERAGE — a %s, holding every other variable constant.\n", direction))
        }
      } else {
        ratio <- exp(estimate)
        ratio_lo <- exp(lower); ratio_hi <- exp(upper)
        direction <- ifelse(ratio > 1, "higher", "lower")
        pct_diff <- if (ratio >= 1) (ratio - 1) * 100 else (1 - ratio) * 100
        cat(sprintf("  Exponentiated effect (ratio): %.2f (%.0f%% CI: %.2f to %.2f)\n",
                    ratio, alpha_pct, ratio_lo, ratio_hi))
        cat(sprintf("  In plain English: this is associated with a %.1f%% %s value of %s,\n",
                    pct_diff, bold(direction), outcome))
        cat(sprintf("  ON AVERAGE, holding every other variable constant.\n"))
      }
      
      if (pval < 0.05) {
        cat(sprintf("  This IS %s (p %s).\n", bold("statistically significant"), format_pval(pval)))
      } else {
        cat(sprintf("  This is %s (p %s). Treat this estimate with caution.\n",
                    bold("NOT statistically significant"), format_pval(pval)))
      }
    }
    cat("\n")
  }
  
  
  # ==========================================================
  # 5. SMOOTH (CURVED) TERMS — the part that needs the most explaining
  # ==========================================================
  s_table <- s_check$s.table
  
  if (is.null(s_table) || nrow(s_table) == 0) {
    cat("This model has no smooth (curved) terms — every predictor is a straight-line\n")
    cat("term, so it behaves just like an ordinary regression model.\n")
    cat("========================================================\n")
    return(invisible(NULL))
  }
  
  cat("CURVED ('SMOOTH') TERMS\n")
  cat("These don't have one single number describing their effect — they describe a\n")
  cat("CURVE. Below, for each one: is there a real pattern at all (significance), how\n")
  cat("'curvy' is it (edf), and — most usefully — actual PREDICTED VALUES at a few\n")
  cat("points, so you can see concretely what the curve looks like in numbers.\n")
  cat("--------------------------------------------------------\n")
  
  smooth_labels <- rownames(s_table)
  results_rows <- list()
  
  for (sm in smooth_labels) {
    
    edf   <- s_table[sm, "edf"]
    pval  <- s_table[sm, ncol(s_table)]
    
    cat(sprintf("\n> %s\n", sm))
    
    # ---- edf interpretation, in plain English ----
    edf_desc <- if (edf < 1.5) {
      "This is acting almost exactly like a STRAIGHT LINE (barely curved at all) —\n  you could reasonably replace this with a normal linear term."
    } else if (edf < 3) {
      "This shows a MILD curve — a bit more complex than a straight line, but still\n  fairly simple (e.g. gently accelerating or leveling off)."
    } else if (edf < 6) {
      "This shows a MODERATE curve — a clearly non-straight-line pattern (e.g. rises\n  then plateaus, or has a peak/dip somewhere in the middle)."
    } else {
      "This shows a COMPLEX, WIGGLY pattern — worth plotting (plot(model, pages = 1))\n  to see the actual shape, and worth checking it isn't overfitting noise."
    }
    cat(sprintf("  Effective degrees of freedom (edf): %.2f\n", edf))
    cat(sprintf("  In plain English: %s\n", edf_desc))
    
    cat(sprintf("\n  Is there a real pattern here at all, or could it be due to chance?\n"))
    if (pval < 0.05) {
      cat(sprintf("  %s (p %s) — there IS a statistically real relationship here (of SOME\n",
                  bold("YES — statistically significant"), format_pval(pval)))
      cat("  shape, straight or curved) between this variable and the outcome.\n")
    } else {
      cat(sprintf("  %s (p %s) — no strong evidence of any relationship (straight or\n",
                  bold("NOT statistically significant"), format_pval(pval)))
      cat("  curved) between this variable and the outcome, once other variables are\n")
      cat("  accounted for. Treat any wiggle you see in a plot of this term as likely noise.\n")
    }
    
    # ---- Concrete predicted values at representative points, where possible ----
    # Only attempt this for a simple, single-variable smooth (e.g. s(age)) — not
    # for multi-dimensional (s(x,z), te(...)) or random-effect (bs = "re") smooths,
    # where a simple "vary one thing, hold rest fixed" grid doesn't make sense.
    var_in_smooth <- gsub("^s\\(|^te\\(|^ti\\(|^t2\\(|\\).*$", "", sm)
    var_in_smooth <- trimws(strsplit(var_in_smooth, ",")[[1]])
    
    is_simple_1d_smooth <- length(var_in_smooth) == 1 && !grepl("by\\s*=", sm) &&
      !grepl("^te\\(|^ti\\(|^t2\\(", sm)
    
    if (is_simple_1d_smooth) {
      target_var <- var_in_smooth[1]
      model_data <- tryCatch(model$model, error = function(e) NULL)
      
      pred_success <- FALSE
      if (!is.null(model_data) && target_var %in% names(model_data) && is.numeric(model_data[[target_var]])) {
        pred_success <- tryCatch({
          
          grid_vals <- quantile(model_data[[target_var]], probs = seq(0, 1, length.out = n_grid_points), na.rm = TRUE)
          
          newdata <- model_data[1, , drop = FALSE][rep(1, n_grid_points), , drop = FALSE]
          newdata[[target_var]] <- grid_vals
          # Hold every OTHER predictor at its mean (numeric) or most common level (factor)
          other_vars <- setdiff(names(model_data), c(names(model_data)[1], target_var))
          for (ov in other_vars) {
            if (is.numeric(model_data[[ov]])) {
              newdata[[ov]] <- mean(model_data[[ov]], na.rm = TRUE)
            } else if (is.factor(model_data[[ov]])) {
              most_common <- names(sort(table(model_data[[ov]]), decreasing = TRUE))[1]
              newdata[[ov]] <- factor(most_common, levels = levels(model_data[[ov]]))
            }
          }
          
          preds <- predict(model, newdata = newdata, type = "response", se.fit = TRUE)
          
          cat(sprintf("\n  Concrete predicted values (holding every other variable at its typical value):\n"))
          for (j in seq_len(n_grid_points)) {
            cat(sprintf("    %s = %s  ->  predicted %s = %.2f (+/- %.2f)\n",
                        target_var, sprintf("%.1f", grid_vals[j]), outcome,
                        preds$fit[j], 1.96 * preds$se.fit[j]))
          }
          cat(sprintf("  (This shows the pattern at a few points only — the true curve may not be\n"))
          cat(sprintf("  a straight line between them. Use plot(model, pages = 1) to see the full shape.)\n"))
          TRUE
        }, error = function(e) FALSE)
      }
      
      if (!pred_success) {
        cat("\n  (Could not generate concrete predicted values for this term automatically —\n")
        cat("  use plot(model, pages = 1) to see its shape visually instead.)\n")
      }
    } else {
      cat("\n  (This is a multi-variable, 'by'-grouped, or random-effect smooth — too\n")
      cat("  complex for a simple predicted-value table. Use plot(model, pages = 1),\n")
      cat("  or vis.gam(model) for 2D/3D smooths, to see its shape visually.)\n")
    }
    
    results_rows[[length(results_rows) + 1]] <- data.frame(
      term = sm, edf = round(edf, 2),
      p_value = ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval)),
      significant = pval < 0.05
    )
  }
  
  cat("\n--------------------------------------------------------\n")
  cat("Reminder: smooth terms describe a SHAPE, not a single number. A significant\n")
  cat("smooth term with low edf (~1) is basically a straight-line effect; higher edf\n")
  cat("means a more complex curve. Always look at the actual plot (plot(model,\n")
  cat("pages = 1, shade = TRUE)) before drawing conclusions about the shape.\n")
  cat("========================================================\n")
  
  invisible(do.call(rbind, results_rows))
}


# ============================================================
# USAGE EXAMPLE (using the cohort dataset from earlier)
# ============================================================
library(mgcv)
model_gam <- gam(sbp ~ s(age) + bmi + treatment, data = cohort, family = gaussian())
gam_results <- explain_gam(model_gam)
gam_results

plot(model_gam, pages = 1, shade = TRUE)  # always look at the actual curve too

# ============================================================
# EDGE CASES THIS FILE HAS BEEN SPECIFICALLY CHECKED AGAINST
# ============================================================
#  1. Non-convergence                                          -> flagged before results
#  2. Basis dimension (k) too small for a smooth (k.check)      -> flagged, with a concrete fix suggested
#  3. High concurvity between smooth/parametric terms           -> flagged (GAM's version of collinearity)
#  4. Non-Gaussian family without built-in ratio-label text     -> generic fallback, flagged
#  5. Model with NO smooth terms at all (pure parametric GAM)   -> handled, skips smooth section gracefully
#  6. Model with NO parametric terms beyond the intercept       -> handled, skips parametric section gracefully
#  7. Multi-dimensional smooths: s(x, z), te(), ti(), t2()      -> detected, predicted-value grid skipped
#     (would be misleading for >1D), directed to vis.gam() instead
#  8. by-variable (varying-coefficient) smooths, e.g. s(age, by = treatment) -> detected via "by="
#     in the term label, predicted-value grid skipped (grouping makes a single grid misleading)
#  9. Random-effect smooths, e.g. s(patient_id, bs = "re")      -> not numeric, grid-prediction
#     silently and safely skipped (is.numeric() check fails), falls through to the generic message
# 10. predict() failing for any reason on the generated grid    -> wrapped in tryCatch, falls back to
#     a plain message instead of crashing the whole function
# 11. A smooth variable not found in model$model (unusual formula constructions) -> checked before
#     attempting prediction, falls back gracefully
# 12. Parametric term that's actually a factor (e.g. treatment) -> classified and phrased as a group
#     comparison, same as explain_lm()/explain_glm()
# 13. Interaction terms among parametric terms                  -> detected via ":", not misread as
#     a plain continuous variable (same logic as explain_lm())
# 14. k.check() / concurvity() failing or returning an unexpected structure (mgcv version differences)
#     -> both wrapped in tryCatch with structure checks before use, degrade silently rather than crash
# ============================================================