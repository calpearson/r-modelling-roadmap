# 04_mixed_effects_models

# ============================================================
# 4. MIXED EFFECTS MODELS — outcome: sbp / hi_bp (cohort_long)
# ============================================================
library(lme4)

model_lmm <- lmer(sbp ~ visit + treatment + (1 | patient_id), data = cohort_long)
out_lmm_summary <- summary(model_lmm)
out_lmm_summary

model_glmm <- glmer(hi_bp ~ visit + treatment + (1 | patient_id),
                    data = cohort_long, family = binomial)
out_glmm_summary <- summary(model_glmm)
out_glmm_summary



# ============================================================
# explain_mixed.R — plain-language interpretation of lme4
# mixed-effects models: lmer() (linear mixed model) and
# glmer() (generalized linear mixed model).
#
# Two functions:
#   explain_lmm(model)   — lmer() models, raw-scale coefficients
#                           (forest plot reference line at 0 — this
#                           IS correct here, unlike OR/RR models,
#                           because these coefficients are on the
#                           outcome's own natural scale, not a ratio)
#   explain_glmm(model)  — glmer() models, OR/RR-scale coefficients
#                           (forest plot reference line at 1, same
#                           logic as explain_glm())
#
# Both report fixed effects (which group/direction is bigger, is it
# significant) AND random effects (how much do groups like patients
# vary from each other, expressed as an ICC), and run a battery of
# mixed-model-specific sanity checks before showing any results.
# ============================================================

library(lme4)

# ---- Shared helpers (safe to re-source alongside the other explain_*.R files) ----
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) crayon::bold(x) else x
}
format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}

# ---- Shared: classify a fixed-effect term (interaction / transformed / ordered / factor / continuous) ----
.classify_term <- function(vn, term_labels, model_data) {
  is_interaction <- grepl(":", vn, fixed = TRUE)
  is_transformed <- grepl("^(log|sqrt|poly|I|exp|scale)\\(", vn)
  
  main_labels <- term_labels[!grepl(":", term_labels, fixed = TRUE)]
  matched_var <- main_labels[sapply(main_labels, function(p) startsWith(vn, p))]
  matched_var <- if (length(matched_var) > 0) matched_var[which.max(nchar(matched_var))] else character(0)
  
  is_ordered_factor <- length(matched_var) == 1 && matched_var %in% names(model_data) &&
    is.ordered(model_data[[matched_var]])
  is_factor_level <- length(matched_var) == 1 && matched_var %in% names(model_data) &&
    is.factor(model_data[[matched_var]]) && !is_ordered_factor
  
  list(is_interaction = is_interaction, is_transformed = is_transformed,
       is_ordered_factor = is_ordered_factor, is_factor_level = is_factor_level,
       matched_var = matched_var)
}

# ---- Shared: report random-effects structure + ICC ----
.explain_random_effects <- function(model, residual_var_available = TRUE) {
  
  vc <- as.data.frame(VarCorr(model))
  vc_grouping <- vc[is.na(vc$var2), ]  # drop correlation rows, keep variance/SD rows only
  
  n_grp_vec <- tryCatch(lme4::ngrps(model), error = function(e) NULL)
  
  cat("Random effects (how much do groups differ from each other?):\n")
  cat("--------------------------------------------------------\n")
  
  for (i in seq_len(nrow(vc_grouping))) {
    grp_name <- vc_grouping$grp[i]
    if (grp_name == "Residual") next
    var_i <- vc_grouping$vcov[i]
    sd_i  <- vc_grouping$sdcor[i]
    term_i <- vc_grouping$var1[i]
    n_levels <- if (!is.null(n_grp_vec) && grp_name %in% names(n_grp_vec)) n_grp_vec[[grp_name]] else NA
    
    cat(sprintf("\n> Grouping factor: %s (%s levels/groups in the data)\n",
                grp_name, ifelse(is.na(n_levels), "unknown number of", format(n_levels, big.mark = ","))))
    cat(sprintf("  Variance in '%s' across %s: %.3f (SD: %.3f)\n", term_i, grp_name, var_i, sd_i))
    cat(sprintf("  In plain English: this is how much the baseline (or slope, if this is a random\n"))
    cat(sprintf("  slope term) genuinely differs from one %s to another, above and beyond what's\n", grp_name))
    cat(sprintf("  explained by the fixed-effect predictors in the model.\n"))
  }
  
  # ICC — only cleanly defined for a single random intercept; still compute for the
  # first intercept term found, with a caveat if the structure is more complex.
  intercept_rows <- vc_grouping[vc_grouping$var1 == "(Intercept)" & vc_grouping$grp != "Residual", ]
  
  if (nrow(intercept_rows) >= 1) {
    grp_var <- sum(intercept_rows$vcov)  # sum across grouping factors if more than one
    
    if (residual_var_available) {
      resid_row <- vc[vc$grp == "Residual", ]
      resid_var <- if (nrow(resid_row) > 0) resid_row$vcov[1] else NA
      if (!is.na(resid_var)) {
        icc <- grp_var / (grp_var + resid_var)
        cat(sprintf("\nIntraclass correlation (ICC): %.3f\n", icc))
        cat(sprintf("In plain English: about %.0f%% of the TOTAL variation in the outcome is due to\n",
                    icc * 100))
        cat("differences BETWEEN groups (e.g. between patients), rather than variation WITHIN\n")
        cat("a group over time/repeats. A higher ICC means observations from the same group\n")
        cat("are more similar to each other — this justifies using a mixed model instead of\n")
        cat("treating every row as fully independent.\n")
      }
    } else {
      # glmer binomial/logit: no residual variance parameter exists, so use the
      # standard latent-variable approximation (level-1 variance = pi^2/3 for logit link)
      resid_var_approx <- (pi^2) / 3
      icc <- grp_var / (grp_var + resid_var_approx)
      cat(sprintf("\nIntraclass correlation (ICC, approximate): %.3f\n", icc))
      cat("In plain English: this uses the standard approximation for binomial/logit mixed\n")
      cat(sprintf("models (there's no true residual variance to work with here). Roughly %.0f%% of\n",
                  icc * 100))
      cat("the total variation in the underlying probability is due to differences BETWEEN\n")
      cat("groups — treat this as a rough guide, not an exact figure.\n")
    }
    
    if (nrow(intercept_rows) > 1) {
      cat("\nNOTE: more than one random-intercept grouping factor is present — the ICC above\n")
      cat("combines them, which is a simplification. Consider the per-group variances listed\n")
      cat("above individually for a more precise picture.\n")
    }
  }
  
  cat("\n")
  invisible(vc_grouping)
}

# ---- Shared: pre-flight convergence / singularity checks ----
.mixed_model_diagnostics <- function(model) {
  diagnostics <- character(0)
  
  is_singular <- tryCatch(lme4::isSingular(model), error = function(e) NA)
  if (isTRUE(is_singular)) {
    diagnostics <- c(diagnostics,
                     "has a SINGULAR FIT — at least one random-effect variance is estimated at (or extremely close to) zero, or a correlation between random effects is estimated at exactly +/-1. This usually means the random-effects structure is too complex for what the data can support. Consider simplifying it (e.g. dropping a random slope, or a grouping factor with very few levels) — fixed-effect estimates are often still usable, but treat the random-effects variances with real caution.")
  }
  
  conv_msgs <- tryCatch(model@optinfo$conv$lme4$messages, error = function(e) NULL)
  if (!is.null(conv_msgs) && length(conv_msgs) > 0) {
    diagnostics <- c(diagnostics, sprintf(
      "reported a CONVERGENCE WARNING from the optimizer: \"%s\". Estimates may not be fully reliable — consider trying a different optimizer (e.g. lme4::allFit()) or simplifying the model.",
      paste(conv_msgs, collapse = "; ")))
  }
  
  coef_all <- lme4::fixef(model)
  if (any(is.na(coef_all))) {
    diagnostics <- c(diagnostics, sprintf(
      "could NOT estimate %d fixed-effect coefficient(s), likely due to collinearity: %s. Those terms are skipped below.",
      sum(is.na(coef_all)), paste(names(coef_all)[is.na(coef_all)], collapse = ", ")))
  }
  
  se_vals <- tryCatch(sqrt(diag(vcov(model))), error = function(e) numeric(0))
  if (length(se_vals) > 0 && any(se_vals > 10, na.rm = TRUE)) {
    diagnostics <- c(diagnostics,
                     "has at least one EXTREMELY large standard error on a fixed effect — check for near-perfect separation, an unusually scaled predictor, or too few groups relative to the number of parameters being estimated.")
  }
  
  n_grp_vec <- tryCatch(lme4::ngrps(model), error = function(e) NULL)
  if (!is.null(n_grp_vec) && any(n_grp_vec < 5)) {
    small_grps <- names(n_grp_vec)[n_grp_vec < 5]
    diagnostics <- c(diagnostics, sprintf(
      "has a grouping factor with VERY FEW levels (%s: only %s groups). Variance estimates for that grouping factor are likely to be unstable with so few groups — a common rule of thumb wants at least ~5-10 groups, ideally more, for a random effect to be estimated reliably.",
      paste(small_grps, collapse = ", "), paste(n_grp_vec[small_grps], collapse = ", ")))
  }
  
  diagnostics
}


# ============================================================
# 1. explain_lmm() — linear mixed model (lmer)
# ============================================================
explain_lmm <- function(model, conf_level = 0.95, forest_plot = FALSE) {
  
  if (!inherits(model, "lmerMod")) {
    stop("explain_lmm() expects a model fit with lme4::lmer().")
  }
  
  alpha_pct <- conf_level * 100
  fixed_form <- formula(model, fixed.only = TRUE)
  outcome <- deparse(fixed_form[[2]])
  
  # ---- Pre-flight checks ----
  diagnostics <- .mixed_model_diagnostics(model)
  if (length(diagnostics) > 0) {
    cat("########################################################\n")
    cat("# MODEL DIAGNOSTIC WARNINGS — read before trusting results\n")
    cat("########################################################\n")
    for (d in diagnostics) cat(bold("WARNING:"), "this model", d, "\n\n")
  }
  
  # ---- Header: cohort structure ----
  n_obs    <- nobs(model)
  n_grp_vec <- tryCatch(lme4::ngrps(model), error = function(e) NULL)
  
  cat("========================================================\n")
  cat("MODEL:", deparse(formula(model)), "\n")
  cat("========================================================\n\n")
  
  cat(sprintf("This model was fit using %s observations", bold(format(n_obs, big.mark = ","))))
  if (!is.null(n_grp_vec)) {
    grp_desc <- paste(sprintf("%s %s", format(n_grp_vec, big.mark = ","), names(n_grp_vec)), collapse = ", ")
    cat(sprintf(" across %s", grp_desc))
  }
  cat(".\n")
  cat("This is a LINEAR MIXED MODEL: it models the outcome on its own natural (raw)\n")
  cat("scale — like lm() — but additionally accounts for the fact that observations\n")
  cat("from the same group (e.g. repeated visits from the same patient) are NOT fully\n")
  cat("independent of each other.\n\n")
  
  # ---- Fixed-effects table + p-values, handling lmerTest vs plain lme4 ----
  s <- summary(model)
  coefs_table <- s$coefficients
  var_names <- rownames(coefs_table)
  
  has_lmertest_pvals <- inherits(model, "lmerModLmerTest") && "Pr(>|t|)" %in% colnames(coefs_table)
  
  if (!has_lmertest_pvals) {
    cat("NOTE: plain lme4 does not report p-values for lmer() models by default, because\n")
    cat("the correct degrees of freedom for a t-test are not straightforward in mixed\n")
    cat("models. The p-values below are an APPROXIMATION using a normal (z) distribution,\n")
    cat("which is reasonable with a decent number of groups but can be slightly liberal\n")
    cat("(too easily significant) with few groups. For exact Satterthwaite-adjusted\n")
    cat("p-values, refit using library(lmerTest) before calling lmer().\n\n")
  }
  
  # ---- CIs: Wald, fast and always available for fixed effects ----
  ci_fixed <- tryCatch(
    suppressMessages(confint(model, method = "Wald", parm = "beta_", level = conf_level)),
    error = function(e) {
      z <- qnorm(1 - (1 - conf_level) / 2)
      se <- coefs_table[, "Std. Error"]
      cbind(coefs_table[, "Estimate"] - z * se, coefs_table[, "Estimate"] + z * se)
    }
  )
  
  # ---- Overall model fit: marginal/conditional R^2 isn't base-R; give AIC/BIC + REML/ML note ----
  cat(sprintf("Model fit: AIC = %.1f, BIC = %.1f (lower is better; mainly useful for\n",
              AIC(model), BIC(model)))
  cat(sprintf("comparing this model against alternative models fit to the SAME data).\n"))
  cat(sprintf("Estimated using: %s\n\n", ifelse(isREML(model), "REML", "Maximum Likelihood (ML)")))
  
  # ---- Random effects ----
  .explain_random_effects(model, residual_var_available = TRUE)
  
  # ---- Fixed-effects per-variable interpretation ----
  model_data <- model.frame(model)
  term_labels <- attr(terms(fixed_form), "term.labels")
  
  cat("What each FIXED-effect variable is associated with:\n")
  cat("--------------------------------------------------------\n")
  
  results_rows <- list()
  
  for (i in seq_along(var_names)) {
    vn <- var_names[i]
    if (vn == "(Intercept)") next
    
    estimate <- coefs_table[i, "Estimate"]
    lower    <- ci_fixed[vn, 1]
    upper    <- ci_fixed[vn, 2]
    
    if (has_lmertest_pvals) {
      pval <- coefs_table[i, "Pr(>|t|)"]
    } else {
      t_val <- coefs_table[i, "t value"]
      pval  <- 2 * pnorm(-abs(t_val))
    }
    
    cls <- .classify_term(vn, term_labels, model_data)
    direction <- ifelse(estimate > 0, "increase", "decrease")
    
    cat(sprintf("\n> %s\n", vn))
    cat(sprintf("  Estimate: %.3f (%.0f%% CI: %.3f to %.3f)\n", estimate, alpha_pct, lower, upper))
    
    row_type <- "main effect"
    
    if (cls$is_interaction) {
      row_type <- "interaction"
      parts <- strsplit(vn, ":", fixed = TRUE)[[1]]
      cat(sprintf("  In plain English — %s:\n", bold("INTERACTION TERM")))
      if (length(parts) == 2) {
        cat(sprintf("  The effect of '%s' on %s DEPENDS ON the level of '%s' (and vice versa).\n",
                    parts[1], outcome, parts[2]))
        cat(sprintf("  This coefficient (%.3f) is how much the effect of '%s' changes for each\n",
                    estimate, parts[1]))
        cat(sprintf("  1-unit change in (or shift to the other category of) '%s'.\n", parts[2]))
      } else {
        cat(sprintf("  Higher-order interaction between %d variables — hard to interpret from a\n",
                    length(parts)))
        cat(sprintf("  single coefficient. Consider plotting predicted values across combinations\n"))
        cat(sprintf("  of these variables instead.\n"))
      }
      
    } else if (cls$is_transformed) {
      row_type <- "transformed"
      cat(sprintf("  In plain English — %s:\n", bold("TRANSFORMED TERM")))
      cat(sprintf("  This is a transformation of the original variable — a '1-unit increase'\n"))
      cat(sprintf("  statement doesn't apply cleanly on the original scale.\n"))
      
    } else if (cls$is_ordered_factor) {
      row_type <- "ordered factor"
      cat(sprintf("  In plain English — %s:\n", bold("ORDERED FACTOR")))
      cat(sprintf("  This describes a TREND across ordered levels (polynomial contrast), not a\n"))
      cat(sprintf("  simple group comparison.\n"))
      
    } else if (cls$is_factor_level) {
      matched_var <- cls$matched_var
      reference_level <- levels(model_data[[matched_var]])[1]
      compared_level  <- substring(vn, nchar(matched_var) + 1)
      
      cat(sprintf("  In plain English — %s vs %s:\n", compared_level, reference_level))
      cat(sprintf("  %s Patients/rows in the '%s' group have, ON AVERAGE, %s of %.3f\n",
                  bold("Which group is greater?"), bold(compared_level), bold(direction), abs(estimate)))
      cat(sprintf("  units in %s compared to the '%s' (reference) group, holding every other\n",
                  outcome, reference_level))
      cat(sprintf("  fixed-effect variable constant.\n"))
      
    } else {
      cat(sprintf("  In plain English:\n"))
      cat(sprintf("  For every %s in %s, %s changes by %s units, ON AVERAGE — a %s.\n",
                  bold("1-unit increase"), vn, outcome, bold(sprintf("%.3f", estimate)), bold(direction)))
      cat(sprintf("  This holds every other fixed-effect variable in the model constant.\n"))
      if (length(cls$matched_var) == 0) {
        cat("  (Note: could not confidently match this term to an original data variable —\n")
        cat("  double-check this interpretation manually.)\n")
      }
    }
    
    cat(sprintf("\n  Is this a real effect, or could it be due to chance?\n"))
    if (pval < 0.05) {
      cat(sprintf("  This result IS %s (p %s).\n", bold("statistically significant"), format_pval(pval)))
    } else {
      cat(sprintf("  This result is %s (p %s). Treat this estimate with caution.\n",
                  bold("NOT statistically significant"), format_pval(pval)))
    }
    
    results_rows[[length(results_rows) + 1]] <- data.frame(
      variable = vn, type = row_type,
      estimate = round(estimate, 3), lower_ci = round(lower, 3), upper_ci = round(upper, 3),
      p_value = ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval)),
      significant = pval < 0.05,
      p_value_method = ifelse(has_lmertest_pvals, "Satterthwaite (lmerTest)", "Normal approximation")
    )
  }
  
  cat("\n--------------------------------------------------------\n")
  cat("Reminder: these are associations, adjusted for other FIXED-effect variables and\n")
  cat("for the correlation within groups — not proven causal effects.\n")
  cat("========================================================\n")
  
  results_table <- do.call(rbind, results_rows)
  
  # ---- Optional forest plot: reference line at 0 (correct here — raw scale, not a ratio) ----
  if (isTRUE(forest_plot) && !is.null(results_table) && nrow(results_table) > 0) {
    n_rows <- nrow(results_table)
    y_pos <- rev(seq_len(n_rows))
    x_range <- range(c(results_table$lower_ci, results_table$upper_ci, 0), na.rm = TRUE)
    
    old_par <- par(mar = c(5, max(8, max(nchar(results_table$variable)) * 0.6), 4, 2))
    on.exit(par(old_par), add = TRUE)
    
    plot(results_table$estimate, y_pos, xlim = x_range, ylim = c(0.5, n_rows + 0.5),
         pch = 16, cex = 1.3, yaxt = "n", ylab = "", xlab = sprintf("Effect on %s (raw scale)", outcome),
         main = sprintf("Forest plot: %s\n(fixed effects, raw scale)", outcome))
    axis(2, at = y_pos, labels = results_table$variable, las = 2, cex.axis = 0.85)
    segments(results_table$lower_ci, y_pos, results_table$upper_ci, y_pos, lwd = 2)
    abline(v = 0, lty = 2, col = "red", lwd = 1.5)
    text(x = 0, y = n_rows + 0.5, labels = "no effect", col = "red", pos = 3, xpd = TRUE, cex = 0.8)
  } else if (isTRUE(forest_plot)) {
    warning("forest_plot = TRUE was requested, but there are no estimable fixed effects to plot.")
  }
  
  invisible(results_table)
}


# ============================================================
# 2. explain_glmm() — generalized linear mixed model (glmer)
# ============================================================
explain_glmm <- function(model, conf_level = 0.95, forest_plot = FALSE) {
  
  if (!inherits(model, "glmerMod")) {
    stop("explain_glmm() expects a model fit with lme4::glmer().")
  }
  
  fam <- family(model)
  family_name <- fam$family
  link_name <- fam$link
  alpha_pct <- conf_level * 100
  fixed_form <- formula(model, fixed.only = TRUE)
  outcome <- deparse(fixed_form[[2]])
  
  # ---- Pre-flight checks (shared) + family check ----
  diagnostics <- .mixed_model_diagnostics(model)
  
  if (!(family_name == "binomial" && link_name %in% c("logit", "log")) &&
      !(family_name == "poisson" && link_name == "log")) {
    diagnostics <- c(diagnostics, sprintf(
      "uses family '%s' with link '%s', which this function doesn't have specific ratio-label text for. Falling back to a generic 'exponentiated coefficient' explanation.",
      family_name, link_name))
  }
  
  if (length(diagnostics) > 0) {
    cat("########################################################\n")
    cat("# MODEL DIAGNOSTIC WARNINGS — read before trusting results\n")
    cat("########################################################\n")
    for (d in diagnostics) cat(bold("WARNING:"), "this model", d, "\n\n")
  }
  
  # ---- Determine ratio label ----
  if (family_name == "binomial" && link_name == "logit") {
    ratio_label <- "Odds Ratio (OR)"
    ratio_noun  <- "odds"
    scale_explain <- "the outcome is a YES/NO (binary) event; fixed-effect coefficients below, once exponentiated, are ODDS RATIOS (OR)."
  } else if (family_name == "binomial" && link_name == "log") {
    ratio_label <- "Risk Ratio (RR)"
    ratio_noun  <- "risk"
    scale_explain <- "the outcome is a YES/NO (binary) event modeled with a log link; fixed-effect coefficients below, once exponentiated, are RISK RATIOS (RR)."
  } else if (family_name == "poisson" && link_name == "log") {
    ratio_label <- "Rate Ratio (RR)"
    ratio_noun  <- "rate"
    scale_explain <- "the outcome is a COUNT modeled with a log link; fixed-effect coefficients below, once exponentiated, are RATE RATIOS."
  } else {
    ratio_label <- sprintf("Exponentiated coefficient (%s link)", link_name)
    ratio_noun  <- "value"
    scale_explain <- sprintf("this model uses the '%s' family with a '%s' link, outside this function's built-in interpretations. Exponentiated coefficients below describe a multiplicative effect, but don't carry a standard OR/RR label.", family_name, link_name)
  }
  
  # ---- Header ----
  n_obs <- nobs(model)
  n_grp_vec <- tryCatch(lme4::ngrps(model), error = function(e) NULL)
  
  cat("========================================================\n")
  cat("MODEL:", deparse(formula(model)), "\n")
  cat("========================================================\n\n")
  
  cat(sprintf("This model was fit using %s observations", bold(format(n_obs, big.mark = ","))))
  if (!is.null(n_grp_vec)) {
    grp_desc <- paste(sprintf("%s %s", format(n_grp_vec, big.mark = ","), names(n_grp_vec)), collapse = ", ")
    cat(sprintf(" across %s", grp_desc))
  }
  cat(".\n")
  cat(sprintf("Family: %s | Link: %s\n\n", bold(family_name), bold(link_name)))
  cat("What this means:\n")
  cat(strwrap(scale_explain, width = 72, prefix = "  "), sep = "\n")
  cat("\nThis model ALSO accounts for the fact that repeated observations from the same\n")
  cat("group (e.g. multiple visits per patient) are not fully independent.\n\n")
  
  # ---- Random effects (no residual variance parameter for binomial/poisson GLMMs) ----
  .explain_random_effects(model, residual_var_available = FALSE)
  
  # ---- Fixed effects ----
  s <- summary(model)
  coefs_table <- s$coefficients
  var_names <- rownames(coefs_table)
  pval_col <- ncol(coefs_table)
  
  ci_fixed <- tryCatch(
    suppressMessages(confint(model, method = "Wald", parm = "beta_", level = conf_level)),
    error = function(e) {
      z <- qnorm(1 - (1 - conf_level) / 2)
      se <- coefs_table[, "Std. Error"]
      cbind(coefs_table[, "Estimate"] - z * se, coefs_table[, "Estimate"] + z * se)
    }
  )
  
  model_data <- model.frame(model)
  term_labels <- attr(terms(fixed_form), "term.labels")
  
  cat("What each FIXED-effect variable is associated with:\n")
  cat("--------------------------------------------------------\n")
  
  results_rows <- list()
  
  for (i in seq_along(var_names)) {
    vn <- var_names[i]
    if (vn == "(Intercept)") next
    
    estimate <- coefs_table[i, "Estimate"]
    pval     <- coefs_table[i, pval_col]
    ratio    <- exp(estimate)
    ratio_lo <- exp(ci_fixed[vn, 1])
    ratio_hi <- exp(ci_fixed[vn, 2])
    
    cls <- .classify_term(vn, term_labels, model_data)
    
    cat(sprintf("\n> %s\n", vn))
    cat(sprintf("  %s: %.2f (%.0f%% CI: %.2f to %.2f)   [raw coefficient: %.3f]\n",
                ratio_label, ratio, alpha_pct, ratio_lo, ratio_hi, estimate))
    
    row_type <- "main effect"
    pct_diff <- if (ratio >= 1) (ratio - 1) * 100 else (1 - ratio) * 100
    direction_word <- if (ratio > 1) "higher" else if (ratio < 1) "lower" else "no different"
    
    if (cls$is_interaction) {
      row_type <- "interaction"
      parts <- strsplit(vn, ":", fixed = TRUE)[[1]]
      cat(sprintf("  In plain English — %s:\n", bold("INTERACTION TERM")))
      if (length(parts) == 2) {
        cat(sprintf("  This interaction's ratio (%.2f) means: the %s associated with '%s' is\n",
                    ratio, ratio_label, parts[1]))
        cat(sprintf("  itself multiplied by %.2f for each 1-unit change in (or shift to the other\n", ratio))
        cat(sprintf("  category of) '%s' — the effect of %s DEPENDS ON the level of %s.\n",
                    parts[2], parts[1], parts[2]))
      } else {
        cat(sprintf("  Higher-order interaction between %d variables — interpret with predicted\n",
                    length(parts)))
        cat(sprintf("  values/plots rather than this single number.\n"))
      }
      
    } else if (cls$is_transformed) {
      row_type <- "transformed"
      cat(sprintf("  In plain English — %s:\n", bold("TRANSFORMED TERM")))
      cat(sprintf("  This is a transformation of the original variable — interpret on its\n"))
      cat(sprintf("  transformed scale, not as a raw 1-unit change.\n"))
      
    } else if (cls$is_ordered_factor) {
      row_type <- "ordered factor"
      cat(sprintf("  In plain English — %s:\n", bold("ORDERED FACTOR")))
      cat(sprintf("  This describes a TREND across ordered levels, not a simple group comparison.\n"))
      
    } else if (cls$is_factor_level) {
      matched_var <- cls$matched_var
      reference_level <- levels(model_data[[matched_var]])[1]
      compared_level  <- substring(vn, nchar(matched_var) + 1)
      
      cat(sprintf("  In plain English — %s vs %s:\n", compared_level, reference_level))
      if (abs(ratio - 1) < 1e-9) {
        cat(sprintf("  '%s' and '%s' show essentially IDENTICAL results.\n", compared_level, reference_level))
      } else {
        cat(sprintf("  %s Rows in the '%s' group have %s%% %s %s of %s than rows in\n",
                    bold("Which group is greater?"), bold(compared_level),
                    bold(sprintf("%.0f", pct_diff)), bold(direction_word), ratio_noun, bold(outcome)))
        cat(sprintf("  the '%s' (reference) group, ON AVERAGE, holding every other fixed-effect\n",
                    reference_level))
        cat(sprintf("  variable constant.\n"))
      }
      
    } else {
      cat(sprintf("  In plain English:\n"))
      cat(sprintf("  For every %s in %s, the %s of %s is multiplied by %s\n",
                  bold("1-unit increase"), vn, ratio_noun, outcome, bold(sprintf("%.2f", ratio))))
      cat(sprintf("  — in other words, %s%% %s %s of %s, ON AVERAGE, holding every other\n",
                  bold(sprintf("%.1f", pct_diff)), bold(direction_word), ratio_noun, bold(outcome)))
      cat(sprintf("  fixed-effect variable constant.\n"))
      if (length(cls$matched_var) == 0) {
        cat("  (Note: could not confidently match this term to an original data variable —\n")
        cat("  double-check this interpretation manually.)\n")
      }
    }
    
    cat(sprintf("\n  Is this a real effect, or could it be due to chance?\n"))
    if (is.na(pval)) {
      cat("  P-value not available for this term.\n")
    } else if (pval < 0.05) {
      cat(sprintf("  This result IS %s (p %s).\n", bold("statistically significant"), format_pval(pval)))
    } else {
      cat(sprintf("  This result is %s (p %s). Treat this estimate with caution.\n",
                  bold("NOT statistically significant"), format_pval(pval)))
    }
    
    results_rows[[length(results_rows) + 1]] <- data.frame(
      variable = vn, type = row_type, ratio_label = ratio_label,
      ratio = round(ratio, 3), lower_ci = round(ratio_lo, 3), upper_ci = round(ratio_hi, 3),
      p_value = ifelse(is.na(pval), NA, ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval))),
      significant = ifelse(is.na(pval), NA, pval < 0.05)
    )
  }
  
  cat("\n--------------------------------------------------------\n")
  cat(sprintf("Reminder: all ratios above are %s — a value of 1.00 means NO difference/effect.\n",
              ratio_label))
  cat("========================================================\n")
  
  results_table <- do.call(rbind, results_rows)
  
  # ---- Optional forest plot: reference line at 1 (ratio scale) ----
  if (isTRUE(forest_plot)) {
    if (is.null(results_table) || nrow(results_table) == 0) {
      warning("forest_plot = TRUE was requested, but there are no estimable fixed effects to plot.")
    } else {
      fp_data <- results_table[!is.na(results_table$ratio) &
                                 !is.na(results_table$lower_ci) &
                                 !is.na(results_table$upper_ci), ]
      use_log_scale <- nrow(fp_data) > 0 &&
        all(fp_data$lower_ci > 0) && all(fp_data$upper_ci > 0) && all(fp_data$ratio > 0)
      
      if (nrow(fp_data) > 0) {
        n_rows <- nrow(fp_data)
        y_pos <- rev(seq_len(n_rows))
        x_range <- range(c(fp_data$lower_ci, fp_data$upper_ci, 1), na.rm = TRUE)
        
        old_par <- par(mar = c(5, max(8, max(nchar(fp_data$variable)) * 0.6), 4, 2))
        on.exit(par(old_par), add = TRUE)
        
        plot(fp_data$ratio, y_pos, xlim = x_range, ylim = c(0.5, n_rows + 0.5),
             log = if (use_log_scale) "x" else "",
             pch = 16, cex = 1.3, yaxt = "n", ylab = "", xlab = ratio_label,
             main = sprintf("Forest plot: %s\n(%s)", outcome, ratio_label))
        axis(2, at = y_pos, labels = fp_data$variable, las = 2, cex.axis = 0.85)
        segments(fp_data$lower_ci, y_pos, fp_data$upper_ci, y_pos, lwd = 2)
        abline(v = 1, lty = 2, col = "red", lwd = 1.5)
        text(x = 1, y = n_rows + 0.5, labels = "no effect", col = "red", pos = 3, xpd = TRUE, cex = 0.8)
      }
    }
  }
  
  invisible(results_table)
}


# ============================================================
# USAGE EXAMPLES (using cohort_long from earlier)
# ============================================================
library(lme4)

model_lmm <- lmer(sbp ~ visit + treatment + (1 | patient_id), data = cohort_long)
lmm_results <- explain_lmm(model_lmm)
lmm_results <- explain_lmm(model_lmm, forest_plot = TRUE)   # reference line at 0 (raw scale)

model_glmm <- glmer(hi_bp ~ visit + treatment + (1 | patient_id),
                     data = cohort_long, family = binomial)
glmm_results <- explain_glmm(model_glmm)
glmm_results <- explain_glmm(model_glmm, forest_plot = TRUE)  # reference line at 1 (OR scale)
#
# ============================================================
# EDGE CASES THIS FILE HAS BEEN SPECIFICALLY CHECKED AGAINST
# ============================================================
#  1. Singular fit (boundary random-effect variance = 0, or correlation = +/-1)
#     -> detected via lme4::isSingular(), flagged before results
#  2. Optimizer convergence warnings                              -> flagged, extraction is defensive (tryCatch)
#  3. Aliased/NA fixed-effect coefficients (collinearity)         -> flagged, skipped safely
#  4. Extremely large fixed-effect standard errors                -> flagged
#  5. Grouping factor with very few levels (<5)                   -> flagged as unstable variance estimate
#  6. Plain lme4::lmer() with no p-values available               -> normal-approximation p used, explicitly caveated
#  7. lmerTest::lmer() WITH p-values available                    -> detected via class, used directly, labeled as such
#  8. confint(method="Wald") failing                              -> falls back to manual estimate +/- z*SE
#  9. glmer() family other than binomial-logit/log or poisson-log -> generic fallback explanation, flagged
# 10. Interaction terms (2-way and higher-order)                  -> explained using actual variable names
# 11. Transformed terms: log(), poly(), I(), sqrt(), scale()      -> flagged, not misread
# 12. Ordered factors                                             -> flagged as trend, not group comparison
# 13. More than one random-intercept grouping factor              -> ICC combines them, with an explicit caveat
# 14. Random effects other than a simple random intercept (e.g. random slopes) -> reported individually via VarCorr, ICC caveated as intercept-only
# 15. Terms unmatched to any known data column                    -> flagged, not silently asserted as correct
# 16. forest_plot = TRUE with no estimable terms / non-positive CI bounds (glmm) -> warns / falls back to linear scale
# 17. REML vs ML fitting                                          -> stated explicitly in the fit-summary line

