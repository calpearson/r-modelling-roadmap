# 02_glm_logistic_poisson

# ============================================================
# 2. GLMs — outcome: event (logistic), n_hosp (Poisson/NegBin)
# ============================================================
model_logit <- glm(event ~ age + sex + treatment + comorbidity_score,
                   data = cohort, family = binomial(link = "logit"))

out_logit_summary <- summary(model_logit)
out_logit_or      <- exp(coef(model_logit))       # odds ratios
out_logit_or_ci   <- exp(confint(model_logit))

out_logit_summary
out_logit_or
out_logit_or_ci

model_pois <- glm(n_hosp ~ age + treatment + offset(log(person_time)),
                  data = cohort, family = poisson(link = "log"))
out_pois_summary <- summary(model_pois)
out_pois_summary

library(MASS)
model_nb <- glm.nb(n_hosp ~ age + treatment + offset(log(person_time)), data = cohort)
out_nb_summary <- summary(model_nb)
out_nb_summary



# ============================================================
# explain_glm() — plain-language interpretation of a glm() /
# MASS::glm.nb() model. Handles logistic (binomial/logit),
# log-binomial (binomial/log), Poisson, and Negative Binomial
# models, with or without an offset, for ANY number/combination
# of predictors — continuous, factor, or interaction terms.
#
# WHAT THIS FUNCTION DOES DIFFERENTLY FROM explain_lm():
#  - States cohort size AND family/link up front, in plain English
#  - Exponentiates coefficients and labels them correctly:
#       binomial + logit  -> Odds Ratio (OR)
#       binomial + log    -> Risk Ratio (RR)
#       poisson/negbin + log + offset(log(...)) -> Incidence Rate
#                                                    Ratio (IRR)
#       poisson/negbin + log, no offset -> Rate Ratio (RR)
#       anything else     -> generic "multiplicative effect"
#  - States which group is bigger/smaller in plain English, using
#    the ACTUAL level names (not just "the reference group")
#  - Explains interaction terms as a "ratio of ratios" concept,
#    using the actual variable names in the interaction
#  - Runs a battery of pre-flight sanity checks before saying
#    anything, and reports on anything that looks like it could
#    make the numbers unreliable (see full list below)
#  - Optionally draws a forest plot of the ratios (OR/RR/IRR) with
#    95% CIs, on the exponentiated (ratio) scale. The reference
#    line for "no effect" is drawn at 1 — NOT 0 — because 1 is the
#    null value for ANY ratio (OR, RR, or IRR); 0 is only the null
#    value on the raw, un-exponentiated log-coefficient scale, which
#    isn't what a forest plot conventionally shows.
# ============================================================

# ---- Bold-text helper (safe if crayon isn't installed) ----
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    crayon::bold(x)
  } else {
    x
  }
}

# ---- p-value formatter ----
format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}


explain_glm <- function(model, conf_level = 0.95, forest_plot = FALSE) {
  
  # ==========================================================
  # 0. VALIDITY CHECK — is this even a model type we can handle?
  # ==========================================================
  if (!inherits(model, "glm")) {
    stop("explain_glm() expects a model fit with glm() or MASS::glm.nb().")
  }
  
  fam <- tryCatch(family(model), error = function(e) NULL)
  if (is.null(fam)) {
    stop("Could not determine this model's family/link — is it a valid glm object?")
  }
  
  is_negbin    <- inherits(model, "negbin")
  family_name  <- if (is_negbin) "Negative Binomial" else fam$family
  link_name    <- fam$link
  alpha_pct    <- conf_level * 100
  
  
  # ==========================================================
  # 1. PRE-FLIGHT SANITY CHECKS
  # Every one of these is a known way this kind of model — or
  # this function's own logic — could mislead you. Anything
  # triggered here is printed as a warning BEFORE any results,
  # so you see it before you see (and trust) a single number.
  # ==========================================================
  diagnostics <- character(0)
  
  # (a) Did the model actually converge?
  if (!is.null(model$converged) && !isTRUE(model$converged)) {
    diagnostics <- c(diagnostics,
                     "did NOT converge. Coefficients, SEs, and p-values below may be unreliable — refit with more iterations (glm(..., control = glm.control(maxit = 100))) or reconsider the model.")
  }
  
  # (b) Aliased / non-estimable coefficients (perfect collinearity)
  coef_all <- coef(model)
  if (any(is.na(coef_all))) {
    diagnostics <- c(diagnostics, sprintf(
      "could NOT estimate %d coefficient(s), almost always because of perfect collinearity between predictors (one variable is a exact linear combination of others): %s. Those terms are dropped from the interpretation below.",
      sum(is.na(coef_all)), paste(names(coef_all)[is.na(coef_all)], collapse = ", ")))
  }
  
  # (c) Suspiciously huge standard errors -> possible (quasi-)complete separation
  s_check <- summary(model)
  se_vals <- s_check$coefficients[, "Std. Error"]
  if (any(se_vals > 10, na.rm = TRUE)) {
    diagnostics <- c(diagnostics,
                     "has at least one EXTREMELY large standard error. This usually signals (quasi-)complete separation — e.g. a predictor perfectly (or near-perfectly) predicts the outcome within some subgroup, which is common with small samples or rare events. The odds/rate ratio and CI for the affected variable(s) may be absurdly wide or unstable — treat them with real caution.")
  }
  
  # (d) Unrecognised family — we'll still try, but flag it
  known_family_pattern <- "^(binomial|poisson|Negative Binomial|quasibinomial|quasipoisson|gaussian|Gamma|inverse\\.gaussian)"
  if (!grepl(known_family_pattern, family_name)) {
    diagnostics <- c(diagnostics, sprintf(
      "uses family '%s', which this function doesn't have specific interpretation text for. Falling back to a generic 'exponentiated coefficient' explanation — treat the ratio labels (OR/RR/IRR) with caution for this family.",
      family_name))
  }
  
  # (e) Response given as a 2-column matrix, e.g. cbind(successes, failures)
  resp <- model.response(model.frame(model))
  response_is_matrix <- is.matrix(resp)
  if (response_is_matrix) {
    diagnostics <- c(diagnostics,
                     "has a response given as cbind(successes, failures) — i.e. each ROW represents a GROUP of trials, not a single patient. Cohort size below reflects the number of GROUPS, not individual patients. Interpretation (OR/RR direction) is still valid.")
  }
  
  # (f) Non-uniform prior weights (aggregated / weighted data)
  has_weights <- !is.null(model$prior.weights) && length(unique(model$prior.weights)) > 1
  if (has_weights) {
    diagnostics <- c(diagnostics,
                     "was fit with non-uniform weights (e.g. aggregated data or survey weights). Coefficients are still valid, but 'cohort size' below reflects the number of ROWS in the data, not necessarily the number of independent patients.")
  }
  
  # (g) Extract offset (if any) directly from the formula text — robust to
  #     both offset(...) inside the formula and an offset= argument
  form_txt <- deparse(formula(model), width.cutoff = 500)
  form_txt <- paste(form_txt, collapse = " ")
  offset_match <- regmatches(form_txt, regexpr("offset\\([^)]*\\)", form_txt))
  offset_expr  <- if (length(offset_match) > 0 && nchar(offset_match) > 0) offset_match else NULL
  if (is.null(offset_expr) && !is.null(model$call$offset)) {
    offset_expr <- paste0("offset(", deparse(model$call$offset), ")")
  }
  has_offset   <- !is.null(offset_expr)
  offset_is_log <- has_offset && grepl("log\\(", offset_expr)
  if (has_offset && !offset_is_log) {
    diagnostics <- c(diagnostics, sprintf(
      "has an offset (%s) that does NOT look log-transformed. Offsets for count models should almost always be log(exposure_time) — an offset on the raw (non-logged) scale will produce coefficients that don't mean what this function assumes. Double-check this model.",
      offset_expr))
  }
  
  # (h) Confidence intervals: glm's default confint() uses profile likelihood,
  #     which can fail to converge on some models. Wrap it and fall back to
  #     a Wald-based interval (estimate +/- z*SE) if it errors, rather than
  #     letting the whole function crash.
  ci_link <- tryCatch(
    suppressMessages(confint(model, level = conf_level)),
    error = function(e) {
      diagnostics <<- c(diagnostics,
                        "profile-likelihood confidence intervals failed to compute (this can happen with small samples, separation, or unusual data). Falling back to approximate Wald CIs (estimate +/- z*SE), which are less accurate in small samples — treat interval widths with extra caution.")
      z <- qnorm(1 - (1 - conf_level) / 2)
      cbind(coef_all - z * s_check$coefficients[, "Std. Error"],
            coef_all + z * s_check$coefficients[, "Std. Error"])
    }
  )
  # confint() can return a plain vector (not a matrix) if there's only 1 non-intercept
  # coefficient — force it back into matrix form so later indexing doesn't break.
  if (is.null(dim(ci_link))) {
    ci_link <- matrix(ci_link, nrow = 1, dimnames = list(names(coef_all)[2], c("lo", "hi")))
  }
  
  # Print all triggered diagnostics now, before any results
  if (length(diagnostics) > 0) {
    cat("########################################################\n")
    cat("# MODEL DIAGNOSTIC WARNINGS — read before trusting results\n")
    cat("########################################################\n")
    for (d in diagnostics) cat(bold("WARNING:"), "this model", d, "\n\n")
  }
  
  
  # ==========================================================
  # 2. WORK OUT THE CORRECT RATIO LABEL (OR / RR / IRR / generic)
  # ==========================================================
  outcome <- deparse(formula(model)[[2]])
  
  if (grepl("^(binomial|quasibinomial)", family_name)) {
    if (link_name == "logit") {
      ratio_label <- "Odds Ratio (OR)"
      scale_explain <- paste0(
        "the outcome is a YES/NO (binary) event, and this model estimates the probability ",
        "of that event using a LOGIT link. Coefficients below, once exponentiated, are ",
        "ODDS RATIOS (OR): how many times higher or lower the ODDS of the event are.")
    } else if (link_name == "log") {
      ratio_label <- "Risk Ratio (RR)"
      scale_explain <- paste0(
        "the outcome is a YES/NO (binary) event, and this model estimates the probability ",
        "of that event using a LOG link. Coefficients below, once exponentiated, are ",
        "RISK RATIOS (RR): how many times higher or lower the PROBABILITY of the event is ",
        "(this is a more directly interpretable number than an odds ratio, but log-binomial ",
        "models can be numerically unstable — check the warnings above).")
    } else {
      ratio_label <- sprintf("Exponentiated coefficient (%s link)", link_name)
      scale_explain <- paste0(
        "the outcome is a YES/NO (binary) event, modeled with a '", link_name, "' link, which ",
        "isn't logit or log. Exponentiating the coefficient does not have a standard OR/RR ",
        "interpretation for this link — treat the numbers below as a rough multiplicative ",
        "guide only.")
    }
  } else if (grepl("^(poisson|quasipoisson)", family_name) || is_negbin) {
    if (link_name == "log") {
      if (has_offset && offset_is_log) {
        ratio_label <- "Incidence Rate Ratio (IRR)"
        scale_explain <- sprintf(paste0(
          "the outcome is a COUNT of events (e.g. hospitalisations), modeled as a RATE per ",
          "unit of follow-up time via the offset %s. Coefficients below, once exponentiated, ",
          "are INCIDENCE RATE RATIOS (IRR): how many times higher or lower the RATE of events ",
          "is, per unit of exposure/follow-up time."), offset_expr)
      } else {
        ratio_label <- "Rate Ratio (RR)"
        scale_explain <- paste0(
          "the outcome is a COUNT of events, modeled with a log link but WITHOUT an offset for ",
          "exposure/follow-up time. Coefficients below, once exponentiated, describe ratios in ",
          "the EXPECTED COUNT — this is only equivalent to a true 'rate' if every patient had ",
          "the same amount of follow-up time, which is worth checking.")
      }
    } else {
      ratio_label <- sprintf("Exponentiated coefficient (%s link)", link_name)
      scale_explain <- paste0(
        "the outcome is a count, modeled with a '", link_name, "' link, which isn't the usual ",
        "log link. Exponentiating does not have the standard IRR interpretation here.")
    }
  } else {
    ratio_label <- sprintf("Exponentiated coefficient (%s link)", link_name)
    scale_explain <- paste0(
      "this model uses the '", family_name, "' family with a '", link_name, "' link, which is ",
      "outside this function's built-in interpretations (typically used for continuous, ",
      "always-positive outcomes). The exponentiated coefficients below describe a ",
      "MULTIPLICATIVE effect on the mean of ", outcome, ", but don't carry a standard OR/RR/IRR ",
      "label.")
  }
  
  
  # ==========================================================
  # 3. HEADER — cohort size, family/link, what's being modeled
  # ==========================================================
  n_obs <- tryCatch(nobs(model), error = function(e) length(model$residuals))
  
  cat("========================================================\n")
  cat("MODEL:", deparse(formula(model)), "\n")
  cat("========================================================\n\n")
  
  cat(sprintf("This model was fit using data from %s %s.\n",
              bold(format(n_obs, big.mark = ",")),
              ifelse(response_is_matrix, "groups of patients (rows)", "patients")))
  cat(sprintf("Family: %s | Link function: %s\n\n", bold(family_name), bold(link_name)))
  cat("What this means:\n")
  cat(strwrap(scale_explain, width = 72, prefix = "  "), sep = "\n")
  cat("\n")
  
  if (is_negbin && !is.null(model$theta)) {
    theta_val <- model$theta
    theta_se  <- model$SE.theta
    cat(sprintf(
      "This is a NEGATIVE BINOMIAL model rather than a plain Poisson model. It estimated an\n"))
    cat(sprintf(
      "overdispersion parameter (theta) of %.2f (SE %.2f). A plain Poisson model assumes the\n",
      theta_val, ifelse(is.null(theta_se), NA, theta_se)))
    cat("variance of the counts equals their mean; needing a finite theta here means the actual\n")
    cat("data show MORE variability than a Poisson model would allow for — a common feature of\n")
    cat("real-world event-count data (e.g. hospitalisations), and the reason Negative Binomial\n")
    cat("was used instead of Poisson.\n\n")
  }
  
  
  # ==========================================================
  # 4. OVERALL MODEL FIT — does this model explain anything?
  # ==========================================================
  null_df  <- model$df.null
  res_df   <- model$df.residual
  lr_df    <- null_df - res_df
  
  if (is.null(lr_df) || lr_df <= 0) {
    cat("NOTE: this model has no predictors to test against a null (intercept-only) model,\n")
    cat("so an overall likelihood-ratio test can't be computed. Skipping to per-variable\n")
    cat("results below.\n\n")
  } else {
    lr_stat <- model$null.deviance - model$deviance
    lr_p    <- pchisq(lr_stat, lr_df, lower.tail = FALSE)
    pseudo_r2 <- tryCatch(1 - model$deviance / model$null.deviance,
                          error = function(e) NA)
    
    cat(sprintf("Overall, this model's predictors explain roughly %.1f%% of the deviance\n",
                pseudo_r2 * 100))
    cat("(a rough, non-linear-regression analogue of R-squared — useful for comparing\n")
    cat("models, less useful as a standalone 'percent explained' figure).\n\n")
    cat(sprintf("Taken together, do these variables actually help explain %s?\n", outcome))
    cat(sprintf("%s (p %s).\n\n",
                ifelse(lr_p < 0.05,
                       "Yes — the variables in this model collectively have a real, non-random relationship with the outcome",
                       "Not clearly — the variables in this model, taken together, don't show a statistically reliable relationship with the outcome"),
                format_pval(lr_p)))
  }
  
  
  # ==========================================================
  # 5. PER-VARIABLE INTERPRETATION
  # ==========================================================
  coefs_table <- s_check$coefficients
  var_names   <- rownames(coefs_table)
  pval_col    <- ncol(coefs_table)   # last column is always the p-value, regardless of
  # whether R labels it Pr(>|z|) or Pr(>|t|)
  
  term_labels <- attr(terms(model), "term.labels")  # excludes response & offset automatically
  
  cat("What each variable is associated with:\n")
  cat("--------------------------------------------------------\n")
  
  results_rows <- list()
  
  for (i in seq_along(var_names)) {
    
    vn <- var_names[i]
    if (vn == "(Intercept)") next
    if (is.na(coef_all[vn])) next  # skip non-estimable (aliased) coefficients — see diagnostic (b)
    
    estimate <- coefs_table[i, "Estimate"]
    pval     <- coefs_table[i, pval_col]
    lower_link <- if (vn %in% rownames(ci_link)) ci_link[vn, 1] else NA
    upper_link <- if (vn %in% rownames(ci_link)) ci_link[vn, 2] else NA
    
    ratio    <- exp(estimate)
    ratio_lo <- exp(lower_link)
    ratio_hi <- exp(upper_link)
    
    # ---- Classify this term ----
    is_interaction <- grepl(":", vn, fixed = TRUE)
    is_transformed <- grepl("^(log|sqrt|poly|I|exp|scale)\\(", vn)
    
    main_labels <- term_labels[!grepl(":", term_labels, fixed = TRUE)]
    matched_var <- main_labels[sapply(main_labels, function(p) startsWith(vn, p))]
    matched_var <- if (length(matched_var) > 0) matched_var[which.max(nchar(matched_var))] else character(0)
    
    is_ordered_factor <- length(matched_var) == 1 && matched_var %in% names(model$model) &&
      is.ordered(model$model[[matched_var]])
    is_factor_level <- length(matched_var) == 1 && matched_var %in% names(model$model) &&
      is.factor(model$model[[matched_var]]) && !is_ordered_factor
    
    cat(sprintf("\n> %s\n", vn))
    cat(sprintf("  %s: %.2f (%.0f%% CI: %.2f to %.2f)   [raw model coefficient: %.3f]\n",
                ratio_label, ratio, alpha_pct, ratio_lo, ratio_hi, estimate))
    
    row_type <- "main effect"
    
    if (is_interaction) {
      row_type <- "interaction"
      parts <- strsplit(vn, ":", fixed = TRUE)[[1]]
      
      cat(sprintf("  In plain English — %s:\n", bold("INTERACTION TERM")))
      if (length(parts) == 2) {
        cat(sprintf("  On this model's scale, effects combine MULTIPLICATIVELY, not additively.\n"))
        cat(sprintf("  This interaction's ratio (%.2f) tells you: the %s associated with\n",
                    ratio, sub(" \\(.*\\)", "", ratio_label)))
        cat(sprintf("  '%s' is itself multiplied by %.2f for each 1-unit change in (or shift\n",
                    parts[1], ratio))
        cat(sprintf("  to the other category of) '%s' — and symmetrically, vice versa.\n", parts[2]))
        cat(sprintf("  In short: the effect of %s DEPENDS ON the level of %s. If this ratio is\n",
                    parts[1], parts[2]))
        cat(sprintf("  close to 1, the two variables act roughly independently of each other;\n"))
        cat(sprintf("  the further from 1, the more they modify each other's effect.\n"))
      } else {
        cat(sprintf("  This is a higher-order interaction between %d variables (%s).\n",
                    length(parts), paste(parts, collapse = ", ")))
        cat(sprintf("  These are notoriously hard to interpret from a single coefficient.\n"))
        cat(sprintf("  Strongly recommend generating predicted values across combinations of\n"))
        cat(sprintf("  these variables and plotting them, rather than reading this number alone.\n"))
      }
      
    } else if (is_transformed) {
      row_type <- "transformed"
      cat(sprintf("  In plain English — %s:\n", bold("TRANSFORMED TERM")))
      cat(sprintf("  This term is a transformation of the original variable (log, polynomial,\n"))
      cat(sprintf("  or similar), not the raw variable. A '1-unit increase' or simple group\n"))
      cat(sprintf("  comparison doesn't apply cleanly on the original scale — interpret this\n"))
      cat(sprintf("  on its transformed scale, or generate predictions at specific values of\n"))
      cat(sprintf("  the original variable instead.\n"))
      
    } else if (is_ordered_factor) {
      row_type <- "ordered factor"
      cat(sprintf("  In plain English — %s:\n", bold("ORDERED FACTOR")))
      cat(sprintf("  '%s' comes from an ordered factor (e.g. Mild < Moderate < Severe). R fits\n", vn))
      cat(sprintf("  this with polynomial contrasts (trend components), not simple group\n"))
      cat(sprintf("  comparisons — this describes a TREND across the ordered levels, not a\n"))
      cat(sprintf("  'group A vs group B' difference. Re-fit with an unordered factor if you\n"))
      cat(sprintf("  want simple, directly comparable group differences instead.\n"))
      
    } else if (is_factor_level) {
      reference_level <- levels(model$model[[matched_var]])[1]
      compared_level  <- substring(vn, nchar(matched_var) + 1)
      
      # Noun to use for "higher/lower ___" so it matches the actual ratio type
      ratio_noun <- if (grepl("Odds", ratio_label)) "odds"
      else if (grepl("Risk", ratio_label)) "risk"
      else if (grepl("Incidence Rate", ratio_label)) "rate"
      else if (grepl("Rate Ratio", ratio_label)) "rate"
      else "value"
      
      # IMPORTANT: always describe the COMPARED level relative to the REFERENCE
      # level, in the SAME direction the model estimated (compared / reference).
      # Do NOT invert the ratio to describe "how much higher the reference is" —
      # percentage differences are not symmetric under inversion (a group with
      # half the odds of another is not the same magnitude as "the other group
      # has double" when expressed as a percentage), so always keep one fixed
      # direction to avoid a wrong number.
      pct_diff <- if (ratio >= 1) (ratio - 1) * 100 else (1 - ratio) * 100
      direction_word <- if (ratio > 1) "higher" else if (ratio < 1) "lower" else "no different"
      
      cat(sprintf("  In plain English — %s vs %s:\n", compared_level, reference_level))
      if (abs(ratio - 1) < 1e-9) {
        cat(sprintf("  '%s' and '%s' show essentially IDENTICAL results — no meaningful\n",
                    compared_level, reference_level))
        cat(sprintf("  difference between the two groups.\n"))
      } else {
        cat(sprintf("  %s Patients in the '%s' group have %s%% %s %s of %s than patients in\n",
                    bold("Which group is greater?"), bold(compared_level),
                    bold(sprintf("%.0f", pct_diff)), bold(direction_word), ratio_noun, bold(outcome)))
        cat(sprintf("  the '%s' (reference) group, ON AVERAGE, holding every other variable\n",
                    reference_level))
        cat(sprintf("  in the model constant.\n"))
      }
      
    } else {
      # Continuous variable
      ratio_noun <- if (grepl("Odds", ratio_label)) "odds"
      else if (grepl("Risk", ratio_label)) "risk"
      else if (grepl("Incidence Rate", ratio_label)) "rate"
      else if (grepl("Rate Ratio", ratio_label)) "rate"
      else "value"
      
      pct_diff <- if (ratio >= 1) (ratio - 1) * 100 else (1 - ratio) * 100
      direction_word <- if (ratio > 1) "higher" else if (ratio < 1) "lower" else "no different"
      
      cat(sprintf("  In plain English:\n"))
      cat(sprintf("  For every %s in %s, the %s of %s is multiplied by %s\n",
                  bold("1-unit increase"), vn, ratio_noun, outcome,
                  bold(sprintf("%.2f", ratio))))
      cat(sprintf("  — in other words, %s%% %s %s of %s, ON AVERAGE, holding every other\n",
                  bold(sprintf("%.1f", pct_diff)), bold(direction_word), ratio_noun, bold(outcome)))
      cat(sprintf("  variable in the model constant.\n"))
      if (length(matched_var) == 0) {
        cat(sprintf("  (Note: could not confidently match this term back to an original\n"))
        cat(sprintf("  variable in your data — double-check this interpretation manually.)\n"))
      }
    }
    
    cat(sprintf("\n  Is this a real effect, or could it be due to chance?\n"))
    if (is.na(pval)) {
      cat("  P-value not available for this term.\n")
    } else if (pval < 0.05) {
      cat(sprintf("  This result IS %s (p %s).\n", bold("statistically significant"), format_pval(pval)))
      cat("  In plain terms: it's unlikely (less than a 5% chance) that we'd see a ratio this\n")
      cat("  far from 1 purely by random chance if there were truly no relationship here.\n")
    } else {
      cat(sprintf("  This result is %s (p %s).\n", bold("NOT statistically significant"), format_pval(pval)))
      cat("  In plain terms: a ratio this far from 1 could plausibly happen just by random\n")
      cat("  chance, even with NO real relationship. Treat this estimate with caution.\n")
    }
    
    results_rows[[length(results_rows) + 1]] <- data.frame(
      variable   = vn,
      type       = row_type,
      ratio_label = ratio_label,
      ratio      = round(ratio, 3),
      lower_ci   = round(ratio_lo, 3),
      upper_ci   = round(ratio_hi, 3),
      p_value    = ifelse(is.na(pval), NA, ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval))),
      significant = ifelse(is.na(pval), NA, pval < 0.05)
    )
  }
  
  cat("\n--------------------------------------------------------\n")
  cat(sprintf("Reminder: all ratios above are %s — a value of 1.00 means NO difference/effect.\n",
              ratio_label))
  cat("These are associations, adjusted for the other variables in the model — not proven\n")
  cat("causal effects, and not adjusted for anything left OUT of the model.\n")
  cat("========================================================\n")
  
  results_table <- do.call(rbind, results_rows)
  
  # ==========================================================
  # 6. OPTIONAL FOREST PLOT
  # Reference/"no effect" line is drawn at 1 — the null value for
  # ANY ratio (OR, RR, or IRR). This is deliberately NOT 0: 0 would
  # only be the null value on the raw, un-exponentiated log-odds /
  # log-rate scale, which this plot does not show.
  # ==========================================================
  if (isTRUE(forest_plot)) {
    
    if (is.null(results_table) || nrow(results_table) == 0) {
      warning("forest_plot = TRUE was requested, but there are no estimable terms to plot.")
    } else {
      
      fp_data <- results_table[!is.na(results_table$ratio) &
                                 !is.na(results_table$lower_ci) &
                                 !is.na(results_table$upper_ci), ]
      
      # Log scale requires strictly positive bounds — guard against a Wald CI
      # (fallback case) that happened to produce a non-positive lower bound.
      use_log_scale <- all(fp_data$lower_ci > 0) && all(fp_data$upper_ci > 0) && all(fp_data$ratio > 0)
      
      if (nrow(fp_data) == 0) {
        warning("forest_plot = TRUE was requested, but no terms had usable ratio/CI values to plot (all were NA or non-positive).")
      } else {
        
        n_rows <- nrow(fp_data)
        y_pos  <- rev(seq_len(n_rows))
        
        x_range <- range(c(fp_data$lower_ci, fp_data$upper_ci, 1), na.rm = TRUE)
        
        old_par <- par(mar = c(5, max(8, max(nchar(fp_data$variable)) * 0.6), 4, 2))
        on.exit(par(old_par), add = TRUE)
        
        plot(fp_data$ratio, y_pos,
             xlim = x_range,
             ylim = c(0.5, n_rows + 0.5),
             log  = if (use_log_scale) "x" else "",
             pch  = 16, cex = 1.3,
             yaxt = "n", ylab = "",
             xlab = ratio_label,
             main = sprintf("Forest plot: %s\n(%s)", outcome, ratio_label))
        
        axis(2, at = y_pos, labels = fp_data$variable, las = 2, cex.axis = 0.85)
        
        segments(fp_data$lower_ci, y_pos, fp_data$upper_ci, y_pos, lwd = 2)
        
        # Reference ("no effect") line — always at 1 for a ratio scale
        abline(v = 1, lty = 2, col = "red", lwd = 1.5)
        text(x = 1, y = n_rows + 0.5, labels = "no effect", col = "red",
             pos = 3, xpd = TRUE, cex = 0.8)
        
        if (!use_log_scale) {
          warning("Some CI bounds were zero or negative (likely from a Wald-CI fallback — see diagnostics above), so this plot uses a LINEAR x-axis instead of the conventional log scale. Interpret spacing between points with extra caution.")
        }
      }
    }
  }
  
  invisible(results_table)
}


# ============================================================
# USAGE EXAMPLES (using the cohort dataset from earlier)
# ============================================================
# --- Logistic regression (binomial, logit link -> Odds Ratios) ---
model_logit <- glm(event ~ age + sex + treatment + comorbidity_score,
                    data = cohort, family = binomial(link = "logit"))
logit_results <- explain_glm(model_logit)
logit_results

#--- Poisson with offset (log link -> Incidence Rate Ratios) ---
model_pois <- glm(n_hosp ~ age + treatment + offset(log(person_time)),
                   data = cohort, family = poisson(link = "log"))
pois_results <- explain_glm(model_pois)
pois_results

#--- Negative Binomial (overdispersed counts -> IRRs, with theta note) ---
library(MASS)
model_nb <- glm.nb(n_hosp ~ age + treatment + offset(log(person_time)), data = cohort)
nb_results <- explain_glm(model_nb)
nb_results

#--- Interaction term example ---
model_int <- glm(event ~ age * treatment, data = cohort, family = binomial)
explain_glm(model_int)

#--- With a forest plot (reference line always at ratio = 1, never 0) ---
explain_glm(model_logit, forest_plot = TRUE)
explain_glm(model_pois, forest_plot = TRUE)

# ============================================================
# EDGE CASES THIS FUNCTION HAS BEEN SPECIFICALLY CHECKED AGAINST
# ============================================================
# 1. Non-convergence                          -> flagged, results still shown with caveat
# 2. Aliased/NA coefficients (collinearity)    -> flagged, those terms skipped safely
# 3. (Quasi-)complete separation (huge SEs)    -> flagged before results
# 4. Unsupported/unusual family                -> generic fallback explanation, flagged
# 5. cbind(success, failure) matrix response   -> flagged, cohort size caveat added
# 6. Non-uniform prior weights                 -> flagged, cohort size caveat added
# 7. Offset present but not log-transformed    -> flagged (IRR interpretation would be wrong)
# 8. Offset via offset() in formula OR offset= argument -> both detected
# 9. confint() failing to converge (profile CI)-> falls back to Wald CI automatically
# 10. Single-predictor models (confint returns a vector, not matrix) -> reshaped safely
# 11. Interaction terms (2-way and higher-order)-> both explained, differently
# 12. Transformed terms: log(), poly(), I(), sqrt(), scale() -> flagged, not misread
# 13. Ordered factors                          -> flagged as trend, not group comparison
# 14. Intercept-only / no predictors to test   -> overall LR test skipped gracefully
# 15. Negative Binomial via MASS::glm.nb       -> detected via class, theta explained
# 16. Terms that can't be matched to original data column -> flagged, not silently wrong
# 17. p-value column naming differences (Pr(>|z|) vs Pr(>|t|)) -> read positionally, not by name
# 18. forest_plot = TRUE with no estimable terms -> warns instead of erroring
# 19. forest_plot = TRUE where a Wald-CI fallback produced a non-positive bound
#     (breaks log-scale plotting) -> auto-switches to linear x-axis, with a warning

