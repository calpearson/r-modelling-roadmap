# 01_linear_regression

# This is a multiple linear regression asking: 
# how do age, BMI, treatment, and comorbidity burden each independently relate to systolic blood pressure (SBP)?


# ============================================================
# 1. LINEAR REGRESSION — outcome: sbp
# ============================================================
model_lm <- lm(sbp ~ age + bmi + treatment + comorbidity_score, data = cohort)

out_lm_summary <- summary(model_lm)
out_lm_confint <- confint(model_lm)
out_lm_predict <- predict(model_lm,
                          newdata = data.frame(age = 60, bmi = 28, treatment = "Drug",
                                               comorbidity_score = 5),
                          interval = "confidence")

out_lm_summary
out_lm_confint
out_lm_predict

# par(mfrow = c(2, 2)); plot(model_lm); par(mfrow = c(1, 1))  # diagnostic plots (visual, not an object)


# ============================================================
# explain_lm() — plain-language interpretation of an lm() model
# Works for ANY number/combination of predictors (continuous
# or categorical), no need to edit the function per-model.
# ============================================================

# ---- Bold-text helper ----
# Uses the crayon package if installed (auto-detects whether your terminal
# supports it and prints plain text if not). Falls back to plain text if
# crayon isn't installed at all, so the function still runs either way.
# install.packages("crayon")  # run once, if you don't have it
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    crayon::bold(x)
  } else {
    x
  }
}

explain_lm <- function(model, conf_level = 0.95) {
  
  if (!inherits(model, "lm")) stop("explain_lm() expects an object from lm()")
  
  outcome_data <- model.response(model.frame(model))
  if (!is.numeric(outcome_data)) {
    warning(paste0(
      "The outcome variable in this model does not look numeric (it's a factor ",
      "or character). lm() will have silently converted it to numeric codes ",
      "(e.g. 1, 2, 3) behind the scenes, which usually isn't what you want for ",
      "a categorical outcome. Consider glm(..., family = binomial) instead if ",
      "this is a yes/no outcome. Proceeding anyway, but treat the numbers below with caution."
    ))
  }
  
  s          <- summary(model)
  ci         <- confint(model, level = conf_level)
  coefs      <- s$coefficients
  var_names  <- rownames(coefs)
  outcome    <- deparse(formula(model)[[2]])  # preserves log(), sqrt(), etc. instead of stripping them
  alpha_pct  <- conf_level * 100
  
  cat("========================================================\n")
  cat("MODEL:", deparse(formula(model)), "\n")
  cat("Outcome variable:", outcome, "\n")
  cat("========================================================\n\n")
  
  # ---- Overall model fit ----
  r2      <- s$r.squared
  adj_r2  <- s$adj.r.squared
  fstat   <- s$fstatistic
  n_obs   <- length(model$residuals)
  
  if (is.null(fstat)) {
    # Happens with intercept-only models (no predictors) or saturated models
    # (residual degrees of freedom = 0, e.g. n observations = n parameters).
    cat("NOTE: this model has no predictors to test, or has zero residual\n")
    cat("degrees of freedom (perfect/saturated fit). An overall model-fit\n")
    cat("test and R-squared can't be meaningfully computed here — skipping\n")
    cat("straight to the per-variable summary below, if any predictors exist.\n\n")
    f_pval <- NA
  } else {
    f_pval  <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
    
    cat(sprintf("Overall, this model explains %.1f%% of why %s varies across patients\n",
                r2 * 100, outcome))
    cat(sprintf("in this dataset (%d patients used). The remaining %.1f%% is due to factors\n",
                n_obs, 100 - r2 * 100))
    cat(sprintf("not captured by this model (or random variation).\n\n"))
    cat(sprintf("Taken together, do these variables actually help explain %s?\n", outcome))
    cat(sprintf("%s (p %s).\n\n",
                ifelse(f_pval < 0.05,
                       "Yes — the variables in this model collectively have a real, non-random relationship with the outcome",
                       "Not clearly — the variables in this model, taken together, don't show a statistically reliable relationship with the outcome"),
                format_pval(f_pval)))
  }
  
  # ---- Per-variable interpretation ----
  cat("What each variable is associated with:\n")
  cat("--------------------------------------------------------\n")
  
  for (i in seq_along(var_names)) {
    
    vn <- var_names[i]
    if (vn == "(Intercept)") next  # skip intercept, rarely of direct interest
    
    estimate <- coefs[i, "Estimate"]
    pval     <- coefs[i, "Pr(>|t|)"]
    lower    <- ci[vn, 1]
    upper    <- ci[vn, 2]
    
    direction <- ifelse(estimate > 0, "increase", "decrease")
    
    # ---- Classify this term: interaction / transformed / ordered factor / factor / continuous ----
    # Checked in this order because a term can only be safely called "a factor
    # group" or "a 1-unit continuous change" once the trickier cases are ruled out.
    
    is_interaction <- grepl(":", vn, fixed = TRUE)
    is_transformed <- grepl("^(log|sqrt|poly|I|exp|scale)\\(", vn)
    
    predictor_names <- names(model$model)[-1]  # drop the outcome column
    matched_var <- predictor_names[sapply(predictor_names, function(p) startsWith(vn, p))]
    matched_var <- if (length(matched_var) > 0) matched_var[which.max(nchar(matched_var))] else character(0)
    
    is_ordered_factor <- length(matched_var) == 1 && is.ordered(model$model[[matched_var]])
    is_factor_level    <- length(matched_var) == 1 && is.factor(model$model[[matched_var]]) &&
      !is_ordered_factor
    
    cat(sprintf("\n> %s\n", vn))
    cat(sprintf("  Raw model coefficient: %.3f (%.0f%% CI: %.3f to %.3f)\n",
                estimate, alpha_pct, lower, upper))
    
    if (is_interaction) {
      cat(sprintf("  In plain English:\n"))
      cat(sprintf("  This is an %s between two (or more) variables — it captures\n", bold("INTERACTION TERM")))
      cat(sprintf("  whether the effect of one variable DEPENDS ON the level of the other.\n"))
      cat(sprintf("  A simple 'X unit increase' or 'group A vs B' statement doesn't apply\n"))
      cat(sprintf("  cleanly here — interpret this alongside its component main-effect terms,\n"))
      cat(sprintf("  ideally by plotting predicted values across combinations of both variables.\n"))
      
    } else if (is_transformed) {
      cat(sprintf("  In plain English:\n"))
      cat(sprintf("  This term is a %s of the original variable (e.g. log, polynomial,\n", bold("TRANSFORMED/NON-LINEAR TERM")))
      cat(sprintf("  or a custom transformation), not the raw variable itself. A '1-unit\n"))
      cat(sprintf("  increase' statement would be misleading, because a 1-unit change in the\n"))
      cat(sprintf("  transformed scale does NOT correspond to a 1-unit change in the original\n"))
      cat(sprintf("  variable. Interpret this coefficient on its transformed scale, or generate\n"))
      cat(sprintf("  predictions at specific values of the original variable instead.\n"))
      
    } else if (is_ordered_factor) {
      cat(sprintf("  In plain English:\n"))
      cat(sprintf("  '%s' comes from an %s (e.g. Mild < Moderate < Severe).\n", vn, bold("ORDERED FACTOR")))
      cat(sprintf("  R fits this with polynomial contrasts (.L = linear trend, .Q = quadratic\n"))
      cat(sprintf("  trend, etc.) rather than simple group comparisons, so this is NOT a\n"))
      cat(sprintf("  'group A vs group B' effect — it describes a trend ACROSS the ordered\n"))
      cat(sprintf("  levels. Consider re-running with the variable as a plain (unordered) factor\n"))
      cat(sprintf("  if you want simple, directly comparable group differences instead.\n"))
      
    } else if (is_factor_level) {
      cat(sprintf("  In plain English:\n"))
      cat(sprintf("  Patients in the '%s' group have, ON AVERAGE, %s of %s in %s\n",
                  vn, bold(direction), bold(sprintf("%.3f", abs(estimate))), outcome))
      cat(sprintf("  compared to patients in the reference group — assuming every other\n"))
      cat(sprintf("  variable in the model (e.g. age, BMI, etc.) is the SAME for both groups.\n"))
      cat(sprintf("  (i.e. this isolates the effect of group membership alone, not age/BMI/etc.)\n"))
    } else {
      cat(sprintf("  In plain English:\n"))
      cat(sprintf("  For every %s in %s (e.g. one extra year, if this is age),\n",
                  bold("1-unit increase"), vn))
      cat(sprintf("  %s changes by %s units, ON AVERAGE — a %s.\n",
                  outcome, bold(sprintf("%.3f", estimate)), bold(direction)))
      cat(sprintf("  This assumes every other variable in the model stays the same\n"))
      cat(sprintf("  (i.e. this isolates the effect of %s alone, not the other variables).\n", vn))
      if (length(matched_var) == 0) {
        cat(sprintf("  (Note: could not confidently match this term back to an original\n"))
        cat(sprintf("  variable in your data — double-check this interpretation manually.)\n"))
      }
    }
    
    cat(sprintf("\n  How confident are we? We are %.0f%% confident the TRUE effect lies\n",
                alpha_pct))
    cat(sprintf("  somewhere between %.3f and %.3f units of %s.\n", lower, upper, outcome))
    
    cat(sprintf("\n  Is this a real effect, or could it be due to chance?\n"))
    if (pval < 0.05) {
      cat(sprintf("  This result IS %s (p %s).\n", bold("statistically significant"), format_pval(pval)))
      cat(sprintf("  In plain terms: it's unlikely (less than a 5%% chance) that we'd see\n"))
      cat(sprintf("  an effect this large purely by random chance if %s had NO real\n", vn))
      cat(sprintf("  relationship with %s. So this is probably a genuine pattern in the data.\n", outcome))
    } else {
      cat(sprintf("  This result is %s (p %s).\n", bold("NOT statistically significant"), format_pval(pval)))
      cat(sprintf("  In plain terms: an effect this size could plausibly happen just by\n"))
      cat(sprintf("  random chance, even if %s has NO real relationship with %s.\n", vn, outcome))
      cat(sprintf("  Treat this estimate with caution — we can't be confident it's a real effect.\n"))
    }
    
  }
  
  cat("\n--------------------------------------------------------\n")
  cat("Reminder: these are associations, not proven causal effects.\n")
  cat("Each coefficient is 'adjusted' for the other variables in the model,\n")
  cat("but not for anything left out of the model.\n")
  cat("========================================================\n")
  
  # Return a tidy data frame silently, in case you want to use it programmatically
  raw_pvals <- coefs[var_names != "(Intercept)", "Pr(>|t|)"]
  invisible(data.frame(
    variable   = var_names[var_names != "(Intercept)"],
    estimate   = round(coefs[var_names != "(Intercept)", "Estimate"], 3),
    lower_ci   = round(ci[var_names != "(Intercept)", 1], 3),
    upper_ci   = round(ci[var_names != "(Intercept)", 2], 3),
    p_value    = sapply(raw_pvals, function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))),
    significant = raw_pvals < 0.05,
    row.names  = NULL
  ))
}

# Helper: format p-values sensibly (avoids "p = 0.000000012")
format_pval <- function(p) {
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}


# ============================================================
# USAGE EXAMPLE (using the cohort dataset / model_lm from before)
# ============================================================
model_lm <- lm(sbp ~ age + bmi + treatment + comorbidity_score, data = cohort)
explain_lm(model_lm)

#Works identically if you add/remove variables, e.g.:
model_lm2 <- lm(sbp ~ age + bmi + treatment + comorbidity_score + sex, data = cohort)
explain_lm(model_lm2)

#You can also capture the tidy summary table it returns:
lm_results_table <- explain_lm(model_lm)
lm_results_table
#
# ============================================================
# EXAMPLES OF THE EDGE CASES THIS FUNCTION NOW HANDLES SAFELY
# ============================================================
#Interaction term — now explicitly flagged, not misread as plain "age":
model_int <- lm(sbp ~ age * treatment, data = cohort)
explain_lm(model_int)

#Transformed term — now explicitly flagged, not given a false "1-unit" claim:
model_log <- lm(log(sbp) ~ age + treatment, data = cohort)
explain_lm(model_log)

#Intercept-only model — now skips the F-test/R2 section instead of erroring:
model_null <- lm(sbp ~ 1, data = cohort)
explain_lm(model_null)

#Ordered factor — now flagged as a trend across levels, not a group comparison:
cohort$severity <- factor(sample(c("Mild","Moderate","Severe"), nrow(cohort), replace = TRUE),
                           ordered = TRUE, levels = c("Mild","Moderate","Severe"))
#model_ord <- lm(sbp ~ severity, data = cohort)
explain_lm(model_ord)
