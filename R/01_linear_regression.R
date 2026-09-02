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

# ============================================================
# explain_lm() — plain-language interpretation of an lm() model
# ============================================================

# ---- Bold-text helper ----
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    crayon::bold(x)
  } else {
    x
  }
}

# ---- Helper: format p-values sensibly ----
format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}

explain_lm <- function(model, conf_level = 0.95) {
  
  if (!inherits(model, "lm")) {
    stop("explain_lm() expects an object from lm()")
  }
  
  outcome_data <- model.response(model.frame(model))
  
  if (!is.numeric(outcome_data)) {
    warning(
      paste0(
        "The outcome variable does not appear numeric. ",
        "lm() may have converted a factor outcome to numeric codes. ",
        "If this is a binary/categorical outcome, consider glm(..., family = binomial)."
      )
    )
  }
  
  s <- summary(model)
  
  ci <- tryCatch(
    confint(model, level = conf_level),
    error = function(e) NULL
  )
  
  coefs <- s$coefficients
  var_names <- rownames(coefs)
  
  outcome <- deparse(formula(model)[[2]])
  alpha_pct <- conf_level * 100
  
  cat("========================================================\n")
  cat("MODEL:", deparse(formula(model)), "\n")
  cat("Outcome variable:", outcome, "\n")
  cat("========================================================\n\n")
  
  # ---------------------------------------------------------
  # Overall model fit
  # ---------------------------------------------------------
  
  r2 <- s$r.squared
  adj_r2 <- s$adj.r.squared
  fstat <- s$fstatistic
  n_obs <- length(model$residuals)
  
  if (is.null(fstat)) {
    
    cat("NOTE: this model has no predictors to test, or has zero residual\n")
    cat("degrees of freedom (perfect/saturated fit).\n")
    cat("Overall model statistics cannot be meaningfully calculated.\n\n")
    
    f_pval <- NA
    
  } else {
    
    f_pval <- pf(
      fstat[1],
      fstat[2],
      fstat[3],
      lower.tail = FALSE
    )
    
    cat(sprintf(
      "Overall, this model accounts for %.1f%% of the observed variation in %s\n",
      r2 * 100,
      outcome
    ))
    
    cat(sprintf(
      "(adjusted R² = %.1f%%) using %d observations.\n\n",
      adj_r2 * 100,
      n_obs
    ))
    
    cat(sprintf(
      "Taken together, do these variables help explain variation in %s?\n",
      outcome
    ))
    
    cat(sprintf(
      "%s (p %s).\n\n",
      ifelse(
        f_pval < 0.05,
        "Yes — the predictors collectively show evidence of association with the outcome",
        "Not clearly — the predictors do not collectively show strong statistical evidence of association with the outcome"
      ),
      format_pval(f_pval)
    ))
  }
  
  # ---------------------------------------------------------
  # Per-variable interpretation
  # ---------------------------------------------------------
  
  cat("What each variable is associated with:\n")
  cat("--------------------------------------------------------\n")
  
  predictor_names <- names(model$model)[-1]
  
  for (i in seq_along(var_names)) {
    
    vn <- var_names[i]
    
    if (vn == "(Intercept)") {
      next
    }
    
    estimate <- coefs[i, "Estimate"]
    pval <- coefs[i, "Pr(>|t|)"]
    
    if (!is.null(ci) && vn %in% rownames(ci)) {
      lower <- ci[vn, 1]
      upper <- ci[vn, 2]
    } else {
      lower <- NA
      upper <- NA
    }
    
    direction <- ifelse(estimate > 0, "increase", "decrease")
    
    # -------------------------------------------------------
    # Term classification
    # -------------------------------------------------------
    
    is_interaction <- grepl(":", vn, fixed = TRUE)
    
    is_transformed <- grepl("\\(", vn)
    
    matched_var <- predictor_names[
      sapply(
        predictor_names,
        function(p) startsWith(vn, p)
      )
    ]
    
    matched_var <- if (length(matched_var) > 0) {
      matched_var[which.max(nchar(matched_var))]
    } else {
      character(0)
    }
    
    is_ordered_factor <-
      length(matched_var) == 1 &&
      is.ordered(model$model[[matched_var]])
    
    is_factor_level <-
      length(matched_var) == 1 &&
      is.factor(model$model[[matched_var]]) &&
      !is_ordered_factor
    
    cat(sprintf("\n> %s\n", vn))
    
    if (!is.na(lower) && !is.na(upper)) {
      cat(sprintf(
        "  Raw model coefficient: %.3f (%.0f%% CI: %.3f to %.3f)\n",
        estimate,
        alpha_pct,
        lower,
        upper
      ))
    } else {
      cat(sprintf(
        "  Raw model coefficient: %.3f\n",
        estimate
      ))
    }
    
    # -------------------------------------------------------
    # Interpretation
    # -------------------------------------------------------
    
    if (is_interaction) {
      
      cat("  In plain English:\n")
      cat(sprintf(
        "  This is an %s between variables.\n",
        bold("INTERACTION TERM")
      ))
      cat("  It describes how the effect of one variable changes\n")
      cat("  depending on the value of another variable.\n")
      cat("  Interpret alongside the corresponding main effects,\n")
      cat("  ideally using predicted values or plots.\n")
      
    } else if (is_transformed) {
      
      cat("  In plain English:\n")
      cat(sprintf(
        "  This is a %s.\n",
        bold("TRANSFORMED/NON-LINEAR TERM")
      ))
      cat("  The coefficient is expressed on the transformed scale,\n")
      cat("  so a simple '1-unit increase' interpretation may be misleading.\n")
      cat("  Consider generating predictions at meaningful values instead.\n")
      
    } else if (is_ordered_factor) {
      
      cat("  In plain English:\n")
      cat(sprintf(
        "  '%s' arises from an %s.\n",
        vn,
        bold("ORDERED FACTOR")
      ))
      cat("  Ordered factors are often represented using trend-based\n")
      cat("  contrasts rather than simple group comparisons.\n")
      cat("  This coefficient therefore reflects an underlying trend\n")
      cat("  across ordered categories rather than a single group contrast.\n")
      
    } else if (is_factor_level) {
      
      cat("  In plain English:\n")
      
      cat(sprintf(
        "  Compared with the reference group, patients in '%s'\n",
        vn
      ))
      
      cat(sprintf(
        "  have, on average, a %s of %s units in %s,\n",
        bold(direction),
        bold(sprintf("%.3f", abs(estimate))),
        outcome
      ))
      
      cat("  assuming all other variables in the model remain the same.\n")
      
    } else {
      
      cat("  In plain English:\n")
      
      cat(sprintf(
        "  For every %s in %s,\n",
        bold("1-unit increase"),
        vn
      ))
      
      cat(sprintf(
        "  %s changes by %s units on average,\n",
        outcome,
        bold(sprintf("%.3f", estimate))
      ))
      
      cat(sprintf(
        "  representing a %s.\n",
        bold(direction)
      ))
      
      cat("  This interpretation assumes all other variables in the model\n")
      cat("  remain unchanged.\n")
      
      if (length(matched_var) == 0) {
        cat("  (Note: this term could not be confidently matched back\n")
        cat("  to an original variable name. Check manually.)\n")
      }
    }
    
    # -------------------------------------------------------
    # Confidence interval
    # -------------------------------------------------------
    
    if (!is.na(lower) && !is.na(upper)) {
      
      cat(sprintf(
        "\n  %.0f%% confidence interval: %.3f to %.3f.\n",
        alpha_pct,
        lower,
        upper
      ))
    }
    
    # -------------------------------------------------------
    # Significance
    # -------------------------------------------------------
    
    cat("\n  Statistical evidence:\n")
    
    if (is.na(pval)) {
      
      cat("  The p-value could not be calculated reliably.\n")
      cat("  This often occurs because of singularities,\n")
      cat("  multicollinearity, or insufficient information.\n")
      
    } else if (pval < 0.05) {
      
      cat(sprintf(
        "  This result is %s (p %s).\n",
        bold("statistically significant"),
        format_pval(pval)
      ))
      
      cat("  This provides evidence that the predictor is associated\n")
      cat("  with the outcome in this dataset.\n")
      
    } else {
      
      cat(sprintf(
        "  This result is %s (p %s).\n",
        bold("NOT statistically significant"),
        format_pval(pval)
      ))
      
      cat("  The data do not provide strong statistical evidence\n")
      cat("  for an association between this predictor and the outcome.\n")
    }
  }
  
  cat("\n--------------------------------------------------------\n")
  cat("Reminder: these are associations, not proof of causation.\n")
  cat("Coefficients are adjusted for other variables included\n")
  cat("in the model, but not for omitted variables.\n")
  cat("========================================================\n")
  
  keep_rows <- var_names != "(Intercept)"
  
  if (!any(keep_rows)) {
    return(invisible(data.frame()))
  }
  
  raw_pvals <- coefs[keep_rows, "Pr(>|t|)"]
  
  invisible(
    data.frame(
      variable = var_names[keep_rows],
      estimate = round(coefs[keep_rows, "Estimate"], 3),
      lower_ci = if (!is.null(ci))
        round(ci[keep_rows, 1], 3)
      else
        NA,
      upper_ci = if (!is.null(ci))
        round(ci[keep_rows, 2], 3)
      else
        NA,
      p_value = sapply(raw_pvals, format_pval),
      significant = raw_pvals < 0.05,
      row.names = NULL
    )
  )
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
model_ord <- lm(sbp ~ severity, data = cohort)
explain_lm(model_ord)
