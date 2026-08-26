# 03_survival_models

# ============================================================
# 3. SURVIVAL MODELS — outcome: time, status
# ============================================================
library(survival)

out_km_fit      <- survfit(Surv(time, status) ~ treatment, data = cohort)
out_km_logrank  <- survdiff(Surv(time, status) ~ treatment, data = cohort)

out_km_fit
out_km_logrank
# plot(out_km_fit, col = c("blue", "red"))  # visual, not an object

model_cox <- coxph(Surv(time, status) ~ age + sex + treatment + comorbidity_score,
                   data = cohort)
out_cox_summary <- summary(model_cox)
out_cox_zph     <- cox.zph(model_cox)   # test proportional hazards assumption

out_cox_summary
out_cox_zph


# ============================================================
# explain_km()
#
# Plain-English interpretation of Kaplan-Meier survival curves
# and optional log-rank tests.
#
# Handles:
#   - survfit()
#   - survdiff()
#   - grouped KM curves
#   - single-group KM curves
#   - median survival extraction
#   - optional survminer plot
#
# ============================================================

# ---- Bold-text helper ----
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

explain_km <- function(
    km_fit,
    logrank = NULL,
    km_plot = FALSE
){
  
  # ==========================================================
  # VALIDITY CHECKS
  # ==========================================================
  
  if (!inherits(km_fit, "survfit")) {
    stop("km_fit must be a survfit object.")
  }
  
  diagnostics <- character(0)
  
  if (any(km_fit$n < 10)) {
    
    diagnostics <- c(
      diagnostics,
      "one or more groups contain very few patients (<10), making survival estimates potentially unstable."
    )
    
  }
  
  if (length(diagnostics) > 0) {
    
    cat("########################################################\n")
    cat("# KAPLAN-MEIER DIAGNOSTIC WARNINGS\n")
    cat("########################################################\n\n")
    
    for (d in diagnostics) {
      cat(bold("WARNING:"), d, "\n\n")
    }
    
  }
  
  # ==========================================================
  # HEADER
  # ==========================================================
  
  cat("========================================================\n")
  cat("KAPLAN-MEIER SURVIVAL ANALYSIS\n")
  cat("========================================================\n\n")
  
  total_n <- sum(km_fit$n)
  
  cat(sprintf(
    "This analysis includes %s patients.\n\n",
    bold(format(total_n, big.mark = ","))
  ))
  
  cat("What this means:\n")
  cat("Kaplan-Meier analysis estimates the probability of remaining\n")
  cat("event-free over time.\n\n")
  
  cat("Interpreting the curve:\n")
  cat("  - vertical drops indicate events occurred\n")
  cat("  - flat sections indicate no events occurred\n")
  cat("  - steeper drops indicate events occurring more rapidly\n\n")
  
  # ==========================================================
  # GROUP SUMMARY
  # ==========================================================
  
  strata_names <- names(km_fit$strata)
  
  if (!is.null(strata_names)) {
    
    cat("Groups analysed:\n")
    cat("--------------------------------------------------------\n")
    
    for (i in seq_along(strata_names)) {
      
      cat(sprintf(
        "%s: %s patients\n",
        strata_names[i],
        format(km_fit$strata[i], big.mark = ",")
      ))
      
    }
    
    cat("\n")
  }
  
  # ==========================================================
  # MEDIAN SURVIVAL TIMES
  # ==========================================================
  
  km_table <- tryCatch(
    summary(km_fit)$table,
    error = function(e) NULL
  )
  
  if (!is.null(km_table)) {
    
    cat("Median survival / event-free times:\n")
    cat("--------------------------------------------------------\n")
    
    if (is.matrix(km_table)) {
      
      for (i in seq_len(nrow(km_table))) {
        
        med <- km_table[i, "median"]
        
        if (is.na(med)) {
          
          cat(sprintf(
            "%s : median survival NOT reached during follow-up\n",
            rownames(km_table)[i]
          ))
          
        } else {
          
          cat(sprintf(
            "%s : %.2f time units\n",
            rownames(km_table)[i],
            med
          ))
          
        }
        
      }
      
    } else {
      
      med <- km_table["median"]
      
      if (is.na(med)) {
        
        cat("Median survival NOT reached during follow-up.\n")
        
      } else {
        
        cat(sprintf(
          "Median survival = %.2f time units\n",
          med
        ))
        
      }
      
    }
    
    cat("\n")
    
  }
  
  # ==========================================================
  # LOG-RANK TEST
  # ==========================================================
  
  if (!is.null(logrank)) {
    
    cat("========================================================\n")
    cat("LOG-RANK TEST\n")
    cat("========================================================\n\n")
    
    chi_sq <- logrank$chisq
    df <- length(logrank$n) - 1
    p <- pchisq(chi_sq, df, lower.tail = FALSE)
    
    cat(sprintf(
      "Chi-square statistic = %.2f\n",
      chi_sq
    ))
    
    cat(sprintf(
      "Degrees of freedom = %d\n",
      df
    ))
    
    cat(sprintf(
      "P-value %s\n\n",
      format_pval(p)
    ))
    
    cat("In plain English:\n")
    
    if (p < 0.05) {
      
      cat(sprintf(
        "%s\n\n",
        bold("The survival curves differ significantly.")
      ))
      
      cat("The groups experience events at different rates\n")
      cat("over follow-up.\n\n")
      
      cat("This means survival differs between at least two\n")
      cat("of the groups being compared.\n\n")
      
    } else {
      
      cat(sprintf(
        "%s\n\n",
        bold("No statistically significant difference detected.")
      ))
      
      cat("Any observed separation between the Kaplan-Meier\n")
      cat("curves could plausibly have arisen by random chance.\n\n")
      
    }
    
  }
  
  # ==========================================================
  # OPTIONAL SURVMINER PLOT
  # ==========================================================
  
  if (isTRUE(km_plot)) {
    
    if (!requireNamespace("survminer", quietly = TRUE)) {
      
      warning(
        "Package 'survminer' is not installed. Install it with install.packages('survminer')."
      )
      
    } else {
      
      print(
        survminer::ggsurvplot(
          km_fit
        )
      )
      
    }
    
  }
  
  # ==========================================================
  # SUMMARY
  # ==========================================================
  
  cat("========================================================\n")
  cat("Reminder:\n")
  cat("Kaplan-Meier analysis is descriptive and unadjusted.\n")
  cat("Differences between curves may reflect underlying\n")
  cat("differences in patient characteristics.\n")
  cat("========================================================\n")
  
  invisible(
    list(
      survfit_object = km_fit,
      summary_table = km_table
    )
  )
  
}


# ============================================================
# explain_cox()
#
# Plain-English interpretation of Cox Proportional Hazards
# models fit using survival::coxph()
#
# Handles:
#  - Continuous variables
#  - Factors
#  - Interactions
#  - Transformed terms
#  - PH assumption checks
#  - Forest plots
#  - Hazard ratios
#
# ============================================================

bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) {
    crayon::bold(x)
  } else {
    x
  }
}

format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}

explain_cox <- function(
    model,
    zph = NULL,
    conf_level = 0.95,
    forest_plot = FALSE
){
  
  # ==========================================================
  # VALIDATION
  # ==========================================================
  
  if (!inherits(model, "coxph")) {
    stop("model must be a survival::coxph object.")
  }
  
  s <- summary(model)
  
  diagnostics <- character(0)
  
  # ==========================================================
  # DIAGNOSTICS
  # ==========================================================
  
  if (!is.null(model$fail)) {
    
    diagnostics <- c(
      diagnostics,
      paste(
        "did not converge correctly:",
        model$fail
      )
    )
  }
  
  coefs <- coef(model)
  
  if (any(abs(coefs) > 10, na.rm = TRUE)) {
    
    diagnostics <- c(
      diagnostics,
      "contains extremely large coefficients. This may indicate separation, sparse data, or unstable estimates."
    )
  }
  
  n_events <- s$nevent
  n_coef <- length(coefs)
  
  if (!is.na(n_events) &&
      n_events > 0 &&
      n_coef > 0 &&
      (n_events / n_coef) < 10) {
    
    diagnostics <- c(
      diagnostics,
      sprintf(
        "only %.1f events per coefficient were available. Hazard ratios may be unstable.",
        n_events / n_coef
      )
    )
  }
  
  if (!is.null(zph)) {
    
    zph_tab <- zph$table
    
    bad_terms <- rownames(zph_tab)[
      zph_tab[, "p"] < 0.05
    ]
    
    bad_terms <- setdiff(
      bad_terms,
      "GLOBAL"
    )
    
    if (length(bad_terms) > 0) {
      
      diagnostics <- c(
        diagnostics,
        paste0(
          "possible proportional hazards violations detected for: ",
          paste(bad_terms, collapse = ", "),
          ". Hazard ratios for these variables may change over time."
        )
      )
    }
    
    if ("GLOBAL" %in% rownames(zph_tab)) {
      
      gp <- zph_tab["GLOBAL", "p"]
      
      if (!is.na(gp) && gp < 0.05) {
        
        diagnostics <- c(
          diagnostics,
          "global proportional hazards test was statistically significant."
        )
      }
    }
  }
  
  if (length(diagnostics) > 0) {
    
    cat("########################################################\n")
    cat("# MODEL DIAGNOSTIC WARNINGS\n")
    cat("########################################################\n\n")
    
    for (d in diagnostics) {
      
      cat(
        bold("WARNING:"),
        "this model",
        d,
        "\n\n"
      )
    }
  }
  
  # ==========================================================
  # HEADER
  # ==========================================================
  
  n_obs <- s$n
  
  cat("========================================================\n")
  cat("COX PROPORTIONAL HAZARDS MODEL\n")
  cat("========================================================\n\n")
  
  cat(
    "MODEL:",
    deparse(formula(model)),
    "\n\n"
  )
  
  cat(sprintf(
    "This model was fit using data from %s patients.\n",
    bold(format(n_obs, big.mark = ","))
  ))
  
  cat(sprintf(
    "%s events occurred during follow-up.\n\n",
    format(n_events, big.mark = ",")
  ))
  
  cat("What this means:\n")
  cat("The outcome is TIME UNTIL AN EVENT OCCURS.\n")
  cat("The model estimates the HAZARD of the event.\n\n")
  
  cat("Hazard Ratio (HR) interpretation:\n")
  cat("  HR = 1.0 -> no difference\n")
  cat("  HR > 1.0 -> higher hazard\n")
  cat("  HR < 1.0 -> lower hazard\n\n")
  
  # ==========================================================
  # OVERALL MODEL TEST
  # ==========================================================
  
  lr <- s$logtest
  
  if (!is.null(lr)) {
    
    cat("========================================================\n")
    cat("OVERALL MODEL FIT\n")
    cat("========================================================\n\n")
    
    cat(sprintf(
      "Likelihood ratio test: Chi-square %.2f on %d df\n",
      lr["test"],
      lr["df"]
    ))
    
    cat(sprintf(
      "P-value %s\n\n",
      format_pval(lr["pvalue"])
    ))
    
    if (lr["pvalue"] < 0.05) {
      
      cat(
        "Taken together, the variables in this model\n",
        "show a statistically significant relationship\n",
        "with time-to-event.\n\n"
      )
      
    } else {
      
      cat(
        "Taken together, the variables do not show a\n",
        "clear statistically significant relationship\n",
        "with time-to-event.\n\n"
      )
    }
  }
  
  # ==========================================================
  # PH ASSUMPTION SECTION
  # ==========================================================
  
  if (!is.null(zph)) {
    
    cat("========================================================\n")
    cat("PROPORTIONAL HAZARDS ASSUMPTION\n")
    cat("========================================================\n\n")
    
    ztab <- zph$table
    
    for (i in seq_len(nrow(ztab))) {
      
      nm <- rownames(ztab)[i]
      
      p <- ztab[i, "p"]
      
      if (nm == "GLOBAL") next
      
      if (p < 0.05) {
        
        cat(sprintf(
          "%s : possible violation (p %s)\n",
          nm,
          format_pval(p)
        ))
        
      } else {
        
        cat(sprintf(
          "%s : no evidence of violation (p %s)\n",
          nm,
          format_pval(p)
        ))
      }
    }
    
    if ("GLOBAL" %in% rownames(ztab)) {
      
      gp <- ztab["GLOBAL", "p"]
      
      cat("\nGLOBAL TEST:\n")
      
      if (gp < 0.05) {
        
        cat(sprintf(
          "Overall PH assumption may be violated (p %s)\n\n",
          format_pval(gp)
        ))
        
      } else {
        
        cat(sprintf(
          "No evidence of overall PH violation (p %s)\n\n",
          format_pval(gp)
        ))
      }
    }
  }
  
  # ==========================================================
  # VARIABLE INTERPRETATION
  # ==========================================================
  
  coef_table <- s$coefficients
  ci_table <- s$conf.int
  
  term_labels <- attr(
    terms(model),
    "term.labels"
  )
  
  results_rows <- list()
  
  cat("========================================================\n")
  cat("VARIABLE INTERPRETATION\n")
  cat("========================================================\n")
  
  for (i in seq_len(nrow(coef_table))) {
    
    vn <- rownames(coef_table)[i]
    
    hr <- ci_table[i, "exp(coef)"]
    lower <- ci_table[i, "lower .95"]
    upper <- ci_table[i, "upper .95"]
    
    p <- coef_table[i, "Pr(>|z|)"]
    
    cat("\n--------------------------------------------------------\n")
    cat(vn, "\n\n")
    
    cat(sprintf(
      "Hazard Ratio: %.2f (%.0f%% CI %.2f to %.2f)\n\n",
      hr,
      conf_level * 100,
      lower,
      upper
    ))
    
    is_interaction <- grepl(":", vn)
    
    is_transformed <- grepl(
      "^(log|sqrt|poly|I|scale|exp)\\(",
      vn
    )
    
    main_labels <- term_labels[
      !grepl(":", term_labels)
    ]
    
    matched_var <- main_labels[
      sapply(main_labels,
             function(x) startsWith(vn, x))
    ]
    
    matched_var <- if(length(matched_var)>0)
      matched_var[which.max(nchar(matched_var))]
    else
      character(0)
    
    if (is_interaction) {
      
      cat(
        bold("INTERACTION TERM"),
        "\n\n"
      )
      
      cat(
        "The effect of one variable depends on the\n",
        "level or value of another variable.\n",
        "Consider plotting predicted survival curves\n",
        "for easier interpretation.\n"
      )
      
    } else if (is_transformed) {
      
      cat(
        bold("TRANSFORMED TERM"),
        "\n\n"
      )
      
      cat(
        "This coefficient relates to a transformed\n",
        "version of the original variable.\n",
        "Interpretation on the raw scale is not direct.\n"
      )
      
    } else {
      
      if (
        length(matched_var) == 1 &&
        matched_var %in% names(model$model) &&
        is.factor(model$model[[matched_var]])
      ) {
        
        ref <- levels(model$model[[matched_var]])[1]
        
        cmp <- substring(
          vn,
          nchar(matched_var) + 1
        )
        
        pct <- abs(hr - 1) * 100
        
        cat(sprintf(
          "%s vs %s\n\n",
          cmp,
          ref
        ))
        
        if (hr > 1) {
          
          cat(sprintf(
            "Patients in '%s' have %.0f%% higher hazard\n",
            cmp,
            pct
          ))
          
          cat(
            "of experiencing the event at any point\n",
            "during follow-up compared with patients\n",
            "in the reference group.\n"
          )
          
        } else {
          
          cat(sprintf(
            "Patients in '%s' have %.0f%% lower hazard\n",
            cmp,
            pct
          ))
          
          cat(
            "of experiencing the event at any point\n",
            "during follow-up compared with patients\n",
            "in the reference group.\n"
          )
        }
        
      } else {
        
        pct <- abs(hr - 1) * 100
        
        cat(
          "Continuous variable interpretation:\n\n"
        )
        
        if (hr > 1) {
          
          cat(sprintf(
            "For every 1-unit increase in %s,\n",
            vn
          ))
          
          cat(sprintf(
            "hazard increases by approximately %.1f%%.\n",
            pct
          ))
          
        } else {
          
          cat(sprintf(
            "For every 1-unit increase in %s,\n",
            vn
          ))
          
          cat(sprintf(
            "hazard decreases by approximately %.1f%%.\n",
            pct
          ))
        }
      }
    }
    
    cat("\n")
    
    if (p < 0.05) {
      
      cat(sprintf(
        "Statistically significant (p %s)\n",
        format_pval(p)
      ))
      
    } else {
      
      cat(sprintf(
        "Not statistically significant (p %s)\n",
        format_pval(p)
      ))
    }
    
    results_rows[[length(results_rows)+1]] <- data.frame(
      variable = vn,
      hr = round(hr,3),
      lower_ci = round(lower,3),
      upper_ci = round(upper,3),
      p_value = ifelse(
        p < 0.001,
        "< 0.001",
        sprintf("%.3f", p)
      ),
      significant = p < 0.05
    )
  }
  
  results_table <- do.call(
    rbind,
    results_rows
  )
  
  # ==========================================================
  # FOREST PLOT
  # ==========================================================
  
  if (isTRUE(forest_plot)) {
    
    fp <- results_table
    
    y_pos <- rev(seq_len(nrow(fp)))
    
    x_range <- range(
      c(
        fp$lower_ci,
        fp$upper_ci,
        1
      ),
      na.rm = TRUE
    )
    
    old_par <- par(
      mar = c(
        5,
        max(
          8,
          max(nchar(fp$variable))*0.6
        ),
        4,
        2
      )
    )
    
    on.exit(par(old_par),
            add = TRUE)
    
    plot(
      fp$hr,
      y_pos,
      pch = 16,
      cex = 1.2,
      log = "x",
      xlim = x_range,
      yaxt = "n",
      ylab = "",
      xlab = "Hazard Ratio (HR)",
      main = "Forest Plot: Hazard Ratios"
    )
    
    axis(
      2,
      at = y_pos,
      labels = fp$variable,
      las = 2
    )
    
    segments(
      fp$lower_ci,
      y_pos,
      fp$upper_ci,
      y_pos,
      lwd = 2
    )
    
    abline(
      v = 1,
      col = "red",
      lty = 2,
      lwd = 2
    )
    
    text(
      x = 1,
      y = max(y_pos)+0.8,
      labels = "No Effect",
      col = "red",
      pos = 3
    )
  }
  
  cat("\n========================================================\n")
  cat("Reminder: Hazard ratios describe associations,\n")
  cat("not proven causal effects.\n")
  cat("========================================================\n")
  
  invisible(results_table)
}


# Example usage

library(survival)

out_km_fit <- survfit(
  Surv(time, status) ~ treatment,
  data = cohort
)

out_km_logrank <- survdiff(
  Surv(time, status) ~ treatment,
  data = cohort
)

explain_km(
  km_fit = out_km_fit,
  logrank = out_km_logrank,
  km_plot = TRUE
)

# another nice way to do this
library(survminer)
ggsurvplot(out_km_fit)



library(survival)

model_cox <- coxph(
  Surv(time, status) ~ age + sex + treatment + comorbidity_score,
  data = cohort
)

out_cox_zph <- cox.zph(model_cox)

explain_cox(
  model = model_cox,
  zph = out_cox_zph,
  forest_plot = TRUE
)
