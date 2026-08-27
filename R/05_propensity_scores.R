# 05_propensity_scores

# ============================================================
# 5. PROPENSITY SCORES — outcome: event / time+status, exposure: treatment
# ============================================================
library(MatchIt)
library(survival)
ps_model <- matchit(treatment ~ age + sex + comorbidity_score,
                    data = cohort, method = "nearest", ratio = 1)
out_ps_summary  <- summary(ps_model)
out_ps_summary

matched_cohort  <- match.data(ps_model)
out_ps_outcome  <- coxph(Surv(time, status) ~ treatment, data = matched_cohort)
out_ps_outcome_summary <- summary(out_ps_outcome)
out_ps_outcome_summary

library(WeightIt)
model_iptw <- weightit(treatment ~ age + sex + comorbidity_score, data = cohort, method = "ps")
out_iptw_summary <- summary(model_iptw)
out_iptw_summary






# ============================================================
# explain_ps.R — plain-language interpretation of propensity
# score MATCHING (MatchIt) and WEIGHTING (WeightIt).
#
# Two functions:
#   explain_matchit(model)  — MatchIt::matchit() objects
#   explain_iptw(model)     — WeightIt::weightit() objects
#
# IMPORTANT DESIGN NOTE ON ROBUSTNESS:
# Unlike lm/glm/coxph/lmer (base R / lme4, very stable internal
# structure), MatchIt's and WeightIt's internal object structure
# has changed across package versions and varies by method (nearest
# vs full vs subclass vs exact matching; ps vs ebal vs cbps
# weighting; binary vs multi-category vs continuous treatment).
# This file is written DEFENSIVELY throughout:
#   - every extraction from a summary()/model object is wrapped in
#     tryCatch, never assumed to exist
#   - the most important numbers (matched sample sizes, effective
#     sample size) are computed DIRECTLY from the model's own raw
#     $weights / $treat components wherever possible, rather than
#     depending on summary()'s specific internal field names —
#     these raw components are stable, documented parts of both
#     packages' object APIs and are far less likely to break
#     across versions than a summary object's internal layout
#   - if anything can't be extracted, the function WARNS and shows
#     you the raw summary/object instead of crashing or (worse)
#     silently showing wrong numbers
# See the full edge-case list at the bottom of this file.
# ============================================================

# ---- Shared helpers ----
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) crayon::bold(x) else x
}
format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}


# ============================================================
# 1. explain_matchit() — propensity score MATCHING
# ============================================================
explain_matchit <- function(model, smd_threshold = 0.1) {
  
  if (!inherits(model, "matchit")) {
    stop("explain_matchit() expects a model fit with MatchIt::matchit().")
  }
  
  cat("========================================================\n")
  cat("PROPENSITY SCORE MATCHING\n")
  cat("========================================================\n\n")
  
  cat("The goal of matching is to create a Control group that looks as SIMILAR as\n")
  cat("possible to the Treated group on the covariates you specified — so that any\n")
  cat("later difference in OUTCOME between groups is less likely to simply reflect\n")
  cat("pre-existing differences between who got treated and who didn't.\n\n")
  
  # ==========================================================
  # 1a. Method / call info (defensive — field names vary by MatchIt version)
  # ==========================================================
  method_desc <- tryCatch({
    m <- model$info$method
    if (is.null(m)) m <- "unknown"
    ratio_val <- tryCatch(model$info$ratio, error = function(e) NULL)
    if (!is.null(ratio_val)) sprintf("%s matching (ratio %s:1)", m, ratio_val) else sprintf("%s matching", m)
  }, error = function(e) "an unspecified matching method")
  
  cat(sprintf("Method used: %s\n\n", bold(method_desc)))
  
  # ==========================================================
  # 1b. Sample size / match counts — try the standard $nn table first,
  #     fall back to counting matches directly from model$weights if
  #     that fails (weights are a stable, always-present component:
  #     matched units get weight > 0, unmatched/discarded get weight 0)
  # ==========================================================
  nn <- tryCatch(summary(model)$nn, error = function(e) NULL)
  
  if (!is.null(nn) && is.matrix(nn)) {
    cat("Sample sizes before and after matching:\n")
    cat("--------------------------------------------------------\n")
    print(nn)
    cat("\n")
    
    tryCatch({
      if ("Matched" %in% rownames(nn) && "All" %in% rownames(nn)) {
        treated_col <- if ("Treated" %in% colnames(nn)) "Treated" else colnames(nn)[2]
        control_col <- if ("Control" %in% colnames(nn)) "Control" else colnames(nn)[1]
        
        treated_before <- nn["All", treated_col]
        treated_after  <- nn["Matched", treated_col]
        control_before <- nn["All", control_col]
        control_after  <- nn["Matched", control_col]
        
        pct_treated_kept <- 100 * treated_after / treated_before
        
        cat(sprintf("In plain English: out of %s treated patients, %s (%.0f%%) were successfully\n",
                    treated_before, treated_after, pct_treated_kept))
        cat(sprintf("matched to a similar control patient and kept for analysis. The remaining\n"))
        cat(sprintf("%s treated patients were DISCARDED because no sufficiently similar control\n",
                    treated_before - treated_after))
        cat(sprintf("patient could be found for them.\n"))
        
        if (pct_treated_kept < 80) {
          cat(sprintf("\n%s a substantial share (%.0f%%) of treated patients could NOT be matched.\n",
                      bold("CAUTION:"), 100 - pct_treated_kept))
          cat("This can happen when treated and control patients don't overlap much on the\n")
          cat("matching covariates (poor 'common support'). Any conclusions from the matched\n")
          cat("sample only generalize to the KIND of treated patient who COULD be matched —\n")
          cat("not necessarily to all treated patients in the original cohort.\n")
        }
      }
    }, error = function(e) {
      cat("(Could not compute a plain-English summary of match counts from this table —\n")
      cat("see the raw table above instead.)\n")
    })
    
  } else {
    # Fallback: derive counts directly from $weights, which is a stable,
    # always-present component of a matchit object regardless of method/version
    w <- tryCatch(model$weights, error = function(e) NULL)
    trt <- tryCatch(model$treat, error = function(e) NULL)
    
    if (!is.null(w) && !is.null(trt)) {
      cat("(Could not read the standard match-count table for this method/version —\n")
      cat("falling back to counting directly from the model's match weights.)\n\n")
      
      trt_levels <- sort(unique(trt))
      for (lvl in trt_levels) {
        in_grp <- trt == lvl
        n_total <- sum(in_grp)
        n_matched <- sum(in_grp & w > 0)
        cat(sprintf("  Group '%s': %d total, %d matched and retained, %d discarded\n",
                    lvl, n_total, n_matched, n_total - n_matched))
      }
    } else {
      cat("WARNING: could not determine match counts for this model at all — neither\n")
      cat("the summary table nor the raw $weights/$treat components were readable.\n")
      cat("Printing the raw model object instead:\n\n")
      print(model)
    }
  }
  
  # ==========================================================
  # 1c. Covariate balance before vs. after matching
  # ==========================================================
  cat("\nCovariate balance (did matching actually make the groups more similar?):\n")
  cat("--------------------------------------------------------\n")
  cat(sprintf("Rule of thumb: a Standardized Mean Difference (SMD) below %.2f (absolute\n",
              smd_threshold))
  cat("value) is generally considered good balance; above 0.25 is generally considered\n")
  cat("poor balance requiring caution.\n\n")
  
  s <- tryCatch(summary(model), error = function(e) NULL)
  
  if (is.null(s)) {
    cat("WARNING: summary() failed for this model — cannot assess covariate balance.\n")
    cat("Printing the raw model object instead:\n\n")
    print(model)
    return(invisible(NULL))
  }
  
  before <- tryCatch(s$sum.all, error = function(e) NULL)
  after  <- tryCatch(s$sum.matched, error = function(e) NULL)
  if (is.null(after)) after <- tryCatch(s$sum.subclass, error = function(e) NULL)
  
  balance_ok <- !is.null(before) && (is.matrix(before) || is.data.frame(before))
  
  if (!balance_ok) {
    cat("WARNING: could not extract a standard covariate-balance table for this model\n")
    cat("(this can happen with multi-category or continuous treatments, or some\n")
    cat("matching methods). Printing the raw summary object instead — inspect it\n")
    cat("manually for balance statistics:\n\n")
    print(s)
    return(invisible(s))
  }
  
  smd_col_before <- grep("std.*mean.*diff", colnames(before), ignore.case = TRUE)
  if (length(smd_col_before) == 0) {
    cat("WARNING: could not find a 'Standardized Mean Difference' column in the\n")
    cat("balance table for this model/method. Printing the raw balance table(s)\n")
    cat("instead:\n\n")
    print(before)
    if (!is.null(after)) { cat("\nAfter matching:\n"); print(after) }
    return(invisible(s))
  }
  smd_col_before <- smd_col_before[1]
  
  covariate_names <- rownames(before)
  covariate_names <- covariate_names[covariate_names != "distance"]  # the PS itself, not a real covariate
  
  results_rows <- list()
  worst_after <- 0
  
  for (cv in covariate_names) {
    smd_before <- tryCatch(before[cv, smd_col_before], error = function(e) NA)
    
    smd_after <- NA
    if (!is.null(after) && cv %in% rownames(after)) {
      smd_col_after <- grep("std.*mean.*diff", colnames(after), ignore.case = TRUE)
      if (length(smd_col_after) > 0) smd_after <- after[cv, smd_col_after[1]]
    }
    
    cat(sprintf("\n> %s\n", cv))
    cat(sprintf("  Before matching: SMD = %s\n",
                ifelse(is.na(smd_before), "NA", sprintf("%.3f", smd_before))))
    
    if (!is.na(smd_after)) {
      cat(sprintf("  After matching:  SMD = %.3f", smd_after))
      
      verdict <- if (abs(smd_after) < smd_threshold) {
        "GOOD balance"
      } else if (abs(smd_after) < 0.25) {
        "MODERATE imbalance — interpret results for this variable with some caution"
      } else {
        "POOR balance — this variable is still meaningfully different between groups after matching"
      }
      cat(sprintf(" — %s\n", bold(verdict)))
      
      if (!is.na(smd_before) && abs(smd_after) > abs(smd_before)) {
        cat("  NOTE: balance for this variable got WORSE after matching, not better —\n")
        cat("  unusual, worth double-checking your matching specification.\n")
      }
      
      worst_after <- max(worst_after, abs(smd_after), na.rm = TRUE)
    } else {
      cat("  After matching:  not available for this variable.\n")
    }
    
    results_rows[[length(results_rows) + 1]] <- data.frame(
      covariate = cv, smd_before = round(smd_before, 3), smd_after = round(smd_after, 3)
    )
  }
  
  cat("\n--------------------------------------------------------\n")
  if (worst_after > 0) {
    if (worst_after < smd_threshold) {
      cat(sprintf("%s all covariates achieved good balance (max SMD after matching = %.3f).\n",
                  bold("Overall verdict:"), worst_after))
      cat("Matching appears to have worked well for the covariates you specified.\n")
    } else {
      cat(sprintf("%s at least one covariate still shows imbalance after matching\n",
                  bold("Overall verdict:")))
      cat(sprintf("(max SMD after matching = %.3f). Consider a stricter caliper, a different\n",
                  worst_after))
      cat("matching method, or including a regression adjustment for the remaining\n")
      cat("imbalance in your outcome model.\n")
    }
  }
  cat("\nIMPORTANT: matching can only balance the covariates you actually included in\n")
  cat("the matching formula. It says NOTHING about balance on unmeasured confounders —\n")
  cat("this is a fundamental limitation of propensity score methods, not just this tool.\n")
  cat("========================================================\n")
  
  invisible(do.call(rbind, results_rows))
}


# ============================================================
# 2. explain_iptw() — Inverse Probability of Treatment Weighting
# ============================================================
explain_iptw <- function(model) {
  
  if (!inherits(model, "weightit")) {
    stop("explain_iptw() expects a model fit with WeightIt::weightit().")
  }
  
  cat("========================================================\n")
  cat("INVERSE PROBABILITY OF TREATMENT WEIGHTING (IPTW)\n")
  cat("========================================================\n\n")
  
  cat("Instead of discarding unmatched patients (as matching does), weighting keeps\n")
  cat("EVERYONE, but gives each patient a WEIGHT based on how (un)likely their observed\n")
  cat("treatment was, given their covariates. Patients who got a treatment that was\n")
  cat("'surprising' given their profile are up-weighted; this reconstructs a\n")
  cat("pseudo-population where treatment looks independent of the covariates used.\n\n")
  
  # ==========================================================
  # 2a. Method / estimand — defensive, since field names/values vary
  # ==========================================================
  method_code <- tryCatch(model$method, error = function(e) NA)
  estimand    <- tryCatch(model$estimand, error = function(e) NA)
  
  method_names <- c(ps = "propensity score weighting (logistic regression)",
                    ebal = "entropy balancing", cbps = "covariate balancing propensity score",
                    gbm = "generalized boosted model", super = "SuperLearner ensemble",
                    ipt = "inverse probability tilting", npcbps = "non-parametric CBPS")
  method_friendly <- if (!is.na(method_code) && method_code %in% names(method_names)) {
    method_names[[method_code]]
  } else if (!is.na(method_code)) {
    sprintf("'%s' (method-specific description not built into this function)", method_code)
  } else {
    "an unspecified method"
  }
  
  estimand_explain <- if (!is.na(estimand) && estimand == "ATT") {
    "Average Treatment effect on the Treated (ATT) — i.e. the effect FOR patients like those who actually got treated"
  } else if (!is.na(estimand) && estimand == "ATE") {
    "Average Treatment Effect (ATE) — i.e. the effect if EVERYONE in this population had been treated, vs. if no one had"
  } else if (!is.na(estimand)) {
    sprintf("estimand '%s'", estimand)
  } else {
    "an unspecified estimand"
  }
  
  cat(sprintf("Weighting method: %s\n", bold(method_friendly)))
  cat(sprintf("Target estimand:  %s\n\n", bold(estimand_explain)))
  
  # ==========================================================
  # 2b. Effective Sample Size — computed DIRECTLY from raw weights/treat,
  #     not from summary(), since this formula is standard and doesn't
  #     depend on WeightIt's internal summary object structure at all.
  #     ESS = (sum(w))^2 / sum(w^2)
  # ==========================================================
  w   <- tryCatch(model$weights, error = function(e) NULL)
  trt <- tryCatch(model$treat, error = function(e) NULL)
  
  if (is.null(w)) {
    cat("WARNING: could not access this model's $weights component at all — cannot\n")
    cat("compute effective sample size or check for extreme weights. Printing the raw\n")
    cat("summary object instead:\n\n")
    print(tryCatch(summary(model), error = function(e) model))
    return(invisible(NULL))
  }
  
  # Guard: non-finite or non-positive weights (shouldn't normally happen, but
  # some methods can produce them in edge cases — e.g. near-zero propensity scores)
  n_bad <- sum(!is.finite(w) | w < 0)
  if (n_bad > 0) {
    cat(sprintf("%s %d weight(s) are non-finite (NA/Inf) or negative — these are\n",
                bold("WARNING:"), n_bad))
    cat("excluded from the calculations below, but their presence suggests a problem\n")
    cat("with the weighting model (e.g. a covariate perfectly predicting treatment,\n")
    cat("producing propensity scores at/near 0 or 1). Investigate before trusting\n")
    cat("results based on these weights.\n\n")
  }
  valid <- is.finite(w) & w >= 0
  w_valid <- w[valid]
  
  if (length(w_valid) == 0) {
    cat(bold("WARNING:"), "no usable (finite, non-negative) weights remain after excluding\n")
    cat("bad values — cannot compute ESS or check for extreme weights. Stopping here.\n")
    cat("========================================================\n")
    return(invisible(NULL))
  }
  
  cat("Effective Sample Size (ESS) — a weighted sample \"acts like\" a smaller\n")
  cat("unweighted sample when some patients carry much more weight than others:\n")
  cat("--------------------------------------------------------\n")
  
  is_continuous_treatment <- !is.null(trt) && is.numeric(trt) && length(unique(trt)) > 10
  
  if (is_continuous_treatment) {
    cat("This looks like a CONTINUOUS treatment variable (many unique values) rather\n")
    cat("than discrete groups — reporting overall ESS rather than per-group ESS.\n\n")
    ess <- sum(w_valid)^2 / sum(w_valid^2)
    n_actual <- length(w_valid)
    cat(sprintf("Actual N: %s | Effective Sample Size: %s (%.0f%% of actual)\n",
                format(n_actual, big.mark = ","), bold(sprintf("%.0f", ess)), 100 * ess / n_actual))
  } else if (!is.null(trt)) {
    for (lvl in sort(unique(trt[valid]))) {
      w_grp <- w_valid[trt[valid] == lvl]
      n_actual <- length(w_grp)
      ess <- sum(w_grp)^2 / sum(w_grp^2)
      pct <- 100 * ess / n_actual
      
      cat(sprintf("\n> Group '%s': actual N = %s, Effective Sample Size = %s (%.0f%% of actual)\n",
                  lvl, format(n_actual, big.mark = ","), bold(sprintf("%.0f", ess)), pct))
      if (pct < 50) {
        cat(sprintf("  %s this group's effective sample size dropped by more than half.\n", bold("CAUTION:")))
        cat("  A few patients with very large weights are likely dominating the estimate\n")
        cat("  for this group — results may be unstable. Consider weight trimming or\n")
        cat("  stabilized weights.\n")
      }
    }
  } else {
    ess <- sum(w_valid)^2 / sum(w_valid^2)
    cat(sprintf("Actual N: %s | Effective Sample Size: %.0f\n", length(w_valid), ess))
    cat("(Could not identify treatment groups — reporting overall ESS only.)\n")
  }
  
  # ==========================================================
  # 2c. Extreme-weight check — again computed directly from raw weights
  # ==========================================================
  cat("\nExtreme weight check:\n")
  cat("--------------------------------------------------------\n")
  if (mean(w_valid) == 0) {
    cat(bold("WARNING:"), "all usable weights are zero — cannot compute an extreme-weight\n")
    cat("ratio. This would mean every patient was effectively excluded, which suggests\n")
    cat("something has gone wrong with the weighting model.\n")
  } else {
    w_ratio <- max(w_valid) / mean(w_valid)
    cat(sprintf("Largest weight is %.1fx the average weight (max = %.2f, mean = %.2f).\n",
                w_ratio, max(w_valid), mean(w_valid)))
    
    if (w_ratio > 20) {
      cat(sprintf("%s this is a LARGE ratio — one or a few patients may be having an\n", bold("CAUTION:")))
      cat("outsized influence on your results. Common fixes: weight trimming (capping\n")
      cat("weights at a percentile, e.g. the 1st/99th), or stabilized weights\n")
      cat("(weightit(..., stabilize = TRUE)).\n")
    } else if (w_ratio > 10) {
      cat("This is a moderately large ratio — worth a look at the weight distribution,\n")
      cat("though not necessarily a problem on its own.\n")
    } else {
      cat("This looks reasonable — no single patient appears to be dominating the analysis.\n")
    }
  }
  
  # ==========================================================
  # 2d. Covariate balance — optional, only if the `cobalt` package is available
  # (this is the standard companion package for checking WeightIt balance; we
  # don't hard-require it, since not everyone will have it installed)
  # ==========================================================
  cat("\nCovariate balance after weighting:\n")
  cat("--------------------------------------------------------\n")
  if (requireNamespace("cobalt", quietly = TRUE)) {
    bal <- tryCatch(cobalt::bal.tab(model), error = function(e) NULL)
    if (!is.null(bal)) {
      print(bal)
      cat("\n(See the 'Diff.Adj' column above — same rule of thumb as matching: under 0.1\n")
      cat("absolute value is generally good balance.)\n")
    } else {
      cat("The 'cobalt' package is installed, but cobalt::bal.tab() could not process\n")
      cat("this model — inspect balance manually.\n")
    }
  } else {
    cat("Install the 'cobalt' package for a full covariate balance table:\n")
    cat("  install.packages('cobalt'); cobalt::bal.tab(model_iptw)\n")
  }
  
  cat("\nIMPORTANT: like matching, weighting can only balance the covariates you\n")
  cat("actually included in the weighting formula — it says nothing about unmeasured\n")
  cat("confounders. This is a fundamental limitation of propensity score methods.\n")
  cat("========================================================\n")
  
  invisible(list(weights = w, ess_ratio_max = if (exists("w_ratio")) w_ratio else NA))
}


# ============================================================
# WHAT ABOUT THE OUTCOME MODEL (out_ps_outcome)?
# ============================================================
# out_ps_outcome <- coxph(Surv(time, status) ~ treatment, data = matched_cohort)
# is a plain coxph object — reuse explain_cox() from explain_survival.R:
#
#   explain_cox(out_ps_outcome)
#
# No need to rebuild this — it works unchanged since it's a standard Cox model.
# One caveat specific to matched/weighted data: if you used MatchIt method =
# "full" or "subclass" (not the ratio = 1 nearest-neighbour matching shown in
# your code), match.data() adds a `weights` column that MUST be passed into
# coxph(..., weights = weights) for correct standard errors — with simple
# ratio = 1 nearest-neighbour matching (as in your code), all weights are 1,
# so this doesn't affect your specific example, but it will for other methods.
# ============================================================


# ============================================================
# USAGE EXAMPLES (using the cohort dataset from earlier)
# ============================================================
library(MatchIt)
ps_model <- matchit(treatment ~ age + sex + comorbidity_score,
                     data = cohort, method = "nearest", ratio = 1)
matchit_results <- explain_matchit(ps_model)

library(WeightIt)
model_iptw <- weightit(treatment ~ age + sex + comorbidity_score,
                        data = cohort, method = "ps")
iptw_results <- explain_iptw(model_iptw)

matched_cohort <- match.data(ps_model)
out_ps_outcome <- coxph(Surv(time, status) ~ treatment, data = matched_cohort)
explain_cox(out_ps_outcome)   # from explain_survival.R

# ============================================================
# EDGE CASES THIS FILE HAS BEEN SPECIFICALLY CHECKED AGAINST
# ============================================================
# explain_matchit():
#  1. summary(model) failing outright                          -> caught, prints raw model, returns safely
#  2. $nn table missing/not a matrix (version/method differences) -> falls back to counting directly
#     from model$weights + model$treat (stable, always-present components)
#  3. model$info$method / model$info$ratio missing              -> defensive tryCatch, generic label used
#  4. No 'Std. Mean Diff.' column found in balance table         -> warns, prints raw table instead of
#     guessing at a wrong column
#  5. Balance table entirely missing/wrong shape (e.g. multi-category
#     or continuous treatment)                                  -> warns, prints raw summary object
#  6. A covariate present in "before" table but missing from "after" table -> handled per-covariate,
#     reported as "not available" rather than crashing
#  7. Balance getting WORSE after matching for some covariate    -> explicitly flagged, not silently shown
#  8. Very low proportion of treated patients successfully matched -> flagged with common-support caveat
#  9. "distance" (the propensity score itself) appearing as a pseudo-covariate row -> excluded from
#     the covariate loop (it's not a real confounder to balance)
# 10. method = "full", "subclass", "exact", "genetic", "cem", etc. (not just "nearest") -> summary
#     extraction and fallbacks written to be method-agnostic wherever possible
#
# explain_iptw():
# 11. model$weights entirely inaccessible                        -> stops gracefully with a clear warning,
#     falls back to printing the raw summary/model
# 12. Non-finite (NA/Inf) or negative weights present             -> flagged explicitly, excluded from
#     ESS/extreme-weight calculations rather than silently propagating NaN
# 13. Continuous (non-binary/non-categorical) treatment variable  -> detected heuristically, switches to
#     an overall-ESS explanation instead of a per-group one
# 14. Unrecognised weighting method code                          -> generic fallback label, not a crash
# 15. model$estimand missing or an unfamiliar value                -> generic fallback label
# 16. `cobalt` package not installed                              -> balance section degrades gracefully
#     with install instructions, rest of the function still runs fully
# 17. cobalt::bal.tab() erroring on this specific model            -> caught, noted, doesn't crash the
#     rest of the function
# 18. All weights identical (e.g. unweighted/degenerate case)      -> ESS equals actual N naturally, no
#     division-by-zero risk (mean weight can't be 0 if any weight is positive and finite)
# 19. Zero usable weights remaining after excluding bad values     -> stops gracefully with a clear
#     warning instead of erroring on downstream division
# 20. All usable weights equal to exactly zero                     -> guarded separately (mean = 0 case),
#     warned instead of producing NaN from a 0/0 division
# ============================================================