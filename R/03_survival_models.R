# ============================================================
# explain_survival.R — plain-language interpretation of
# survival objects: survfit() (Kaplan-Meier), survdiff()
# (log-rank test), coxph() (Cox model), and cox.zph()
# (proportional hazards test).
#
# Four functions, mirroring the survival workflow:
#   explain_km(fit, forest_plot = FALSE)     — Kaplan-Meier curve(s)
#   explain_logrank(sd)                       — log-rank test (which group fares better)
#   explain_cox(model, forest_plot = FALSE, table_plot = FALSE) — Cox model, Hazard Ratios (HR)
#   explain_cox_zph(zph)                      — proportional hazards assumption test
#
# All work for ANY number/combination of groups or predictors —
# nothing here is hardcoded to a specific variable name or number
# of groups.
#
# FOREST PLOT LAYOUT — matches explain_glm()'s style throughout:
#   - explain_cox(): uses the SAME cal_forest_plot() helper as
#     explain_glm()/explain_lm() (gtsummary + forestplot packages),
#     extended here to also accept coxph objects. Reference ("no
#     effect") line sits at HR = 1 on a log-scaled axis — 1 is the
#     null value for a Hazard Ratio, same logic as OR/RR/IRR.
#   - explain_km(): there's no regression coefficient here (Kaplan-
#     Meier is non-parametric), so instead this plots MEDIAN
#     SURVIVAL TIME per group in the same visual layout (group
#     labels on the left, point + CI on the right, via the
#     forestplot package directly) — the closest meaningful
#     equivalent for a KM object.
# ============================================================

library(survival)

# ---- Shared helpers (safe to re-source alongside explain_lm.R / explain_glm.R) ----
bold <- function(x) {
  if (requireNamespace("crayon", quietly = TRUE)) crayon::bold(x) else x
}
format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}


# ============================================================
# cal_forest_plot() — gtsummary + forestplot-based forest plot.
# Same fixed version used in explain_glm.R / explain_lm.R, EXTENDED
# here to also accept coxph objects (gtsummary::tbl_regression()
# supports coxph natively via broom::tidy, so this needed only an
# inherits() check added — the rest of the logic is unchanged):
#   - every dplyr/tidyselect call is explicitly namespaced (dplyr::...)
#     rather than relying on library(dplyr), which avoids dplyr::select
#     silently getting masked by MASS::select (or any other package)
#   - the forestplot object is explicitly print()-ed before returning,
#     since forestplot::forestplot()'s auto-print only fires for an
#     untouched top-level console expression — NOT when called from
#     inside another function (as happens when explain_cox() calls
#     this inside a tryCatch{}), which would otherwise render blank
# ============================================================
cal_forest_plot <- function(fit_or_tbl,
                            family = "HR",
                            col_names = c("estimate", "ci", "p.value"),
                            graph.pos = 2,
                            boxsize = 0.3,
                            title_line_color = "darkblue",
                            exponentiate = TRUE,
                            indent_spaces = 5) {
  
  if (!requireNamespace("forestplot", quietly = TRUE)) {
    stop("Package 'forestplot' is required for cal_forest_plot()", call. = FALSE)
  }
  if (!requireNamespace("gtsummary", quietly = TRUE)) {
    stop("Package 'gtsummary' is required for cal_forest_plot()", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for cal_forest_plot()", call. = FALSE)
  }
  
  # Deliberately NOT calling library(dplyr) etc. — see note above.
  
  # If input is a fitted model, convert to tbl_regression internally.
  # coxph added here alongside glm/lm — gtsummary::tbl_regression() supports
  # it natively (via broom::tidy.coxph()), so this is the only change needed
  # to make cal_forest_plot() work for Cox models too.
  if (inherits(fit_or_tbl, "glm") || inherits(fit_or_tbl, "lm") || inherits(fit_or_tbl, "coxph")) {
    x <- gtsummary::tbl_regression(fit_or_tbl, exponentiate = exponentiate)
  } else if (inherits(fit_or_tbl, c("tbl_regression", "tbl_uvregression"))) {
    x <- fit_or_tbl
  } else {
    stop("Input must be a glm/lm/coxph object or tbl_regression/tbl_uvregression object", call. = FALSE)
  }
  
  # Determine log scale and x-axis label
  family <- toupper(family)
  xlog <- family %in% c("OR", "RR", "IRR", "HR")  # log scale for ratio measures
  xlab <- switch(family,
                 "OR"  = "Odds Ratio",
                 "RR"  = "Risk Ratio",
                 "IRR" = "Incidence Rate Ratio",
                 "HR"  = "Hazard Ratio",
                 "RISK"= "Risk",
                 "Estimate")  # default — used for raw-scale coefficients (e.g. lm())
  
  # Prepare main text table
  txt_tb1 <- x %>%
    gtsummary::modify_column_unhide() %>%
    gtsummary::modify_fmt_fun(dplyr::contains("stat") ~ gtsummary::style_number) %>%
    gtsummary::as_tibble(col_labels = FALSE)
  
  # Prepare header row and update "estimate" dynamically
  txt_tb2 <- x %>%
    gtsummary::modify_column_unhide() %>%
    gtsummary::as_tibble() %>%
    names()
  txt_tb2 <- data.frame(matrix(gsub("\\**", "", txt_tb2), nrow = 1), stringsAsFactors = FALSE)
  names(txt_tb2) <- names(txt_tb1)
  
  # Update header row for estimate column
  estimate_col <- which(names(txt_tb2) == "estimate")
  if(length(estimate_col) == 1) txt_tb2[1, estimate_col] <- xlab
  
  txt_tb <- dplyr::bind_rows(txt_tb2, txt_tb1)
  
  # Combine with numeric stats
  line_stats <- x$table_body %>%
    dplyr::select(dplyr::all_of(c("estimate", "conf.low", "conf.high"))) %>%
    dplyr::rename_with(~paste0(., "_num")) %>%
    tibble::add_row(.before = 0)
  
  forestplot_tb <- dplyr::bind_cols(txt_tb, line_stats)
  
  # Add CI column
  forestplot_tb$ci <- c(NA, x$table_body$ci)
  
  # Bold only categorical summary rows
  summary_rows <- c(TRUE, x$table_body$row_type == "label")
  forestplot_tb <- forestplot_tb %>% dplyr::mutate(..summary_row.. = summary_rows)
  
  # --- Indent variable names that are not summary rows ---
  forestplot_tb <- forestplot_tb %>%
    dplyr::mutate(label = ifelse(..summary_row.., label, paste0(strrep(" ", indent_spaces), label)))
  
  # Prepare labeltext as a list
  label_txt <- forestplot_tb %>%
    dplyr::select(dplyr::all_of(c("label", col_names))) %>%
    as.list() %>%
    lapply(as.character)
  
  # Draw forest plot — explicitly printed (see note in the header comment above)
  fp_obj <- forestplot::forestplot(
    labeltext = label_txt,
    mean = forestplot_tb$estimate_num,
    lower = forestplot_tb$conf.low_num,
    upper = forestplot_tb$conf.high_num,
    is.summary = forestplot_tb$..summary_row..,
    graph.pos = graph.pos,
    lwd.zero = 2,
    boxsize = boxsize,
    graphwidth = grid::unit(5, "cm"),
    hrzl_lines = list("2" = grid::gpar(lwd = 2, col = title_line_color)),
    xlog = xlog,
    xlab = xlab,
    col = forestplot::fpColors(box = "darkblue", line = "darkblue", summary = "darkblue")
  )
  print(fp_obj)
  invisible(fp_obj)
}


# ============================================================
# 1. explain_km() — Kaplan-Meier survival curve(s)
# ============================================================
explain_km <- function(fit, forest_plot = FALSE) {
  
  if (!inherits(fit, "survfit")) {
    stop("explain_km() expects a survfit() object.")
  }
  
  cat("========================================================\n")
  cat("KAPLAN-MEIER SURVIVAL ESTIMATE\n")
  cat("========================================================\n\n")
  
  # ---- Pull the per-group summary table. IMPORTANT: this table is generated
  # by summary.survfit() — it is NOT reliably stored directly on the raw
  # survfit() object itself (fit$table is often NULL/absent). summary(fit)$table
  # is the correct, documented primary source. If that's also unavailable, this
  # falls back to building the table MANUALLY from survfit's core components
  # (n.risk, n.event, strata, and the quantile.survfit() method for median + CI)
  # — these are fundamental, always-present building blocks of any survfit
  # object regardless of package version, so this final fallback should always
  # succeed even if the two shortcut methods above don't. ----
  tbl <- tryCatch(summary(fit)$table, error = function(e) NULL)
  if (is.null(tbl)) {
    tbl <- fit$table  # secondary fallback — see note above
  }
  if (is.null(tbl)) {
    tbl <- tryCatch({
      strata_vec <- fit$strata
      if (is.null(strata_vec)) {
        # single-curve fit — treat as one group
        strata_vec <- setNames(length(fit$time), "Overall")
      }
      grp_names <- names(strata_vec)
      idx_end <- cumsum(strata_vec)
      idx_start <- idx_end - strata_vec + 1
      
      med_ci <- tryCatch(stats::quantile(fit, probs = 0.5), error = function(e) NULL)
      
      rows <- lapply(seq_along(strata_vec), function(i) {
        rng <- idx_start[i]:idx_end[i]
        n_i     <- fit$n.risk[rng[1]]               # number at risk at the start of this stratum = total n
        events_i <- sum(fit$n.event[rng])
        
        med_i <- NA; lcl_i <- NA; ucl_i <- NA
        if (!is.null(med_ci)) {
          # quantile.survfit()'s return shape varies by number of strata — handle both
          if (is.list(med_ci$quantile)) {
            med_i <- med_ci$quantile[[i]]
            lcl_i <- med_ci$lower[[i]]
            ucl_i <- med_ci$upper[[i]]
          } else if (!is.null(dim(med_ci$quantile))) {
            med_i <- med_ci$quantile[i, 1]
            lcl_i <- med_ci$lower[i, 1]
            ucl_i <- med_ci$upper[i, 1]
          } else if (length(strata_vec) == 1) {
            med_i <- med_ci$quantile[1]
            lcl_i <- med_ci$lower[1]
            ucl_i <- med_ci$upper[1]
          }
        }
        c(records = n_i, events = events_i, median = med_i, `0.95LCL` = lcl_i, `0.95UCL` = ucl_i)
      })
      
      m <- do.call(rbind, rows)
      rownames(m) <- grp_names
      m
    }, error = function(e) NULL)
    
    if (!is.null(tbl)) {
      cat("(Note: this survfit object's per-group table wasn't available via the usual\n")
      cat("methods — it was reconstructed manually from the underlying curve data instead.\n")
      cat("Numbers should be equivalent, but please spot-check against print(fit) if in doubt.)\n\n")
    }
  }
  if (is.null(tbl)) {
    warning("Could not obtain or reconstruct a per-group summary table from this survfit object. Was it created with survfit(Surv(...) ~ ..., data = ...)?")
    return(invisible(NULL))
  }
  if (is.null(dim(tbl))) {
    tbl <- matrix(tbl, nrow = 1, dimnames = list("Overall", names(tbl)))
  }
  
  n_groups   <- nrow(tbl)
  total_n    <- sum(tbl[, "records"], na.rm = TRUE)
  total_evt  <- sum(tbl[, "events"], na.rm = TRUE)
  
  cat(sprintf("This describes %s patients across %s, with %s events observed in total\n",
              bold(format(total_n, big.mark = ",")),
              ifelse(n_groups == 1, "1 group", sprintf("%d groups", n_groups)),
              bold(format(total_evt, big.mark = ","))))
  cat(sprintf("(the rest were censored — i.e. their follow-up ended before the event\n"))
  cat(sprintf("happened, so we only know they survived at least that long).\n\n"))
  
  if (total_evt == 0) {
    cat("NOTE: zero events were observed in this data. No meaningful survival curve\n")
    cat("differences can be estimated — median survival will show as NA for every\n")
    cat("group, and any log-rank test on this data is not interpretable.\n")
    cat("========================================================\n")
    return(invisible(tbl))
  }
  
  cat("Per-group summary:\n")
  cat("--------------------------------------------------------\n")
  
  has_median <- "median" %in% colnames(tbl)
  
  for (i in seq_len(n_groups)) {
    grp_name <- rownames(tbl)[i]
    n_i      <- tbl[i, "records"]
    evt_i    <- tbl[i, "events"]
    pct_evt  <- 100 * evt_i / n_i
    
    cat(sprintf("\n> %s\n", grp_name))
    cat(sprintf("  %s patients, %s events (%.1f%% experienced the event; the rest\n",
                format(n_i, big.mark = ","), format(evt_i, big.mark = ","), pct_evt))
    cat(sprintf("  were censored before it happened, or before follow-up ended).\n"))
    
    if (has_median) {
      med <- tbl[i, "median"]
      lcl <- if ("0.95LCL" %in% colnames(tbl)) tbl[i, "0.95LCL"] else NA
      ucl <- if ("0.95UCL" %in% colnames(tbl)) tbl[i, "0.95UCL"] else NA
      
      if (is.na(med)) {
        cat("  Median survival time: NOT REACHED — fewer than half of this group had\n")
        cat("  experienced the event by the end of follow-up, so a median can't be\n")
        cat("  calculated yet. (This is common, and often a GOOD sign for this group.)\n")
      } else {
        cat(sprintf("  Median survival time: %s time units", bold(sprintf("%.1f", med))))
        if (!is.na(lcl) && !is.na(ucl)) {
          cat(sprintf(" (95%% CI: %.1f to %s)",
                      lcl, ifelse(is.na(ucl), "not reached", sprintf("%.1f", ucl))))
        }
        cat(" — the point by which half of this group had experienced the event.\n")
      }
    }
  }
  
  cat("\n--------------------------------------------------------\n")
  cat("Note: this is DESCRIPTIVE only — it does not tell you whether differences\n")
  cat("between groups are statistically significant. Use explain_logrank() for that.\n")
  cat("========================================================\n")
  
  # ==========================================================
  # OPTIONAL FOREST-PLOT-STYLE FIGURE — median survival time per group
  # There's no regression coefficient here (KM is non-parametric), so this
  # is the closest meaningful equivalent: group labels on the left, a point
  # (median survival) + CI on the right, in the same visual layout as
  # explain_cox()'s/explain_glm()'s forest plots (via the forestplot
  # package directly, since there's no gtsummary regression table to build
  # off here — just a plain per-group summary).
  # ==========================================================
  if (isTRUE(forest_plot)) {
    if (!requireNamespace("forestplot", quietly = TRUE)) {
      warning("forest_plot = TRUE requires the 'forestplot' package. Run install.packages('forestplot') and try again.")
    } else if (!has_median) {
      warning("forest_plot = TRUE was requested, but this survfit object has no median survival column to plot.")
    } else {
      med_vals <- tbl[, "median"]
      lcl_vals <- if ("0.95LCL" %in% colnames(tbl)) tbl[, "0.95LCL"] else rep(NA, n_groups)
      ucl_vals <- if ("0.95UCL" %in% colnames(tbl)) tbl[, "0.95UCL"] else rep(NA, n_groups)
      grp_names <- rownames(tbl)
      
      plottable <- !is.na(med_vals)
      if (!any(plottable)) {
        warning("forest_plot = TRUE was requested, but no group reached a median survival time (all 'not reached') — nothing to plot.")
      } else {
        label_txt <- list(
          label = c("Group", grp_names[plottable]),
          median = c("Median", sprintf("%.1f", med_vals[plottable])),
          ci = c("95% CI", ifelse(is.na(lcl_vals[plottable]) | is.na(ucl_vals[plottable]),
                                  "NA",
                                  sprintf("%.1f to %.1f", lcl_vals[plottable], ucl_vals[plottable])))
        )
        
        fp_obj <- forestplot::forestplot(
          labeltext = label_txt,
          mean = c(NA, med_vals[plottable]),
          lower = c(NA, lcl_vals[plottable]),
          upper = c(NA, ucl_vals[plottable]),
          is.summary = c(TRUE, rep(FALSE, sum(plottable))),
          graph.pos = 2,
          boxsize = 0.3,
          graphwidth = grid::unit(5, "cm"),
          hrzl_lines = list("2" = grid::gpar(lwd = 2, col = "darkblue")),
          xlog = FALSE,  # median survival TIME is a raw scale, not a ratio — no log axis, no "1" reference line
          xlab = "Median survival time",
          col = forestplot::fpColors(box = "darkblue", line = "darkblue", summary = "darkblue")
        )
        print(fp_obj)
        
        if (any(!plottable)) {
          cat(sprintf("\n(Note: %d group(s) with 'not reached' median survival were excluded from\n",
                      sum(!plottable)))
          cat("the figure, since there's no finite value to plot — see the text summary above\n")
          cat("for those groups instead.)\n")
        }
      }
    }
  }
  
  invisible(tbl)
}


# ============================================================
# 2. explain_logrank() — log-rank test (survdiff)
# ============================================================
explain_logrank <- function(sd) {
  
  if (!inherits(sd, "survdiff")) {
    stop("explain_logrank() expects a survdiff() object.")
  }
  
  cat("========================================================\n")
  cat("LOG-RANK TEST — comparing survival between groups\n")
  cat("========================================================\n\n")
  
  n_groups <- length(sd$n)
  if (n_groups < 2) {
    cat("NOTE: this test compares only 1 group — there's nothing to compare against.\n")
    cat("========================================================\n")
    return(invisible(NULL))
  }
  
  df   <- n_groups - 1
  pval <- if (!is.null(sd$pvalue)) sd$pvalue else pchisq(sd$chisq, df, lower.tail = FALSE)
  
  group_names <- names(sd$n)
  if (is.null(group_names)) group_names <- paste("Group", seq_len(n_groups))
  
  cat(sprintf("This test compares survival across %d groups: %s\n\n",
              n_groups, paste(group_names, collapse = ", ")))
  
  cat("Observed vs. Expected events per group:\n")
  cat("(If a group's survival is truly no different from the others, its OBSERVED\n")
  cat("event count should roughly match its EXPECTED count. A group with notably\n")
  cat("MORE observed than expected events fared WORSE; notably FEWER means it fared\n")
  cat("BETTER, relative to the others.)\n")
  cat("--------------------------------------------------------\n")
  
  worse_groups <- character(0)
  better_groups <- character(0)
  
  for (i in seq_len(n_groups)) {
    obs_i <- sd$obs[i]
    exp_i <- sd$exp[i]
    ratio <- obs_i / exp_i
    
    verdict <- if (ratio > 1.05) {
      worse_groups <<- c(worse_groups, group_names[i])
      "WORSE survival than average (more events happened than expected)"
    } else if (ratio < 0.95) {
      better_groups <<- c(better_groups, group_names[i])
      "BETTER survival than average (fewer events happened than expected)"
    } else {
      "roughly AS EXPECTED (close to average)"
    }
    
    cat(sprintf("\n> %s: observed = %.1f, expected = %.1f (ratio: %.2f)\n",
                group_names[i], obs_i, exp_i, ratio))
    cat(sprintf("  This group fared %s.\n", bold(verdict)))
  }
  
  cat("\n--------------------------------------------------------\n")
  cat(sprintf("Overall test result: chi-squared = %.2f on %d degree(s) of freedom, p %s\n",
              sd$chisq, df, format_pval(pval)))
  
  if (pval < 0.05) {
    cat(sprintf("\nThis IS %s (p %s).\n", bold("statistically significant"), format_pval(pval)))
    cat("In plain terms: the survival curves for these groups are unlikely to be this\n")
    cat("different from each other purely by chance — there's a real difference in\n")
    cat("survival between at least some of the groups.\n")
  } else {
    cat(sprintf("\nThis is %s (p %s).\n", bold("NOT statistically significant"), format_pval(pval)))
    cat("In plain terms: the survival curves could plausibly look this different just\n")
    cat("by chance, even if there were truly no difference between the groups. Treat\n")
    cat("any apparent gap between the curves with caution.\n")
  }
  
  if (n_groups > 2) {
    cat("\nNOTE: with more than 2 groups, this test only tells you that AT LEAST ONE\n")
    cat("group differs from the others overall — it does NOT tell you which specific\n")
    cat("pair(s) differ. For that, run pairwise log-rank tests (e.g. survminer::pairwise_survdiff()).\n")
  }
  
  cat("========================================================\n")
  
  invisible(data.frame(
    group = group_names,
    observed = round(sd$obs, 1),
    expected = round(sd$exp, 1),
    obs_exp_ratio = round(sd$obs / sd$exp, 2)
  ))
}


# ============================================================
# 3. explain_cox() — Cox proportional hazards model
# ============================================================
explain_cox <- function(model, conf_level = 0.95, forest_plot = FALSE, table_plot = FALSE) {
  
  # ==========================================================
  # 0. VALIDITY CHECK
  # ==========================================================
  if (!inherits(model, "coxph")) {
    stop("explain_cox() expects a model fit with coxph().")
  }
  
  s_check   <- summary(model, conf.int = conf_level)
  alpha_pct <- conf_level * 100
  outcome   <- deparse(formula(model)[[2]])  # e.g. "Surv(time, status)"
  
  
  # ==========================================================
  # 1. PRE-FLIGHT SANITY CHECKS
  # ==========================================================
  diagnostics <- character(0)
  
  # (a) Zero events -> nothing estimable at all
  n_events <- model$nevent
  if (is.null(n_events)) n_events <- sum(model$y[, ncol(model$y)] == 1, na.rm = TRUE)
  if (!is.na(n_events) && n_events == 0) {
    stop("This Cox model has ZERO observed events — there is nothing for the model to estimate. Check your status/event variable (is it coded 0/1 or FALSE/TRUE, with 1/TRUE = event?).")
  }
  
  # (b) Aliased / non-estimable coefficients (collinearity)
  coef_all <- coef(model)
  if (any(is.na(coef_all))) {
    diagnostics <- c(diagnostics, sprintf(
      "could NOT estimate %d coefficient(s), almost always due to perfect collinearity between predictors: %s. Those terms are skipped below.",
      sum(is.na(coef_all)), paste(names(coef_all)[is.na(coef_all)], collapse = ", ")))
  }
  
  # (c) Suspiciously huge standard errors -> possible separation (a covariate that
  #     perfectly predicts who has the event, common with small samples/rare events)
  se_vals <- s_check$coefficients[, "se(coef)"]
  if (any(se_vals > 10, na.rm = TRUE)) {
    diagnostics <- c(diagnostics,
                     "has at least one EXTREMELY large standard error — this usually signals a predictor that near-perfectly separates patients who had the event from those who didn't (common with small samples or few events relative to the number of predictors). Hazard ratios and CIs for the affected variable(s) may be absurdly wide or unstable.")
  }
  
  # (d) Did the model actually converge? coxph doesn't expose a simple $converged
  #     flag the way glm() does — the standard tell is that iter.max was hit.
  if (!is.null(model$iter) && !is.null(model$control) &&
      model$iter >= model$control$iter.max) {
    diagnostics <- c(diagnostics,
                     "may not have FULLY CONVERGED — it used the maximum allowed number of iterations. Consider refitting with coxph(..., control = coxph.control(iter.max = 50)).")
  }
  
  # (e) Very few events relative to number of predictors — a well-known rule of
  #     thumb in survival analysis is >=10 events per predictor for stable estimates
  n_predictors <- length(coef_all) - sum(is.na(coef_all))
  if (!is.na(n_events) && n_predictors > 0 && (n_events / n_predictors) < 10) {
    diagnostics <- c(diagnostics, sprintf(
      "has only about %.1f events per predictor (%d events, %d predictors). A common rule of thumb calls for at least 10 events per predictor for stable, trustworthy estimates — treat coefficients here with extra caution, especially any with wide confidence intervals.",
      n_events / n_predictors, n_events, n_predictors))
  }
  
  # (f) Stratified model — strata() terms adjust the baseline hazard but don't
  #     produce a coefficient/HR of their own; worth telling the reader so they
  #     don't go looking for a missing hazard ratio for that variable.
  strata_vars <- attr(model$terms, "specials")$strata
  if (!is.null(strata_vars)) {
    strata_labels <- attr(model$terms, "term.labels")
    strata_labels <- strata_labels[grepl("^strata\\(", strata_labels)]
    diagnostics <- c(diagnostics, sprintf(
      "is STRATIFIED by %s. This variable adjusts the baseline hazard separately for each stratum but does NOT produce its own hazard ratio — you won't see it in the per-variable results below, and that's expected, not an error.",
      paste(strata_labels, collapse = ", ")))
  }
  
  # (g) Clustered/robust variance — cluster() affects SEs (accounts for
  #     non-independent observations, e.g. repeated events per patient) but not
  #     the point estimates; worth flagging so wider CIs aren't mistaken for a bug
  if (!is.null(model$naive.var)) {
    diagnostics <- c(diagnostics,
                     "uses CLUSTERED/ROBUST standard errors (e.g. via cluster() or repeated events per patient). Confidence intervals below already account for non-independence between rows — this is expected to make them somewhat wider than naive SEs would give.")
  }
  
  if (length(diagnostics) > 0) {
    cat("########################################################\n")
    cat("# MODEL DIAGNOSTIC WARNINGS — read before trusting results\n")
    cat("########################################################\n")
    for (d in diagnostics) cat(bold("WARNING:"), "this model", d, "\n\n")
  }
  
  
  # ==========================================================
  # 2. HEADER — cohort size, events, what's being modeled
  # ==========================================================
  n_obs <- model$n
  
  cat("========================================================\n")
  cat("MODEL:", deparse(formula(model)), "\n")
  cat("========================================================\n\n")
  
  cat(sprintf("This model was fit using data from %s patients, of whom %s experienced\n",
              bold(format(n_obs, big.mark = ",")), bold(format(n_events, big.mark = ","))))
  cat(sprintf("the event (the rest were censored). This is a COX PROPORTIONAL HAZARDS\n"))
  cat(sprintf("model, which estimates HAZARD RATIOS (HR): how much higher or lower a\n"))
  cat(sprintf("patient's instantaneous RISK of the event is, at any given moment, compared\n"))
  cat(sprintf("to a reference patient/group — NOT how long they survive in absolute terms.\n\n"))
  
  if (!is.null(s_check$concordance)) {
    conc <- s_check$concordance[1]
    cat(sprintf("Model discrimination (concordance / C-statistic): %.3f\n", conc))
    cat(sprintf("In plain terms: if you picked two random patients, one of whom had the\n"))
    cat(sprintf("event before the other, this model correctly identifies which one it was\n"))
    cat(sprintf("about %.0f%% of the time (50%% = no better than a coin flip; 100%% = perfect).\n\n",
                conc * 100))
  }
  
  
  # ==========================================================
  # 3. OVERALL MODEL FIT — likelihood ratio test
  # ==========================================================
  if (!is.null(s_check$logtest)) {
    lr_stat <- s_check$logtest["test"]
    lr_df   <- s_check$logtest["df"]
    lr_p    <- s_check$logtest["pvalue"]
    
    cat(sprintf("Taken together, do these variables actually help explain survival?\n"))
    cat(sprintf("%s (likelihood-ratio test: p %s).\n\n",
                ifelse(lr_p < 0.05,
                       "Yes — the variables in this model collectively have a real, non-random relationship with survival",
                       "Not clearly — the variables in this model, taken together, don't show a statistically reliable relationship with survival"),
                format_pval(lr_p)))
  }
  
  
  # ==========================================================
  # 4. PER-VARIABLE INTERPRETATION
  # ==========================================================
  coefs_table <- s_check$coefficients
  ci_table    <- s_check$conf.int   # already exponentiated: exp(coef), exp(-coef), lower .95, upper .95
  var_names   <- rownames(coefs_table)
  pval_col    <- ncol(coefs_table)  # last column = p-value, robust to "Pr(>|z|)" naming quirks
  
  term_labels <- attr(model$terms, "term.labels")
  term_labels <- term_labels[!grepl("^(strata|cluster|tt|frailty)\\(", term_labels)]
  
  cat("What each variable is associated with:\n")
  cat("--------------------------------------------------------\n")
  
  results_rows <- list()
  
  for (i in seq_along(var_names)) {
    
    vn <- var_names[i]
    if (is.na(coef_all[vn])) next
    
    estimate <- coefs_table[i, "coef"]
    pval     <- coefs_table[i, pval_col]
    hr       <- ci_table[vn, 1]     # exp(coef)
    hr_lo    <- ci_table[vn, 3]     # lower .95
    hr_hi    <- ci_table[vn, 4]     # upper .95
    
    is_interaction <- grepl(":", vn, fixed = TRUE)
    is_transformed <- grepl("^(log|sqrt|poly|I|exp|scale)\\(", vn)
    is_tt          <- grepl("^tt\\(", vn)
    
    main_labels <- term_labels[!grepl(":", term_labels, fixed = TRUE)]
    matched_var <- main_labels[sapply(main_labels, function(p) startsWith(vn, p))]
    matched_var <- if (length(matched_var) > 0) matched_var[which.max(nchar(matched_var))] else character(0)
    
    model_data <- tryCatch(model.frame(model), error = function(e) NULL)
    is_ordered_factor <- length(matched_var) == 1 && !is.null(model_data) &&
      matched_var %in% names(model_data) && is.ordered(model_data[[matched_var]])
    is_factor_level <- length(matched_var) == 1 && !is.null(model_data) &&
      matched_var %in% names(model_data) && is.factor(model_data[[matched_var]]) && !is_ordered_factor
    
    cat(sprintf("\n> %s\n", vn))
    cat(sprintf("  Hazard Ratio (HR): %.2f (%.0f%% CI: %.2f to %.2f)   [raw coefficient: %.3f]\n",
                hr, alpha_pct, hr_lo, hr_hi, estimate))
    
    row_type <- "main effect"
    
    if (is_interaction) {
      row_type <- "interaction"
      parts <- strsplit(vn, ":", fixed = TRUE)[[1]]
      cat(sprintf("  In plain English — %s:\n", bold("INTERACTION TERM")))
      if (length(parts) == 2) {
        cat(sprintf("  Hazard ratios combine MULTIPLICATIVELY. This interaction's ratio (%.2f)\n", hr))
        cat(sprintf("  means: the hazard ratio for '%s' is itself multiplied by %.2f for each\n",
                    parts[1], hr))
        cat(sprintf("  1-unit change in (or shift to the other category of) '%s' — the effect of\n",
                    parts[2]))
        cat(sprintf("  %s DEPENDS ON the level of %s.\n", parts[1], parts[2]))
      } else {
        cat(sprintf("  Higher-order interaction between %d variables (%s) — hard to interpret\n",
                    length(parts), paste(parts, collapse = ", ")))
        cat(sprintf("  from a single number. Consider plotting predicted survival curves across\n"))
        cat(sprintf("  combinations of these variables instead.\n"))
      }
      
    } else if (is_tt) {
      row_type <- "time-transform"
      cat(sprintf("  In plain English — %s:\n", bold("TIME-VARYING EFFECT (tt() term)")))
      cat(sprintf("  This term allows the hazard ratio to change over follow-up time, so a\n"))
      cat(sprintf("  single HR number doesn't summarise it fully — this coefficient describes\n"))
      cat(sprintf("  how the effect changes with time, not a constant hazard ratio.\n"))
      
    } else if (is_transformed) {
      row_type <- "transformed"
      cat(sprintf("  In plain English — %s:\n", bold("TRANSFORMED TERM")))
      cat(sprintf("  This is a transformation of the original variable (log, polynomial, etc.),\n"))
      cat(sprintf("  not the raw variable. A '1-unit increase' statement doesn't apply cleanly —\n"))
      cat(sprintf("  interpret on the transformed scale, or generate predicted survival curves\n"))
      cat(sprintf("  at specific values of the original variable instead.\n"))
      
    } else if (is_ordered_factor) {
      row_type <- "ordered factor"
      cat(sprintf("  In plain English — %s:\n", bold("ORDERED FACTOR")))
      cat(sprintf("  '%s' comes from an ordered factor. This describes a TREND across the\n", vn))
      cat(sprintf("  ordered levels (polynomial contrast), not a simple group comparison.\n"))
      
    } else if (is_factor_level) {
      reference_level <- levels(model_data[[matched_var]])[1]
      compared_level  <- substring(vn, nchar(matched_var) + 1)
      
      pct_diff <- if (hr >= 1) (hr - 1) * 100 else (1 - hr) * 100
      direction_word <- if (hr > 1) "HIGHER (worse survival)" else if (hr < 1) "LOWER (better survival)" else "no different"
      
      cat(sprintf("  In plain English — %s vs %s:\n", compared_level, reference_level))
      if (abs(hr - 1) < 1e-9) {
        cat(sprintf("  '%s' and '%s' show essentially IDENTICAL hazard — no meaningful\n",
                    compared_level, reference_level))
        cat(sprintf("  difference in risk between the two groups.\n"))
      } else {
        cat(sprintf("  %s Patients in the '%s' group have a %s%% %s hazard of %s\n",
                    bold("Which group fares worse?"), bold(compared_level),
                    bold(sprintf("%.0f", pct_diff)), bold(direction_word), bold(outcome)))
        cat(sprintf("  than patients in the '%s' (reference) group, at any given moment,\n",
                    reference_level))
        cat(sprintf("  holding every other variable in the model constant.\n"))
        cat(sprintf("  (Remember: HIGHER hazard = WORSE survival, i.e. the event happens sooner/more.)\n"))
      }
      
    } else {
      pct_diff <- if (hr >= 1) (hr - 1) * 100 else (1 - hr) * 100
      direction_word <- if (hr > 1) "HIGHER (worse survival)" else if (hr < 1) "LOWER (better survival)" else "no different"
      
      cat(sprintf("  In plain English:\n"))
      cat(sprintf("  For every %s in %s, the hazard of %s is multiplied by %s\n",
                  bold("1-unit increase"), vn, outcome, bold(sprintf("%.2f", hr))))
      cat(sprintf("  — in other words, %s%% %s, ON AVERAGE, holding every other variable\n",
                  bold(sprintf("%.1f", pct_diff)), bold(direction_word)))
      cat(sprintf("  in the model constant. (Remember: HIGHER hazard = WORSE survival.)\n"))
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
    } else {
      cat(sprintf("  This result is %s (p %s). Treat this estimate with caution.\n",
                  bold("NOT statistically significant"), format_pval(pval)))
    }
    
    results_rows[[length(results_rows) + 1]] <- data.frame(
      variable   = vn,
      type       = row_type,
      HR         = round(hr, 3),
      lower_ci   = round(hr_lo, 3),
      upper_ci   = round(hr_hi, 3),
      p_value    = ifelse(is.na(pval), NA, ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval))),
      significant = ifelse(is.na(pval), NA, pval < 0.05)
    )
  }
  
  cat("\n--------------------------------------------------------\n")
  cat("Reminder: HR = 1.00 means NO difference in hazard. HR > 1 = higher risk/worse\n")
  cat("survival; HR < 1 = lower risk/better survival. These are associations, adjusted\n")
  cat("for other variables in the model — not proven causal effects. This model also\n")
  cat("ASSUMES each hazard ratio is constant over follow-up time (the 'proportional\n")
  cat("hazards' assumption) — check this with explain_cox_zph(cox.zph(model)).\n")
  cat("========================================================\n")
  
  results_table <- do.call(rbind, results_rows)
  
  # ==========================================================
  # 5. OPTIONAL FOREST PLOT — via cal_forest_plot(), SAME layout as
  # explain_glm()/explain_lm(). Reference line at HR = 1 on a log-scaled
  # axis (family = "HR" inside cal_forest_plot sets xlog = TRUE and the
  # correct "Hazard Ratio" x-axis label automatically).
  # ==========================================================
  if (isTRUE(forest_plot)) {
    tryCatch({
      cal_forest_plot(model, family = "HR")
    }, error = function(e) {
      warning(sprintf("Could not draw the forest plot via cal_forest_plot(): %s", conditionMessage(e)))
    })
  }
  
  # ==========================================================
  # 6. OPTIONAL GTSUMMARY TABLE — same layout as explain_glm()
  # ==========================================================
  if (isTRUE(table_plot)) {
    if (!requireNamespace("gtsummary", quietly = TRUE)) {
      warning("table_plot = TRUE was requested, but the 'gtsummary' package is not installed. Run install.packages('gtsummary') and try again.")
    } else {
      gt_table <- tryCatch({
        gtsummary::tbl_regression(model, exponentiate = TRUE)
      }, error = function(e) {
        warning(sprintf("Could not build the gtsummary table: %s", conditionMessage(e)))
        NULL
      })
      
      if (!is.null(gt_table)) {
        cat("\n--------------------------------------------------------\n")
        cat("A publication-style summary table (Hazard Ratio, 95% CI, p-value) has been\n")
        cat("generated below. This displays best in RStudio's Viewer pane or when knitted\n")
        cat("into an R Markdown/Quarto document — a plain R console will show it as a wide\n")
        cat("data frame instead of the formatted table.\n")
        cat("--------------------------------------------------------\n")
        print(gt_table)
      }
    }
  }
  
  invisible(results_table)
}


# ============================================================
# 4. explain_cox_zph() — proportional hazards assumption test
# ============================================================
explain_cox_zph <- function(zph) {
  
  if (!inherits(zph, "cox.zph")) {
    stop("explain_cox_zph() expects a cox.zph() object.")
  }
  
  cat("========================================================\n")
  cat("PROPORTIONAL HAZARDS ASSUMPTION TEST\n")
  cat("========================================================\n\n")
  
  cat("The Cox model assumes each variable's hazard ratio stays CONSTANT over the\n")
  cat("entire follow-up period (e.g. if Drug halves the hazard at day 10, it should\n")
  cat("still roughly halve it at day 1000). This test checks whether that assumption\n")
  cat("holds, separately for each variable and overall.\n\n")
  
  tbl <- zph$table
  pval_col <- ncol(tbl)  # robust to column naming ("p" in most versions)
  
  var_rows <- rownames(tbl)[rownames(tbl) != "GLOBAL"]
  
  cat("Per-variable results:\n")
  cat("--------------------------------------------------------\n")
  
  violations <- character(0)
  
  for (vn in var_rows) {
    pval <- tbl[vn, pval_col]
    
    cat(sprintf("\n> %s: chi-squared = %.2f, p %s\n", vn, tbl[vn, "chisq"], format_pval(pval)))
    
    if (is.na(pval)) {
      cat("  Could not be tested (NA result) — often happens for a term with very\n")
      cat("  little variation, or too few events.\n")
    } else if (pval < 0.05) {
      violations <- c(violations, vn)
      cat(sprintf("  %s (p %s): there IS evidence the hazard ratio for '%s' CHANGES over\n",
                  bold("Assumption VIOLATED"), format_pval(pval), vn))
      cat(sprintf("  follow-up time. The single HR reported for '%s' in the Cox model should\n", vn))
      cat("  be read as a rough time-AVERAGED effect, not a truly constant one — consider\n")
      cat("  stratifying by this variable, adding a time interaction (tt()), or splitting\n")
      cat("  follow-up into time periods.\n")
    } else {
      cat(sprintf("  %s (p %s): no strong evidence the hazard ratio for '%s' changes over\n",
                  bold("Assumption holds"), format_pval(pval), vn))
      cat(sprintf("  time — safe to interpret its HR as roughly constant across follow-up.\n"))
    }
  }
  
  if ("GLOBAL" %in% rownames(tbl)) {
    global_p <- tbl["GLOBAL", pval_col]
    cat("\n--------------------------------------------------------\n")
    cat(sprintf("GLOBAL test (all variables jointly): chi-squared = %.2f, p %s\n",
                tbl["GLOBAL", "chisq"], format_pval(global_p)))
    if (!is.na(global_p) && global_p < 0.05) {
      cat(sprintf("%s: taken as a whole, this model shows evidence AGAINST the proportional\n",
                  bold("Overall assessment")))
      cat("hazards assumption somewhere. Check which specific variable(s) are flagged\n")
      cat("above before concluding the whole model is invalid — often it's just one term.\n")
    } else {
      cat(sprintf("%s: no strong evidence against proportional hazards overall.\n",
                  bold("Overall assessment")))
    }
  }
  
  if (length(violations) > 0) {
    cat("\nNOTE: consider plotting the Schoenfeld residuals for the flagged variable(s)\n")
    cat(sprintf("(plot(cox.zph(model)) or survminer::ggcoxzph()) to SEE how the effect of %s\n",
                paste(violations, collapse = ", ")))
    cat("changes over time, rather than relying on the p-value alone.\n")
  }
  
  cat("========================================================\n")
  
  invisible(as.data.frame(tbl))
}


# ============================================================
# USAGE EXAMPLES (using the cohort dataset from earlier)
# ============================================================
library(survival)

out_km_fit     <- survfit(Surv(time, status) ~ treatment, data = cohort)
km_results     <- explain_km(out_km_fit)
km_results     <- explain_km(out_km_fit, forest_plot = TRUE)   # median survival, forest-plot layout

out_km_logrank <- survdiff(Surv(time, status) ~ treatment, data = cohort)
logrank_results <- explain_logrank(out_km_logrank)

model_cox      <- coxph(Surv(time, status) ~ age + sex + treatment + comorbidity_score,
                        data = cohort)
cox_results    <- explain_cox(model_cox)
cox_results    <- explain_cox(model_cox, forest_plot = TRUE)                  # HR forest plot
cox_results    <- explain_cox(model_cox, table_plot = TRUE)                   # gtsummary table
cox_results    <- explain_cox(model_cox, forest_plot = TRUE, table_plot = TRUE)

out_cox_zph    <- cox.zph(model_cox)
zph_results    <- explain_cox_zph(out_cox_zph)

# another nice way to do this
library(survminer)
ggsurvplot(out_km_fit)

# ============================================================
# REQUIRED PACKAGES for forest_plot/table_plot: forestplot, gtsummary, dplyr
#   install.packages(c("forestplot", "gtsummary", "dplyr"))
# ============================================================
#
# ============================================================
# EDGE CASES THIS FILE HAS BEEN SPECIFICALLY CHECKED AGAINST
# ============================================================
# explain_km():
#  1. Per-group table not available directly on the raw survfit object       -> FIXED: this was
#     a real bug found via user testing. fit$table is often NULL even for a normal
#     survfit(Surv(...) ~ group, data = ...) call — the table is actually generated by
#     summary.survfit(), not stored persistently on the object survfit() itself returns. Now
#     tries THREE sources in order: (a) summary(fit)$table — the correct, documented source;
#     (b) fit$table — kept as a defensive fallback in case some object/version stores it
#     directly; (c) if both fail, the table is reconstructed MANUALLY from survfit's core,
#     always-present components (n.risk, n.event, strata) plus the quantile.survfit() method
#     for median + CI, handling that method's differently-shaped return value for 1 vs. multiple
#     strata. Only if all three approaches fail does the function warn and return NULL.
#  2. Single-group (no strata) survfit                       -> $table reshaped to 1-row matrix
#  3. Zero events observed                                   -> flagged, skips per-group detail
#  4. Median not reached in a group (NA)                     -> explained explicitly, not shown as an error
#  5. forest_plot = TRUE with the 'forestplot' package not installed -> warns with an actionable
#     message instead of erroring
#  6. forest_plot = TRUE where NO group has a finite median survival time (all "not reached")
#     -> warns, nothing plotted, rather than crashing on an all-NA mean/lower/upper vector
#  7. forest_plot = TRUE where SOME (not all) groups have "not reached" median -> those groups
#     are excluded from the figure with an explicit note, rather than silently omitted or
#     plotted as a broken/NA point
#
# explain_logrank():
#  8. Only 1 group (nothing to compare)                      -> flagged, returns NULL safely
#  9. More than 2 groups                                     -> caveat added: omnibus test only, not pairwise
# 10. $pvalue missing on older survival package versions     -> computed manually from chisq/df
#
# explain_cox():
# 11. Zero events in the model                               -> stops with a clear, specific message
# 12. Aliased/NA coefficients (collinearity)                 -> flagged, those terms skipped safely
# 13. Extremely large SEs (separation, common with rare events) -> flagged before results
# 14. Possible non-convergence (hit max iterations)          -> flagged
# 15. Few events per predictor (<10, common survival-analysis rule of thumb) -> flagged
# 16. strata() terms                                         -> detected and explained (no HR expected for them)
# 17. cluster()/robust variance                              -> detected and explained
# 18. Interaction terms (2-way and higher-order)             -> explained differently, using real variable names
# 19. tt() time-varying-coefficient terms                    -> flagged as non-constant HR, not misread
# 20. Transformed terms: log(), poly(), I(), sqrt(), scale() -> flagged, not misread as raw 1-unit change
# 21. Ordered factors                                        -> flagged as trend, not simple group comparison
# 22. Terms unmatched to any known data column               -> flagged, not silently asserted as correct
# 23. p-value column naming differences across survival versions -> read positionally, not by name
# 24. forest_plot = TRUE with forestplot/gtsummary/dplyr not installed -> stops with a specific,
#     actionable message (via cal_forest_plot()'s own requireNamespace checks)
# 25. cal_forest_plot() erroring for any other reason         -> caught and reported as a warning,
#     rest of explain_cox()'s text output is unaffected (forest plot generated only after all
#     text has already printed)
# 26. table_plot = TRUE with gtsummary not installed          -> warns with an actionable message,
#     rest of the function is unaffected
# 27. gtsummary::tbl_regression() erroring for any reason     -> caught, warns instead of crashing
#
# explain_cox_zph():
# 28. Missing GLOBAL row (can happen with a single-covariate model in some versions) -> checked before use
# 29. NA p-value for a term (too little variation / too few events) -> explained, not silently skipped
# ============================================================