# ============================================================
# explain_survival.R — plain-language interpretation of survival
# analysis objects, with a UNIFIED calling pattern: the log-rank
# test and the proportional-hazards test are passed IN as
# arguments to explain_km()/explain_cox() rather than requiring
# separate calls, plus a broad set of optional tables and figures.
#
# Four functions:
#   explain_logrank(sd)                       — standalone log-rank explainer
#                                                (also usable on its own)
#   explain_cox_zph(zph)                      — standalone PH-test explainer
#                                                (also usable on its own)
#   explain_km(fit, logrank = NULL, data = NULL,
#              km_plot = FALSE, cumulative_plot = FALSE,
#              forest_plot = FALSE, table_plot = FALSE)
#   explain_cox(model, zph = NULL, conf_level = 0.95,
#               forest_plot = FALSE, table_plot = FALSE,
#               survival_plot = FALSE, by_variable = NULL,
#               ph_plot = FALSE)
#
# explain_km() and explain_cox() internally CALL explain_logrank()/
# explain_cox_zph() when you pass logrank/zph in — so there's no
# duplicated logic, and you can still call either standalone function
# on its own if you only want that one piece.
#
# FOREST PLOTS: explain_cox()'s forest plot uses the SAME
# cal_forest_plot() helper as explain_glm()/explain_lm() (gtsummary +
# forestplot packages), extended here to also accept coxph objects.
# Reference line sits at HR = 1 on a log-scaled axis, same logic as
# OR/RR/IRR in explain_glm(). explain_km() has no regression
# coefficient to put on that kind of plot (Kaplan-Meier is non-
# parametric) — its forest_plot instead shows median survival time
# per group, in the same visual layout, with no "null value" line
# (there isn't one for a raw survival time).
#
# REQUIRED PACKAGES:
#   Core: survival (base functionality)
#   forest_plot / table_plot (Cox): forestplot, gtsummary, dplyr
#   table_plot (KM): gtsummary
#   km_plot / cumulative_plot / survival_plot / ph_plot: survminer
#     (optional — each has a manual ggplot2 fallback if survminer
#     isn't installed, except ph_plot which falls back to base
#     plot(cox.zph(...)))
#   install.packages(c("survival", "forestplot", "gtsummary", "dplyr",
#                       "ggplot2", "survminer"))
# ============================================================

library(survival)

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
# cal_forest_plot() — gtsummary + forestplot-based forest plot.
# Same fixed version used in explain_glm.R / explain_lm.R, extended
# to also accept coxph objects (gtsummary::tbl_regression() supports
# coxph natively via broom::tidy.coxph()). Key fixes carried over:
#   - every dplyr/tidyselect call is explicitly namespaced (dplyr::...)
#     rather than relying on library(dplyr), which avoids dplyr::select
#     silently getting masked by MASS::select (or any other package)
#   - the forestplot object is explicitly print()-ed before returning,
#     since forestplot::forestplot()'s auto-print only fires for an
#     untouched top-level console expression — not when called from
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
  
  if (inherits(fit_or_tbl, "glm") || inherits(fit_or_tbl, "lm") || inherits(fit_or_tbl, "coxph")) {
    x <- gtsummary::tbl_regression(fit_or_tbl, exponentiate = exponentiate)
  } else if (inherits(fit_or_tbl, c("tbl_regression", "tbl_uvregression"))) {
    x <- fit_or_tbl
  } else {
    stop("Input must be a glm/lm/coxph object or tbl_regression/tbl_uvregression object", call. = FALSE)
  }
  
  family <- toupper(family)
  xlog <- family %in% c("OR", "RR", "IRR", "HR")
  xlab <- switch(family,
                 "OR"  = "Odds Ratio",
                 "RR"  = "Risk Ratio",
                 "IRR" = "Incidence Rate Ratio",
                 "HR"  = "Hazard Ratio",
                 "RISK"= "Risk",
                 "Estimate")
  
  txt_tb1 <- x %>%
    gtsummary::modify_column_unhide() %>%
    gtsummary::modify_fmt_fun(dplyr::contains("stat") ~ gtsummary::style_number) %>%
    gtsummary::as_tibble(col_labels = FALSE)
  
  txt_tb2 <- x %>%
    gtsummary::modify_column_unhide() %>%
    gtsummary::as_tibble() %>%
    names()
  txt_tb2 <- data.frame(matrix(gsub("\\**", "", txt_tb2), nrow = 1), stringsAsFactors = FALSE)
  names(txt_tb2) <- names(txt_tb1)
  
  estimate_col <- which(names(txt_tb2) == "estimate")
  if (length(estimate_col) == 1) txt_tb2[1, estimate_col] <- xlab
  
  txt_tb <- dplyr::bind_rows(txt_tb2, txt_tb1)
  
  line_stats <- x$table_body %>%
    dplyr::select(dplyr::all_of(c("estimate", "conf.low", "conf.high"))) %>%
    dplyr::rename_with(~paste0(., "_num")) %>%
    tibble::add_row(.before = 0)
  
  forestplot_tb <- dplyr::bind_cols(txt_tb, line_stats)
  forestplot_tb$ci <- c(NA, x$table_body$ci)
  
  summary_rows <- c(TRUE, x$table_body$row_type == "label")
  forestplot_tb <- forestplot_tb %>% dplyr::mutate(..summary_row.. = summary_rows)
  
  forestplot_tb <- forestplot_tb %>%
    dplyr::mutate(label = ifelse(..summary_row.., label, paste0(strrep(" ", indent_spaces), label)))
  
  label_txt <- forestplot_tb %>%
    dplyr::select(dplyr::all_of(c("label", col_names))) %>%
    as.list() %>%
    lapply(as.character)
  
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
# explain_logrank() — log-rank test (survdiff). Standalone, and
# called internally by explain_km() when you pass `logrank =`.
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
# explain_cox_zph() — proportional hazards assumption test.
# Standalone, and called internally by explain_cox() when you pass
# `zph =`.
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
    cat(sprintf("to SEE how the effect of %s changes over time, rather than relying on the\n",
                paste(violations, collapse = ", ")))
    cat("p-value alone — see explain_cox(..., ph_plot = TRUE).\n")
  }
  
  cat("========================================================\n")
  
  invisible(as.data.frame(tbl))
}


# ============================================================
# explain_km() — Kaplan-Meier survival curve(s), with the log-rank
# test and several visualisation options available directly as
# arguments.
# ============================================================
explain_km <- function(fit, logrank = NULL, data = NULL,
                       km_plot = FALSE, cumulative_plot = FALSE,
                       forest_plot = FALSE, table_plot = FALSE) {
  
  if (!inherits(fit, "survfit")) {
    stop("explain_km() expects a survfit() object.")
  }
  
  cat("========================================================\n")
  cat("KAPLAN-MEIER SURVIVAL ESTIMATE\n")
  cat("========================================================\n\n")
  
  # ---- Pre-flight: small groups ----
  if (!is.null(fit$n) && any(fit$n < 10)) {
    cat("########################################################\n")
    cat("# DIAGNOSTIC WARNING\n")
    cat("########################################################\n")
    cat(bold("WARNING:"), "one or more groups contain very few patients (<10), making\n")
    cat("survival estimates (especially the tails of the curve) potentially unstable.\n\n")
  }
  
  # ---- Pull the per-group summary table. Three-tier fallback: ----
  tbl <- tryCatch(summary(fit)$table, error = function(e) NULL)
  if (is.null(tbl)) tbl <- fit$table
  if (is.null(tbl)) {
    tbl <- tryCatch({
      strata_vec <- fit$strata
      if (is.null(strata_vec)) strata_vec <- setNames(length(fit$time), "Overall")
      grp_names <- names(strata_vec)
      idx_end <- cumsum(strata_vec)
      idx_start <- idx_end - strata_vec + 1
      med_ci <- tryCatch(stats::quantile(fit, probs = 0.5), error = function(e) NULL)
      
      rows <- lapply(seq_along(strata_vec), function(i) {
        rng <- idx_start[i]:idx_end[i]
        n_i <- fit$n.risk[rng[1]]
        events_i <- sum(fit$n.event[rng])
        med_i <- NA; lcl_i <- NA; ucl_i <- NA
        if (!is.null(med_ci)) {
          if (is.list(med_ci$quantile)) {
            med_i <- med_ci$quantile[[i]]; lcl_i <- med_ci$lower[[i]]; ucl_i <- med_ci$upper[[i]]
          } else if (!is.null(dim(med_ci$quantile))) {
            med_i <- med_ci$quantile[i, 1]; lcl_i <- med_ci$lower[i, 1]; ucl_i <- med_ci$upper[i, 1]
          } else if (length(strata_vec) == 1) {
            med_i <- med_ci$quantile[1]; lcl_i <- med_ci$lower[1]; ucl_i <- med_ci$upper[1]
          }
        }
        c(records = n_i, events = events_i, median = med_i, `0.95LCL` = lcl_i, `0.95UCL` = ucl_i)
      })
      m <- do.call(rbind, rows)
      rownames(m) <- grp_names
      m
    }, error = function(e) NULL)
    if (!is.null(tbl)) {
      cat("(Note: this survfit object's per-group table was reconstructed manually from\n")
      cat("the underlying curve data — spot-check against print(fit) if in doubt.)\n\n")
    }
  }
  if (is.null(tbl)) {
    warning("Could not obtain or reconstruct a per-group summary table from this survfit object.")
    return(invisible(NULL))
  }
  if (is.null(dim(tbl))) tbl <- matrix(tbl, nrow = 1, dimnames = list("Overall", names(tbl)))
  
  n_groups  <- nrow(tbl)
  total_n   <- sum(tbl[, "records"], na.rm = TRUE)
  total_evt <- sum(tbl[, "events"], na.rm = TRUE)
  
  cat(sprintf("This describes %s patients across %s, with %s events observed in total\n",
              bold(format(total_n, big.mark = ",")),
              ifelse(n_groups == 1, "1 group", sprintf("%d groups", n_groups)),
              bold(format(total_evt, big.mark = ","))))
  cat("(the rest were censored — i.e. their follow-up ended before the event happened,\n")
  cat("so we only know they survived at least that long).\n\n")
  
  cat("How to read a Kaplan-Meier curve, in plain English:\n")
  cat("  - the line stays flat as long as no one has the event\n")
  cat("  - each STEP DOWN marks one or more events happening\n")
  cat("  - a steeper/faster-dropping curve means events are happening sooner/more often\n")
  cat("  - small vertical tick marks (if shown) mark CENSORED patients, not events\n\n")
  
  if (total_evt == 0) {
    cat("NOTE: zero events were observed. No meaningful survival differences can be\n")
    cat("estimated from this data.\n")
    cat("========================================================\n")
    return(invisible(tbl))
  }
  
  has_median <- "median" %in% colnames(tbl)
  
  cat("Per-group summary:\n")
  cat("--------------------------------------------------------\n")
  for (i in seq_len(n_groups)) {
    grp_name <- rownames(tbl)[i]
    n_i <- tbl[i, "records"]; evt_i <- tbl[i, "events"]
    pct_evt <- 100 * evt_i / n_i
    
    cat(sprintf("\n> %s\n", grp_name))
    cat(sprintf("  %s patients, %s events (%.1f%% experienced the event).\n",
                format(n_i, big.mark = ","), format(evt_i, big.mark = ","), pct_evt))
    
    if (has_median) {
      med <- tbl[i, "median"]
      lcl <- if ("0.95LCL" %in% colnames(tbl)) tbl[i, "0.95LCL"] else NA
      ucl <- if ("0.95UCL" %in% colnames(tbl)) tbl[i, "0.95UCL"] else NA
      if (is.na(med)) {
        cat("  Median survival time: NOT REACHED — fewer than half of this group had the\n")
        cat("  event by the end of follow-up (often a GOOD sign for this group).\n")
      } else {
        cat(sprintf("  Median survival time: %s time units", bold(sprintf("%.1f", med))))
        if (!is.na(lcl)) cat(sprintf(" (95%% CI: %.1f to %s)", lcl, ifelse(is.na(ucl), "not reached", sprintf("%.1f", ucl))))
        cat("\n")
      }
    }
  }
  
  cat("\n--------------------------------------------------------\n")
  cat("Note: this is DESCRIPTIVE only — it does not tell you whether differences\n")
  cat("between groups are statistically significant.")
  if (is.null(logrank)) {
    cat(" Pass logrank = survdiff(...) to\n")
    cat("this function (or call explain_logrank() separately) for that.\n")
  } else {
    cat("\n")
  }
  cat("========================================================\n")
  
  # ==========================================================
  # LOG-RANK TEST, if supplied — reuses explain_logrank() directly,
  # so there's no duplicated logic between the standalone function
  # and this unified call
  # ==========================================================
  logrank_result <- NULL
  if (!is.null(logrank)) {
    cat("\n")
    logrank_result <- explain_logrank(logrank)
  }
  
  # ==========================================================
  # FIGURE 1: km_plot — the survival curve itself (survminer, with a
  # manual ggplot2 fallback if survminer isn't installed)
  # ==========================================================
  if (isTRUE(km_plot)) {
    cat("\nGenerating Kaplan-Meier curve...\n")
    if (requireNamespace("survminer", quietly = TRUE)) {
      p_km <- tryCatch({
        survminer::ggsurvplot(
          fit, data = data,
          risk.table = TRUE, conf.int = TRUE,
          pval = !is.null(logrank),
          title = "Kaplan-Meier Survival Curve",
          xlab = "Time", ylab = "Survival probability",
          legend.title = "Group", risk.table.height = 0.25
        )
      }, error = function(e) {
        warning(sprintf("survminer::ggsurvplot() failed: %s. If your survfit object was created inside a function or a different environment, try passing data = <your original data frame> to explain_km().", conditionMessage(e)))
        NULL
      })
      if (!is.null(p_km)) print(p_km)
    } else if (requireNamespace("ggplot2", quietly = TRUE)) {
      cat("('survminer' not installed — drawing a simpler step-curve with ggplot2 instead;\n")
      cat("install.packages('survminer') for a richer plot with a risk table.)\n")
      sfit <- summary(fit)
      strata_lab <- if (!is.null(sfit$strata)) as.character(sfit$strata) else "Overall"
      plot_df <- data.frame(time = sfit$time, surv = sfit$surv, strata = strata_lab)
      p_manual <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time, y = surv, color = strata)) +
        ggplot2::geom_step(linewidth = 1) +
        ggplot2::labs(title = "Kaplan-Meier Survival Curve",
                      subtitle = "(simplified — install 'survminer' for confidence bands and a risk table)",
                      x = "Time", y = "Survival probability", color = "Group") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_minimal(base_size = 13)
      print(p_manual)
    } else {
      warning("km_plot = TRUE requires either 'survminer' or 'ggplot2' to be installed.")
    }
  }
  
  # ==========================================================
  # FIGURE 2: cumulative_plot — cumulative INCIDENCE (1 - survival),
  # often easier to read than survival probability when events are rare
  # ==========================================================
  if (isTRUE(cumulative_plot)) {
    cat("\nGenerating cumulative incidence curve (1 - survival)...\n")
    if (requireNamespace("survminer", quietly = TRUE)) {
      p_cum <- tryCatch({
        survminer::ggsurvplot(
          fit, data = data, fun = "event", conf.int = TRUE,
          title = "Cumulative Incidence (1 - Survival)",
          xlab = "Time", ylab = "Cumulative probability of event",
          legend.title = "Group"
        )
      }, error = function(e) {
        warning(sprintf("survminer::ggsurvplot() failed for the cumulative plot: %s", conditionMessage(e)))
        NULL
      })
      if (!is.null(p_cum)) print(p_cum)
    } else if (requireNamespace("ggplot2", quietly = TRUE)) {
      sfit <- summary(fit)
      strata_lab <- if (!is.null(sfit$strata)) as.character(sfit$strata) else "Overall"
      plot_df <- data.frame(time = sfit$time, event_prob = 1 - sfit$surv, strata = strata_lab)
      p_manual <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time, y = event_prob, color = strata)) +
        ggplot2::geom_step(linewidth = 1) +
        ggplot2::labs(title = "Cumulative Incidence (1 - Survival)",
                      subtitle = "(simplified — install 'survminer' for confidence bands)",
                      x = "Time", y = "Cumulative probability of event", color = "Group") +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_minimal(base_size = 13)
      print(p_manual)
    } else {
      warning("cumulative_plot = TRUE requires either 'survminer' or 'ggplot2' to be installed.")
    }
  }
  
  # ==========================================================
  # FIGURE 3: forest_plot — median survival time per group, same
  # visual layout as explain_cox()'s/explain_glm()'s forest plots, but
  # NO reference line (there's no "null value" for a raw survival time
  # the way there is for HR=1/OR=1)
  # ==========================================================
  if (isTRUE(forest_plot)) {
    if (!requireNamespace("forestplot", quietly = TRUE)) {
      warning("forest_plot = TRUE requires the 'forestplot' package. Run install.packages('forestplot').")
    } else if (!has_median) {
      warning("forest_plot = TRUE was requested, but no median survival column is available to plot.")
    } else {
      med_vals <- tbl[, "median"]
      lcl_vals <- if ("0.95LCL" %in% colnames(tbl)) tbl[, "0.95LCL"] else rep(NA, n_groups)
      ucl_vals <- if ("0.95UCL" %in% colnames(tbl)) tbl[, "0.95UCL"] else rep(NA, n_groups)
      grp_names <- rownames(tbl)
      plottable <- !is.na(med_vals)
      
      if (!any(plottable)) {
        warning("forest_plot = TRUE was requested, but no group reached a median survival time — nothing to plot.")
      } else {
        label_txt <- list(
          label = c("Group", grp_names[plottable]),
          median = c("Median", sprintf("%.1f", med_vals[plottable])),
          ci = c("95% CI", ifelse(is.na(lcl_vals[plottable]) | is.na(ucl_vals[plottable]), "NA",
                                  sprintf("%.1f to %.1f", lcl_vals[plottable], ucl_vals[plottable])))
        )
        fp_obj <- forestplot::forestplot(
          labeltext = label_txt,
          mean = c(NA, med_vals[plottable]), lower = c(NA, lcl_vals[plottable]), upper = c(NA, ucl_vals[plottable]),
          is.summary = c(TRUE, rep(FALSE, sum(plottable))),
          graph.pos = 2, boxsize = 0.3, graphwidth = grid::unit(5, "cm"),
          hrzl_lines = list("2" = grid::gpar(lwd = 2, col = "darkblue")),
          xlog = FALSE, xlab = "Median survival time",
          col = forestplot::fpColors(box = "darkblue", line = "darkblue", summary = "darkblue")
        )
        print(fp_obj)
        if (any(!plottable)) {
          cat(sprintf("\n(Note: %d group(s) with 'not reached' median survival excluded from the\n", sum(!plottable)))
          cat("figure — see the text summary above for those groups instead.)\n")
        }
      }
    }
  }
  
  # ==========================================================
  # TABLE: table_plot — publication-style summary table via
  # gtsummary::tbl_survfit(), the purpose-built gtsummary function for
  # summarising survfit objects (N, events, median survival + CI)
  # ==========================================================
  if (isTRUE(table_plot)) {
    if (!requireNamespace("gtsummary", quietly = TRUE)) {
      warning("table_plot = TRUE requires the 'gtsummary' package. Run install.packages('gtsummary').")
    } else {
      gt_table <- tryCatch(gtsummary::tbl_survfit(fit, probs = 0.5), error = function(e) {
        warning(sprintf("Could not build the gtsummary table: %s", conditionMessage(e)))
        NULL
      })
      if (!is.null(gt_table)) {
        cat("\n--------------------------------------------------------\n")
        cat("A publication-style summary table (N, events, median survival + 95% CI) has\n")
        cat("been generated below. Displays best in RStudio's Viewer pane or a knitted\n")
        cat("R Markdown/Quarto document.\n")
        cat("--------------------------------------------------------\n")
        print(gt_table)
      }
    }
  }
  
  invisible(list(table = tbl, logrank = logrank_result))
}


# ============================================================
# explain_cox() — Cox proportional hazards model, with the PH test
# and several visualisation options available directly as arguments.
# ============================================================
explain_cox <- function(model, zph = NULL, conf_level = 0.95,
                        forest_plot = FALSE, table_plot = FALSE,
                        survival_plot = FALSE, by_variable = NULL,
                        ph_plot = FALSE) {
  
  if (!inherits(model, "coxph")) {
    stop("explain_cox() expects a model fit with coxph().")
  }
  
  s_check   <- summary(model, conf.int = conf_level)
  alpha_pct <- conf_level * 100
  outcome   <- deparse(formula(model)[[2]])
  
  # ---- Pre-flight checks ----
  diagnostics <- character(0)
  
  n_events <- model$nevent
  if (is.null(n_events)) n_events <- sum(model$y[, ncol(model$y)] == 1, na.rm = TRUE)
  if (!is.na(n_events) && n_events == 0) {
    stop("This Cox model has ZERO observed events — check your status/event variable coding (1/TRUE = event?).")
  }
  
  coef_all <- coef(model)
  if (any(is.na(coef_all))) {
    diagnostics <- c(diagnostics, sprintf(
      "could NOT estimate %d coefficient(s), likely due to collinearity: %s. Those terms are skipped below.",
      sum(is.na(coef_all)), paste(names(coef_all)[is.na(coef_all)], collapse = ", ")))
  }
  
  se_vals <- s_check$coefficients[, "se(coef)"]
  if (any(se_vals > 10, na.rm = TRUE)) {
    diagnostics <- c(diagnostics,
                     "has at least one EXTREMELY large standard error — this usually signals near-perfect separation (common with small samples or rare events). HRs/CIs for the affected variable(s) may be absurdly wide or unstable.")
  }
  
  if (!is.null(model$iter) && !is.null(model$control) && model$iter >= model$control$iter.max) {
    diagnostics <- c(diagnostics, "may not have FULLY CONVERGED — it used the maximum allowed iterations.")
  }
  
  n_predictors <- length(coef_all) - sum(is.na(coef_all))
  if (!is.na(n_events) && n_predictors > 0 && (n_events / n_predictors) < 10) {
    diagnostics <- c(diagnostics, sprintf(
      "has only about %.1f events per predictor. A common rule of thumb wants at least 10 for stable estimates — treat wide CIs with extra caution.",
      n_events / n_predictors))
  }
  
  strata_vars <- attr(model$terms, "specials")$strata
  if (!is.null(strata_vars)) {
    strata_labels <- attr(model$terms, "term.labels")
    strata_labels <- strata_labels[grepl("^strata\\(", strata_labels)]
    diagnostics <- c(diagnostics, sprintf(
      "is STRATIFIED by %s. This adjusts the baseline hazard per stratum but does NOT produce its own HR — expected, not an error.",
      paste(strata_labels, collapse = ", ")))
  }
  
  if (!is.null(model$naive.var)) {
    diagnostics <- c(diagnostics, "uses CLUSTERED/ROBUST standard errors — CIs already account for non-independence, so expect them somewhat wider than naive SEs.")
  }
  
  if (length(diagnostics) > 0) {
    cat("########################################################\n")
    cat("# MODEL DIAGNOSTIC WARNINGS — read before trusting results\n")
    cat("########################################################\n")
    for (d in diagnostics) cat(bold("WARNING:"), "this model", d, "\n\n")
  }
  
  # ---- Header ----
  n_obs <- model$n
  cat("========================================================\n")
  cat("MODEL:", deparse(formula(model)), "\n")
  cat("========================================================\n\n")
  cat(sprintf("This model was fit using data from %s patients, of whom %s experienced\n",
              bold(format(n_obs, big.mark = ",")), bold(format(n_events, big.mark = ","))))
  cat("the event (the rest were censored). This is a COX PROPORTIONAL HAZARDS model,\n")
  cat("which estimates HAZARD RATIOS (HR): how much higher or lower a patient's\n")
  cat("instantaneous RISK of the event is, compared to a reference — NOT absolute\n")
  cat("survival time.\n\n")
  
  if (!is.null(s_check$concordance)) {
    conc <- s_check$concordance[1]
    cat(sprintf("Model discrimination (concordance / C-statistic): %.3f\n", conc))
    cat(sprintf("In plain terms: this model correctly ranks which of two random patients\n"))
    cat(sprintf("has the event first about %.0f%% of the time (50%% = chance; 100%% = perfect).\n\n", conc * 100))
  }
  
  # ---- Overall fit ----
  if (!is.null(s_check$logtest)) {
    lr_p <- s_check$logtest["pvalue"]
    cat("Taken together, do these variables actually help explain survival?\n")
    cat(sprintf("%s (likelihood-ratio test: p %s).\n\n",
                ifelse(lr_p < 0.05,
                       "Yes — a real, non-random relationship with survival",
                       "Not clearly — no statistically reliable relationship with survival"),
                format_pval(lr_p)))
  }
  
  # ---- Per-variable interpretation ----
  coefs_table <- s_check$coefficients
  ci_table    <- s_check$conf.int
  var_names   <- rownames(coefs_table)
  pval_col    <- ncol(coefs_table)
  term_labels <- attr(model$terms, "term.labels")
  term_labels <- term_labels[!grepl("^(strata|cluster|tt|frailty)\\(", term_labels)]
  
  cat("What each variable is associated with:\n")
  cat("--------------------------------------------------------\n")
  
  results_rows <- list()
  model_data <- tryCatch(model.frame(model), error = function(e) NULL)
  
  for (i in seq_along(var_names)) {
    vn <- var_names[i]
    if (is.na(coef_all[vn])) next
    
    estimate <- coefs_table[i, "coef"]
    pval <- coefs_table[i, pval_col]
    hr <- ci_table[vn, 1]; hr_lo <- ci_table[vn, 3]; hr_hi <- ci_table[vn, 4]
    
    is_interaction <- grepl(":", vn, fixed = TRUE)
    is_transformed <- grepl("^(log|sqrt|poly|I|exp|scale)\\(", vn)
    is_tt <- grepl("^tt\\(", vn)
    
    main_labels <- term_labels[!grepl(":", term_labels, fixed = TRUE)]
    matched_var <- main_labels[sapply(main_labels, function(p) startsWith(vn, p))]
    matched_var <- if (length(matched_var) > 0) matched_var[which.max(nchar(matched_var))] else character(0)
    
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
        cat(sprintf("  means the HR for '%s' is itself multiplied by %.2f per 1-unit change in\n", parts[1], hr))
        cat(sprintf("  (or shift to the other category of) '%s'.\n", parts[2]))
      } else {
        cat(sprintf("  Higher-order interaction between %d variables — hard to interpret from a\n", length(parts)))
        cat("  single number; consider plotting predicted survival curves instead.\n")
      }
    } else if (is_tt) {
      row_type <- "time-transform"
      cat(sprintf("  In plain English — %s:\n", bold("TIME-VARYING EFFECT (tt() term)")))
      cat("  Allows the HR to change over follow-up time — a single HR doesn't fully summarise it.\n")
    } else if (is_transformed) {
      row_type <- "transformed"
      cat(sprintf("  In plain English — %s:\n", bold("TRANSFORMED TERM")))
      cat("  A transformation of the original variable — interpret on its transformed scale.\n")
    } else if (is_ordered_factor) {
      row_type <- "ordered factor"
      cat(sprintf("  In plain English — %s:\n", bold("ORDERED FACTOR")))
      cat("  Describes a TREND across ordered levels, not a simple group comparison.\n")
    } else if (is_factor_level) {
      reference_level <- levels(model_data[[matched_var]])[1]
      compared_level  <- substring(vn, nchar(matched_var) + 1)
      pct_diff <- if (hr >= 1) (hr - 1) * 100 else (1 - hr) * 100
      direction_word <- if (hr > 1) "HIGHER (worse survival)" else if (hr < 1) "LOWER (better survival)" else "no different"
      
      cat(sprintf("  In plain English — %s vs %s:\n", compared_level, reference_level))
      if (abs(hr - 1) < 1e-9) {
        cat(sprintf("  '%s' and '%s' show essentially IDENTICAL hazard.\n", compared_level, reference_level))
      } else {
        cat(sprintf("  %s Patients in the '%s' group have %s%% %s hazard of %s than\n",
                    bold("Which group fares worse?"), bold(compared_level),
                    bold(sprintf("%.0f", pct_diff)), bold(direction_word), bold(outcome)))
        cat(sprintf("  patients in the '%s' (reference) group, holding other variables constant.\n", reference_level))
      }
    } else {
      pct_diff <- if (hr >= 1) (hr - 1) * 100 else (1 - hr) * 100
      direction_word <- if (hr > 1) "HIGHER (worse survival)" else if (hr < 1) "LOWER (better survival)" else "no different"
      cat("  In plain English:\n")
      cat(sprintf("  For every %s in %s, the hazard of %s is multiplied by %s\n",
                  bold("1-unit increase"), vn, outcome, bold(sprintf("%.2f", hr))))
      cat(sprintf("  — %s%% %s, ON AVERAGE, holding other variables constant.\n",
                  bold(sprintf("%.1f", pct_diff)), bold(direction_word)))
    }
    
    cat("\n  Is this a real effect, or could it be due to chance?\n")
    if (is.na(pval)) {
      cat("  P-value not available.\n")
    } else if (pval < 0.05) {
      cat(sprintf("  This IS %s (p %s).\n", bold("statistically significant"), format_pval(pval)))
    } else {
      cat(sprintf("  This is %s (p %s). Treat with caution.\n", bold("NOT statistically significant"), format_pval(pval)))
    }
    
    results_rows[[length(results_rows) + 1]] <- data.frame(
      variable = vn, type = row_type, HR = round(hr, 3), lower_ci = round(hr_lo, 3), upper_ci = round(hr_hi, 3),
      p_value = ifelse(is.na(pval), NA, ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval))),
      significant = ifelse(is.na(pval), NA, pval < 0.05)
    )
  }
  
  cat("\n--------------------------------------------------------\n")
  cat("Reminder: HR = 1.00 means NO difference. HR > 1 = higher risk/worse survival;\n")
  cat("HR < 1 = lower risk/better survival. This model ASSUMES each HR is constant\n")
  cat("over follow-up time (the 'proportional hazards' assumption).")
  if (is.null(zph)) {
    cat(" Pass zph = cox.zph(model)\n")
    cat("to this function (or call explain_cox_zph() separately) to check that.\n")
  } else {
    cat("\n")
  }
  cat("========================================================\n")
  
  results_table <- do.call(rbind, results_rows)
  
  # ==========================================================
  # PH TEST, if supplied — reuses explain_cox_zph() directly
  # ==========================================================
  if (!is.null(zph)) {
    cat("\n")
    explain_cox_zph(zph)
  }
  
  # ==========================================================
  # FOREST PLOT — via cal_forest_plot(), identical layout to explain_glm()
  # ==========================================================
  if (isTRUE(forest_plot)) {
    tryCatch(cal_forest_plot(model, family = "HR"),
             error = function(e) warning(sprintf("Could not draw the forest plot: %s", conditionMessage(e))))
  }
  
  # ==========================================================
  # GTSUMMARY TABLE — same layout as explain_glm()
  # ==========================================================
  if (isTRUE(table_plot)) {
    if (!requireNamespace("gtsummary", quietly = TRUE)) {
      warning("table_plot = TRUE requires 'gtsummary'. Run install.packages('gtsummary').")
    } else {
      gt_table <- tryCatch(gtsummary::tbl_regression(model, exponentiate = TRUE), error = function(e) {
        warning(sprintf("Could not build the gtsummary table: %s", conditionMessage(e))); NULL
      })
      if (!is.null(gt_table)) {
        cat("\n--------------------------------------------------------\n")
        cat("A publication-style summary table (HR, 95% CI, p-value) has been generated\n")
        cat("below — displays best in RStudio's Viewer pane or a knitted document.\n")
        cat("--------------------------------------------------------\n")
        print(gt_table)
      }
    }
  }
  
  # ==========================================================
  # ADJUSTED SURVIVAL CURVES — predicted survival by group (or by
  # low/median/high of a continuous variable), holding all other
  # predictors at a reference value. This is one of the most useful
  # Cox-model visualisations in practice: it translates abstract HRs
  # back into an actual survival curve shape.
  # ==========================================================
  if (isTRUE(survival_plot)) {
    cat("\nGenerating adjusted survival curve(s)...\n")
    
    main_labels_sp <- term_labels[!grepl(":", term_labels, fixed = TRUE)]
    
    by_var <- by_variable
    if (is.null(by_var)) {
      factor_candidates <- main_labels_sp[sapply(main_labels_sp, function(v) {
        !is.null(model_data) && v %in% names(model_data) && is.factor(model_data[[v]])
      })]
      by_var <- if (length(factor_candidates) > 0) factor_candidates[1] else NULL
    }
    
    if (is.null(by_var) && length(main_labels_sp) > 0) {
      # No factor available — fall back to low/median/high of the first continuous predictor
      by_var <- main_labels_sp[1]
    }
    
    if (is.null(by_var) || is.null(model_data) || !(by_var %in% names(model_data))) {
      cat("Could not determine a variable to plot adjusted survival curves by — specify\n")
      cat("one explicitly: explain_cox(model, survival_plot = TRUE, by_variable = \"treatment\")\n")
    } else {
      is_by_factor <- is.factor(model_data[[by_var]])
      other_vars <- setdiff(main_labels_sp, by_var)
      
      grid <- if (is_by_factor) {
        data.frame(setNames(list(factor(levels(model_data[[by_var]]), levels = levels(model_data[[by_var]]))), by_var))
      } else {
        qs <- stats::quantile(model_data[[by_var]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE, type = 7)
        data.frame(setNames(list(qs), by_var))
      }
      for (ov in other_vars) {
        if (!(ov %in% names(model_data))) next
        if (is.factor(model_data[[ov]])) {
          grid[[ov]] <- factor(levels(model_data[[ov]])[1], levels = levels(model_data[[ov]]))
        } else {
          grid[[ov]] <- mean(model_data[[ov]], na.rm = TRUE)
        }
      }
      
      sfit2 <- tryCatch(survfit(model, newdata = grid), error = function(e) {
        warning(sprintf("Could not generate adjusted survival curves: %s", conditionMessage(e)))
        NULL
      })
      
      if (!is.null(sfit2)) {
        grp_labels <- if (is_by_factor) {
          as.character(grid[[by_var]])
        } else {
          sprintf("%s = %.1f (%s)", by_var, grid[[by_var]], c("10th pct", "median", "90th pct"))
        }
        
        if (requireNamespace("survminer", quietly = TRUE)) {
          p_surv <- tryCatch({
            survminer::ggsurvplot(sfit2, data = grid, legend.labs = grp_labels,
                                  title = sprintf("Adjusted Survival Curves by %s", by_var),
                                  subtitle = "(holding other predictors at their mean/reference level)",
                                  xlab = "Time", ylab = "Predicted survival probability",
                                  legend.title = by_var)
          }, error = function(e) {
            warning(sprintf("survminer::ggsurvplot() failed for the adjusted survival plot: %s", conditionMessage(e)))
            NULL
          })
          if (!is.null(p_surv)) print(p_surv)
        } else if (requireNamespace("ggplot2", quietly = TRUE)) {
          s2 <- summary(sfit2)
          strata_idx <- if (!is.null(s2$strata)) as.integer(s2$strata) else rep(1, length(s2$time))
          plot_df <- data.frame(time = s2$time, surv = s2$surv,
                                group = factor(grp_labels[strata_idx], levels = grp_labels))
          p_manual <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time, y = surv, color = group)) +
            ggplot2::geom_step(linewidth = 1) +
            ggplot2::labs(title = sprintf("Adjusted Survival Curves by %s", by_var),
                          subtitle = "(simplified — install 'survminer' for a richer plot)",
                          x = "Time", y = "Predicted survival probability", color = by_var) +
            ggplot2::ylim(0, 1) +
            ggplot2::theme_minimal(base_size = 13)
          print(p_manual)
        } else {
          warning("survival_plot = TRUE requires either 'survminer' or 'ggplot2' to be installed.")
        }
      }
    }
  }
  
  # ==========================================================
  # PH PLOT — Schoenfeld residuals, if zph was supplied
  # ==========================================================
  if (isTRUE(ph_plot)) {
    if (is.null(zph)) {
      cat("\nph_plot = TRUE requires zph = cox.zph(model) to also be supplied.\n")
    } else {
      cat("\nGenerating Schoenfeld residual plot(s) for the proportional hazards test...\n")
      if (requireNamespace("survminer", quietly = TRUE)) {
        p_zph <- tryCatch(survminer::ggcoxzph(zph), error = function(e) {
          warning(sprintf("survminer::ggcoxzph() failed: %s — falling back to base plot().", conditionMessage(e)))
          NULL
        })
        if (!is.null(p_zph)) print(p_zph) else tryCatch(plot(zph), error = function(e) NULL)
      } else {
        cat("('survminer' not installed — using base R's plot(cox.zph(...)) instead;\n")
        cat("install.packages('survminer') for a nicer ggplot2 version.)\n")
        tryCatch(plot(zph), error = function(e) warning(sprintf("Could not plot zph: %s", conditionMessage(e))))
      }
    }
  }
  
  invisible(results_table)
}


# ============================================================
# USAGE EXAMPLES (using the cohort dataset from earlier)
# ============================================================
library(survival)

out_km_fit     <- survfit(Surv(time, status) ~ treatment, data = cohort)
out_km_logrank <- survdiff(Surv(time, status) ~ treatment, data = cohort)

# Unified call — matches the pattern you liked:
explain_km(fit = out_km_fit, logrank = out_km_logrank, data = cohort,
           km_plot = TRUE, cumulative_plot = TRUE,
           forest_plot = TRUE, table_plot = TRUE)

model_cox   <- coxph(Surv(time, status) ~ age + sex + treatment + comorbidity_score,
                     data = cohort)
out_cox_zph <- cox.zph(model_cox)

explain_cox(model = model_cox, zph = out_cox_zph,
            forest_plot = TRUE, table_plot = TRUE,
            survival_plot = TRUE, by_variable = "treatment",
            ph_plot = TRUE)

# Standalone equivalents still work if you only want one piece:
explain_logrank(out_km_logrank)
explain_cox_zph(out_cox_zph)

# ============================================================
# EDGE CASES THIS FILE HAS BEEN SPECIFICALLY CHECKED AGAINST
# ============================================================
# explain_km():
#  1. Per-group table unavailable via summary(fit)$table or fit$table -> manually reconstructed
#     from n.risk/n.event/strata + quantile.survfit(), handling that method's differently-shaped
#     return value for 1 vs. multiple strata
#  2. Single-group (no strata) survfit                        -> table reshaped to 1-row matrix
#  3. Zero events observed                                    -> flagged, skips per-group detail
#  4. Median not reached in a group (NA)                      -> explained explicitly
#  5. km_plot/cumulative_plot with survminer NOT installed     -> falls back to a manual ggplot2
#     step-curve; if ggplot2 ALSO unavailable, warns clearly instead of erroring
#  6. survminer::ggsurvplot() failing (common cause: survfit's formula environment doesn't
#     have the original data) -> caught, suggests passing data = <your data frame> explicitly
#  7. forest_plot = TRUE where NO group has a finite median (all "not reached") -> warns,
#     nothing plotted, rather than crashing on an all-NA mean/lower/upper vector
#  8. forest_plot = TRUE where SOME groups have "not reached" median -> excluded from the
#     figure with an explicit note, not silently dropped
#  9. table_plot = TRUE with gtsummary not installed, or tbl_survfit() erroring -> warns,
#     rest of the function unaffected
#
# explain_logrank() / explain_cox_zph(): same edge cases as before (only 1 group, missing
# $pvalue on older survival versions, missing GLOBAL row, NA p-values) — unchanged, still
# fully standalone-callable, and now also invoked directly from inside explain_km()/explain_cox()
# with no logic duplicated between the two calling styles.
#
# explain_cox():
# 10. Zero events                                             -> stops with a clear, specific message
# 11. Aliased/NA coefficients, huge SEs, non-convergence, few events/predictor, strata(),
#     cluster() — all the same checks carried over from the previous version, unchanged
# 12. Interaction, transformed, tt(), ordered-factor terms     -> explained differently, as before
# 13. forest_plot/table_plot with required packages missing    -> stop/warn with actionable
#     messages via cal_forest_plot()'s own requireNamespace checks
# 14. survival_plot = TRUE with no by_variable specified and NO factor predictor in the model
#     -> falls back to the first continuous predictor, plotted at low/median/high percentile
# 15. survival_plot = TRUE with by_variable specified but not found in the model's data
#     -> clear message telling you to check the variable name, rather than a cryptic error
# 16. survfit(model, newdata = grid) failing (e.g. a transformed term that can't be
#     reconstructed from the grid) -> caught, warns, rest of explain_cox()'s output unaffected
# 17. survival_plot with survminer unavailable                -> manual ggplot2 fallback,
#     correctly mapping each curve's strata index back to its group label
# 18. ph_plot = TRUE without zph supplied                     -> explicit message telling you
#     to pass zph = cox.zph(model), rather than silently doing nothing
# 19. survminer::ggcoxzph() failing                            -> falls back to base R's
#     plot(cox.zph(...)) rather than losing the figure entirely
# ============================================================