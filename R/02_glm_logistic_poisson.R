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
#    95% CIs, via cal_forest_plot() (gtsummary + forestplot packages).
#    The reference line for "no effect" sits at 1 on a log-scaled axis
#    for ratio measures — 1 is the null value for ANY ratio (OR, RR, or
#    IRR); this is different from 0, which is only the null value on
#    the raw, un-exponentiated log-coefficient scale.
#  - Optionally builds a publication-style gtsummary table (table_plot = TRUE)
#  - Optionally plots how the predicted outcome diverges across a 2-way
#    interaction term (show_interaction = TRUE) — e.g. age on the x-axis,
#    one line per sex, showing how the two variables' effects combine
#
# REQUIRED PACKAGES for forest_plot/table_plot: forestplot, gtsummary, dplyr
#   install.packages(c("forestplot", "gtsummary", "dplyr"))
# REQUIRED PACKAGE for show_interaction: ggplot2
#   install.packages("ggplot2")
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


# ============================================================
# cal_forest_plot() — gtsummary + forestplot-based forest plot.
# This is your function, used as-is, with two additions:
#   - explicit requireNamespace checks for dplyr/gtsummary (matching
#     the check you already had for forestplot), so a missing
#     package gives a clear message here rather than a cryptic
#     "could not find function" error from inside the call
#   - everything else (logic, defaults, argument names) is unchanged
# ============================================================
cal_forest_plot <- function(fit_or_tbl,
                            family = "OR",
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
  
  # NOTE: deliberately NOT calling library(dplyr)/library(forestplot)/library(gtsummary)
  # here. Attaching dplyr is exactly what caused dplyr::select to get masked by
  # MASS::select (or vice versa) when this ran inside a script that also loads MASS
  # (e.g. for glm.nb() earlier in the workflow). Every call below is explicitly
  # namespaced (dplyr::..., gtsummary::...) instead, so this function works correctly
  # regardless of what else is attached to the search path or in what order.
  
  # If input is a GLM/fit, convert to tbl_regression internally
  if (inherits(fit_or_tbl, "glm") || inherits(fit_or_tbl, "lm")) {
    x <- gtsummary::tbl_regression(fit_or_tbl, exponentiate = exponentiate)
  } else if (inherits(fit_or_tbl, c("tbl_regression", "tbl_uvregression"))) {
    x <- fit_or_tbl
  } else {
    stop("Input must be a glm/lm object or tbl_regression/tbl_uvregression object", call. = FALSE)
  }
  
  # Determine log scale and x-axis label
  family <- toupper(family)
  xlog <- family %in% c("OR", "RR", "IRR")  # log scale for ratio measures
  xlab <- switch(family,
                 "OR"  = "Odds Ratio",
                 "RR"  = "Risk Ratio",
                 "IRR" = "Incidence Rate Ratio",
                 "RISK"= "Risk",
                 "Estimate")  # default
  
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
  
  # Draw forest plot
  # IMPORTANT: forestplot::forestplot() returns an object that only actually
  # renders via R's auto-print mechanism — which fires for an untouched
  # top-level console expression, but NOT when this function is called from
  # inside another function (e.g. from explain_glm(), inside a tryCatch{}).
  # In that nested context, the plot object was silently built and discarded
  # with nothing ever drawn — which is why the plot appeared blank. Printing
  # explicitly here guarantees it draws regardless of how/where this is called,
  # and returning invisibly avoids a redundant second draw if you also print
  # the returned object yourself (e.g. `fp <- cal_forest_plot(...); fp`).
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
# .build_glm_interaction_plot() — internal helper for show_interaction.
# Builds a prediction grid across two interacting variables (holding
# every OTHER predictor at a reference value: mean for continuous,
# reference/first level for factors), gets predicted values on the
# RESPONSE scale (probability/rate — whatever the model's family
# implies), and plots how they diverge. Handles three cases
# differently, since one plot shape doesn't fit all of them:
#   continuous x factor   -> line plot, x = continuous, one line per
#                             factor level (e.g. age on x, one line
#                             per sex) — this is the classic case
#   factor x factor       -> grouped bar chart
#   continuous x continuous -> line plot at low/median/high (10th/
#                             50th/90th percentile) of the second variable
# Returns a ggplot object, or NULL (with a warning explaining why) if
# it can't be built.
# ============================================================
.build_glm_interaction_plot <- function(model, var1, var2, term_labels,
                                        ratio_noun, outcome, family_name, link_name,
                                        has_offset, offset_expr, offset_is_log) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("show_interaction = TRUE requires the 'ggplot2' package. Run install.packages('ggplot2') and try again.")
    return(NULL)
  }
  
  md <- model$model
  if (!(var1 %in% names(md)) || !(var2 %in% names(md))) {
    warning(sprintf("Could not build an interaction figure for %s:%s — one or both variables could not be found in the model's data.", var1, var2))
    return(NULL)
  }
  
  main_labels <- term_labels[!grepl(":", term_labels, fixed = TRUE)]
  is_factor1 <- is.factor(md[[var1]])
  is_factor2 <- is.factor(md[[var2]])
  
  # ---- Reference values for every OTHER main-effect predictor ----
  other_vars <- setdiff(main_labels, c(var1, var2))
  
  set_reference_values <- function(df) {
    for (ov in other_vars) {
      if (!(ov %in% names(md))) next  # e.g. a transformed term whose underlying variable can't be cleanly resolved — skipped, held at whatever predict()'s own defaulting does
      if (is.factor(md[[ov]])) {
        df[[ov]] <- factor(levels(md[[ov]])[1], levels = levels(md[[ov]]))
      } else {
        df[[ov]] <- mean(md[[ov]], na.rm = TRUE)
      }
    }
    df
  }
  
  # ---- Offset handling: hold the underlying exposure-time variable at its
  # median, so the plot shows predicted rates for a "typical" follow-up time
  # rather than an arbitrary/unset value that could distort the numbers ----
  offset_varname <- NA
  offset_note <- ""
  if (has_offset) {
    inner <- sub("^offset\\((.*)\\)$", "\\1", offset_expr)
    inner <- sub("^log\\((.*)\\)$", "\\1", inner)
    offset_varname <- trimws(inner)
    if (!(offset_varname %in% names(md))) offset_varname <- NA
  }
  set_offset_reference <- function(df) {
    if (!is.na(offset_varname)) {
      ref_val <- stats::median(md[[offset_varname]], na.rm = TRUE)
      df[[offset_varname]] <- ref_val
      offset_note <<- sprintf(" (assuming %s = %s, the median follow-up time)",
                              offset_varname, signif(ref_val, 3))
    }
    df
  }
  
  # ---- Upper bound for the CI ribbon: probabilities can't exceed 1 ----
  is_probability_scale <- grepl("^(binomial|quasibinomial)", family_name) &&
    link_name %in% c("logit", "log")
  upper_bound <- if (is_probability_scale) 1 else Inf
  
  y_lab <- sprintf("Predicted %s of %s", ratio_noun, outcome)
  
  get_predictions <- function(grid) {
    tryCatch(predict(model, newdata = grid, type = "response", se.fit = TRUE),
             error = function(e) NULL)
  }
  
  # ==========================================================
  # CASE 1: continuous x factor (either order) — line plot, one line per level
  # ==========================================================
  if (xor(is_factor1, is_factor2)) {
    cont_var  <- if (is_factor1) var2 else var1
    fact_var  <- if (is_factor1) var1 else var2
    
    cont_seq <- seq(min(md[[cont_var]], na.rm = TRUE), max(md[[cont_var]], na.rm = TRUE), length.out = 50)
    fact_lvls <- levels(md[[fact_var]])
    
    grid <- expand.grid(setNames(list(cont_seq, fact_lvls), c(cont_var, fact_var)), stringsAsFactors = FALSE)
    grid[[fact_var]] <- factor(grid[[fact_var]], levels = fact_lvls)
    grid <- set_reference_values(grid)
    grid <- set_offset_reference(grid)
    
    preds <- get_predictions(grid)
    if (is.null(preds)) {
      warning(sprintf("Could not generate predictions for the %s:%s interaction figure.", var1, var2))
      return(NULL)
    }
    grid$fit <- preds$fit
    grid$lo  <- pmax(0, preds$fit - 1.96 * preds$se.fit)
    grid$hi  <- pmin(upper_bound, preds$fit + 1.96 * preds$se.fit)
    
    p <- ggplot2::ggplot(grid, ggplot2::aes(x = .data[[cont_var]], y = fit,
                                            color = .data[[fact_var]], fill = .data[[fact_var]])) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::labs(title = sprintf("Interaction: %s x %s", var1, var2),
                    subtitle = sprintf("Predicted %s of %s across %s, separately by %s%s (shaded band = approx. 95%% CI)",
                                       ratio_noun, outcome, cont_var, fact_var, offset_note),
                    x = cont_var, y = y_lab, color = fact_var, fill = fact_var) +
      ggplot2::theme_minimal(base_size = 13)
    return(p)
  }
  
  # ==========================================================
  # CASE 2: factor x factor — grouped bar chart
  # ==========================================================
  if (is_factor1 && is_factor2) {
    grid <- expand.grid(setNames(list(levels(md[[var1]]), levels(md[[var2]])), c(var1, var2)),
                        stringsAsFactors = FALSE)
    grid[[var1]] <- factor(grid[[var1]], levels = levels(md[[var1]]))
    grid[[var2]] <- factor(grid[[var2]], levels = levels(md[[var2]]))
    grid <- set_reference_values(grid)
    grid <- set_offset_reference(grid)
    
    preds <- get_predictions(grid)
    if (is.null(preds)) {
      warning(sprintf("Could not generate predictions for the %s:%s interaction figure.", var1, var2))
      return(NULL)
    }
    grid$fit <- preds$fit
    grid$lo  <- pmax(0, preds$fit - 1.96 * preds$se.fit)
    grid$hi  <- pmin(upper_bound, preds$fit + 1.96 * preds$se.fit)
    
    p <- ggplot2::ggplot(grid, ggplot2::aes(x = .data[[var1]], y = fit, fill = .data[[var2]])) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi),
                             position = ggplot2::position_dodge(width = 0.8), width = 0.2) +
      ggplot2::labs(title = sprintf("Interaction: %s x %s", var1, var2),
                    subtitle = sprintf("Predicted %s of %s by %s and %s%s (error bars = approx. 95%% CI)",
                                       ratio_noun, outcome, var1, var2, offset_note),
                    x = var1, y = y_lab, fill = var2) +
      ggplot2::theme_minimal(base_size = 13)
    return(p)
  }
  
  # ==========================================================
  # CASE 3: continuous x continuous — line plot at low/median/high
  # (10th/50th/90th percentile) of var2, since a full 2D surface is
  # harder to read at a glance than a handful of representative lines
  # ==========================================================
  seq1 <- seq(min(md[[var1]], na.rm = TRUE), max(md[[var1]], na.rm = TRUE), length.out = 50)
  qs <- stats::quantile(md[[var2]], probs = c(0.1, 0.5, 0.9), na.rm = TRUE, type = 7)
  qs_labels_full <- c("10th pct", "median", "90th pct")
  
  if (length(unique(qs)) < length(qs)) {
    # Skewed/discrete distributions can collapse two percentiles to the same value
    # (e.g. many tied values in comorbidity_score) — factor(levels = qs) would
    # otherwise error on duplicate levels, so de-duplicate defensively here.
    keep <- !duplicated(qs)
    qs <- qs[keep]
    qs_labels_full <- qs_labels_full[keep]
  }
  
  grid <- expand.grid(setNames(list(seq1, qs), c(var1, var2)))
  grid <- set_reference_values(grid)
  grid <- set_offset_reference(grid)
  
  preds <- get_predictions(grid)
  if (is.null(preds)) {
    warning(sprintf("Could not generate predictions for the %s:%s interaction figure.", var1, var2))
    return(NULL)
  }
  grid$fit <- preds$fit
  grid$lo  <- pmax(0, preds$fit - 1.96 * preds$se.fit)
  grid$hi  <- pmin(upper_bound, preds$fit + 1.96 * preds$se.fit)
  grid$group <- factor(grid[[var2]],
                       levels = qs,
                       labels = sprintf("%s = %.1f (%s)", var2, qs, qs_labels_full))
  
  p <- ggplot2::ggplot(grid, ggplot2::aes(x = .data[[var1]], y = fit, color = group, fill = group)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(title = sprintf("Interaction: %s x %s", var1, var2),
                  subtitle = sprintf("Predicted %s of %s across %s, at low/median/high %s%s (shaded band = approx. 95%% CI)",
                                     ratio_noun, outcome, var1, var2, offset_note),
                  x = var1, y = y_lab, color = NULL, fill = NULL) +
    ggplot2::theme_minimal(base_size = 13)
  p
}


explain_glm <- function(model, conf_level = 0.95, forest_plot = FALSE, table_plot = FALSE, show_interaction = FALSE) {
  
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
  interaction_pairs <- list()  # collected here, plotted later if show_interaction = TRUE
  
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
        
        # Resolve each piece back to its base data-column variable name (a factor
        # piece like "sexM" needs to resolve to "sex"; a continuous piece like
        # "age" already equals its own base name), same matching approach used
        # for main-effect terms above — so show_interaction knows exactly which
        # two original variables to build a prediction grid across.
        resolve_base_var <- function(piece) {
          m <- main_labels[sapply(main_labels, function(p) startsWith(piece, p))]
          if (length(m) == 0) return(NA_character_)
          m[which.max(nchar(m))]
        }
        base_var1 <- resolve_base_var(parts[1])
        base_var2 <- resolve_base_var(parts[2])
        
        if (!is.na(base_var1) && !is.na(base_var2)) {
          interaction_pairs[[length(interaction_pairs) + 1]] <- list(var1 = base_var1, var2 = base_var2, term = vn)
        } else {
          cat(sprintf("  (Note: could not automatically resolve this interaction's variables for a\n"))
          cat(sprintf("  figure — if show_interaction = TRUE, this term will be skipped.)\n"))
        }
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
  # 6. OPTIONAL INTERACTION FIGURE(S)
  # For each clean 2-way interaction found above, plots how the predicted
  # outcome (on the response scale — probability/rate, matching what the
  # rest of this function reports) diverges across the two variables,
  # holding every other predictor at a reference value (mean for
  # continuous, reference level for factors).
  # ==========================================================
  if (isTRUE(show_interaction)) {
    if (length(interaction_pairs) == 0) {
      cat("\n(show_interaction = TRUE was requested, but this model has no 2-way\n")
      cat("interaction terms to plot.)\n")
    } else {
      for (ip in interaction_pairs) {
        cat(sprintf("\nGenerating interaction figure for %s...\n", ip$term))
        p_int <- tryCatch(
          .build_glm_interaction_plot(model, ip$var1, ip$var2, term_labels,
                                      ratio_noun, outcome, family_name, link_name,
                                      has_offset, offset_expr, offset_is_log),
          error = function(e) {
            warning(sprintf("Could not build the interaction figure for %s: %s", ip$term, conditionMessage(e)))
            NULL
          }
        )
        if (!is.null(p_int)) print(p_int)
      }
    }
  }
  
  # ==========================================================
  # 7. OPTIONAL FOREST PLOT (via cal_forest_plot() — gtsummary + forestplot)
  # Reference/"no effect" line is drawn at 1 on the log-scaled axis for
  # ratio measures (OR/RR/IRR) — 1 is the null value for ANY ratio; the
  # log axis (xlog = TRUE inside cal_forest_plot) is what makes that
  # reference point sit visually centered rather than at 0.
  # ==========================================================
  if (isTRUE(forest_plot)) {
    
    # Map the ratio_label already determined earlier in this function to the
    # `family` argument cal_forest_plot() expects, so you don't have to
    # specify it separately — it stays consistent with the text output above.
    fp_family <- if (grepl("^Odds Ratio", ratio_label)) {
      "OR"
    } else if (grepl("^Incidence Rate Ratio", ratio_label)) {
      "IRR"
    } else if (grepl("^Risk Ratio", ratio_label) || grepl("^Rate Ratio", ratio_label)) {
      "RR"
    } else {
      "Estimate"  # generic/unrecognised link
    }
    
    # NOTE on exponentiate: this is intentionally left at cal_forest_plot()'s own
    # default (exponentiate = TRUE) for EVERY case, including the generic
    # "Estimate" fallback. This matches the PER-VARIABLE TEXT OUTPUT above (section
    # 5), which also unconditionally exponentiates every coefficient regardless of
    # family/link — the plot and the text need to agree with each other above all
    # else. For an unusual link (e.g. a plain gaussian/identity glm, which would
    # fall into this generic branch), exponentiating a raw coefficient does NOT
    # have a standard "ratio" interpretation — this is a known limitation carried
    # through from the original design of explain_glm()'s text output, not
    # something newly introduced by the forest plot. If you use explain_glm() on
    # a genuinely raw-scale (identity link) model, treat both the text AND the
    # plot's numbers with real caution for that specific case.
    if (fp_family == "RR" && grepl("^Rate Ratio", ratio_label)) {
      cat("\nNote: this model's ratio is a Poisson/Negative Binomial RATE RATIO, but\n")
      cat("cal_forest_plot()'s x-axis will be labeled 'Risk Ratio' (it does not\n")
      cat("distinguish rate ratios from risk ratios internally) — the numbers and log\n")
      cat("scale are correct, only that axis label's wording is imprecise for a count\n")
      cat("outcome. Edit cal_forest_plot()'s internal `xlab` switch statement if you\n")
      cat("want it to say 'Rate Ratio' specifically for this case.\n\n")
    }
    
    tryCatch({
      cal_forest_plot(model, family = fp_family)
    }, error = function(e) {
      warning(sprintf("Could not draw the forest plot via cal_forest_plot(): %s", conditionMessage(e)))
    })
  }
  
  # ==========================================================
  # 8. OPTIONAL GTSUMMARY TABLE
  # A publication-style regression table (variable, N, OR/RR/IRR, 95% CI,
  # p-value) built with gtsummary::tbl_regression(). Displays best in
  # RStudio's Viewer pane or when knitted in an R Markdown/Quarto document;
  # in a plain console it prints as a wide data frame.
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
        cat(sprintf("A publication-style summary table (%s, 95%% CI, p-value) has been\n", ratio_label))
        cat("generated below. This displays best in RStudio's Viewer pane or when\n")
        cat("knitted into an R Markdown/Quarto document — a plain R console will show\n")
        cat("it as a wide data frame instead of the formatted table.\n")
        cat("--------------------------------------------------------\n")
        print(gt_table)
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
#
# --- Poisson with offset (log link -> Incidence Rate Ratios) ---
model_pois <- glm(n_hosp ~ age + treatment + offset(log(person_time)),
                  data = cohort, family = poisson(link = "log"))
pois_results <- explain_glm(model_pois)
pois_results

# --- Negative Binomial (overdispersed counts -> IRRs, with theta note) ---
library(MASS)
model_nb <- glm.nb(n_hosp ~ age + treatment + offset(log(person_time)), data = cohort)
nb_results <- explain_glm(model_nb)
nb_results
#
# --- Interaction term example ---
model_int <- glm(event ~ age * treatment, data = cohort, family = binomial)
explain_glm(model_int)
#
# --- With the interaction figure (age on x-axis, one line per treatment group) ---
explain_glm(model_int, show_interaction = TRUE)
#
# --- A continuous x continuous interaction (age * bmi) — plotted at low/median/high bmi ---
model_int2 <- glm(event ~ age * bmi, data = cohort, family = binomial)
explain_glm(model_int2, show_interaction = TRUE)
#
# --- A factor x factor interaction (sex * treatment) — plotted as a grouped bar chart ---
model_int3 <- glm(event ~ sex * age, data = cohort, family = binomial)
explain_glm(model_int3, show_interaction = TRUE)
#
# --- With a forest plot (via cal_forest_plot(); reference line at ratio = 1) ---
explain_glm(model_logit, forest_plot = TRUE)
explain_glm(model_pois, forest_plot = TRUE)
#
# --- With a publication-style gtsummary table ---
explain_glm(model_logit, table_plot = TRUE)
#
# --- Both together ---
explain_glm(model_logit, forest_plot = TRUE, table_plot = TRUE)
#
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
# 18. forest_plot = TRUE with forestplot/gtsummary/dplyr not installed -> stops with a specific,
#     actionable message (via cal_forest_plot()'s own requireNamespace checks) rather than a
#     cryptic "could not find function" error
# 19. cal_forest_plot() erroring for any other reason (e.g. an unusual model structure it
#     doesn't handle) -> caught and reported as a warning, rest of explain_glm()'s text output
#     is unaffected since the forest plot is generated only after all text has already printed
# 20. table_plot = TRUE with gtsummary not installed -> warns with an actionable message,
#     rest of the function is unaffected
# 21. gtsummary::tbl_regression() erroring for any reason -> caught, warns instead of crashing
# 22. fp_family mapping: every ratio_label string this function can produce (OR/RR/IRR/generic)
#     is explicitly mapped to a valid cal_forest_plot() family code, so the forest plot always
#     gets a recognised family value consistent with the text output above it — never silently
#     mismatched (e.g. text calling something an IRR while the plot draws it as a plain OR)
# 23. show_interaction = TRUE with no 2-way interaction terms in the model -> a specific message
#     is printed instead of silently doing nothing
# 24. show_interaction = TRUE with ggplot2 not installed -> warns with an actionable message,
#     the rest of explain_glm()'s output is unaffected
# 25. An interaction piece (e.g. "sexM") that can't be resolved back to its original data column
#     -> flagged in the text output at the point it's detected, and that specific interaction is
#     skipped for plotting rather than silently building a wrong or empty grid
# 26. Every combination of variable types in the interaction is handled explicitly and
#     differently: continuous x factor (line plot, one line per level — the classic case),
#     factor x factor (grouped bar chart), continuous x continuous (line plot at low/median/high
#     percentile of the second variable) — never assumed to be one shape for all three
# 27. A model with an offset (e.g. Poisson rate models) -> the underlying exposure-time variable
#     is held at its MEDIAN in the prediction grid (not left unset, which would silently break
#     predict() or produce a meaningless reference value), with the assumed value stated directly
#     in the plot's subtitle
# 28. Predicted-probability confidence ribbons for binomial models are clamped to [0, 1] (a
#     probability can't exceed 1 or go below 0), rather than a symmetric Wald interval spilling
#     outside the valid range on the plot
# 29. predict() failing on the constructed grid for any reason (e.g. a transformed/derived term
#     that can't be reconstructed from the grid's raw variables) -> caught, warns with the
#     specific interaction term named, that one figure is skipped rather than crashing the
#     whole function
# 30. A transformed term (e.g. log(bmi)) appearing among the "other" predictors being held at a
#     reference value -> its underlying variable can't always be cleanly resolved by name, so it
#     is left for predict()'s own default handling rather than the function guessing incorrectly