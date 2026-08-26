

# ============================================================
# ALL RESULTS IN ONE PLACE
# Every out_* object collected into a single named list so you
# can browse them via View(results) or results$out_cox_summary etc.
# ============================================================
results <- list(
  cohort_head        = out_cohort_head,
  cohort_long_head   = out_cohort_long_head,
  lm_summary         = out_lm_summary,
  lm_confint         = out_lm_confint,
  lm_predict         = out_lm_predict,
  logit_summary      = out_logit_summary,
  logit_or           = out_logit_or,
  logit_or_ci        = out_logit_or_ci,
  poisson_summary    = out_pois_summary,
  negbin_summary     = out_nb_summary,
  km_fit             = out_km_fit,
  km_logrank         = out_km_logrank,
  cox_summary        = out_cox_summary,
  cox_zph            = out_cox_zph,
  lmm_summary        = out_lmm_summary,
  glmm_summary       = out_glmm_summary,
  ps_summary         = out_ps_summary,
  ps_outcome_summary = out_ps_outcome_summary,
  iptw_summary       = out_iptw_summary,
  gam_summary        = out_gam_summary,
  lasso_coefs        = out_lasso_coefs,
  lasso_predict      = out_lasso_predict,
  tree_printcp       = out_tree_printcp,
  rf_print           = out_rf_print,
  rf_importance      = out_rf_importance,
  xgb_importance     = out_xgb_importance,
  tidy_predictions   = out_tidy_predictions,
  tidy_metrics       = out_tidy_metrics,
  tidy_cv_metrics    = out_tidy_cv_metrics
)

# Browse everything at once:
# View(results)
# names(results)
# results$cox_summary