# 05_propensity_scores

# ============================================================
# 5. PROPENSITY SCORES — outcome: event / time+status, exposure: treatment
# ============================================================
library(MatchIt)

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