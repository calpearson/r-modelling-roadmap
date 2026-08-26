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
