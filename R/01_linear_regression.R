# 01_linear_regression

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
