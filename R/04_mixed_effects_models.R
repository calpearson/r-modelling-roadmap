# 04_mixed_effects_models

# ============================================================
# 4. MIXED EFFECTS MODELS — outcome: sbp / hi_bp (cohort_long)
# ============================================================
library(lme4)

model_lmm <- lmer(sbp ~ visit + treatment + (1 | patient_id), data = cohort_long)
out_lmm_summary <- summary(model_lmm)
out_lmm_summary

model_glmm <- glmer(hi_bp ~ visit + treatment + (1 | patient_id),
                    data = cohort_long, family = binomial)
out_glmm_summary <- summary(model_glmm)
out_glmm_summary