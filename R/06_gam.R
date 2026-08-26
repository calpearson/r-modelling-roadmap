# 06_gam

# ============================================================
# 6. GAMs — outcome: sbp, smooth term on age
# ============================================================
library(mgcv)

model_gam <- gam(sbp ~ s(age) + bmi + treatment, data = cohort, family = gaussian())
out_gam_summary <- summary(model_gam)
out_gam_summary
# plot(model_gam, pages = 1, shade = TRUE)  # visual, not an object
