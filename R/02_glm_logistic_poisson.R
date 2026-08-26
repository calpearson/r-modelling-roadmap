# 02_glm_logistic_poisson

# ============================================================
# 2. GLMs — outcome: event (logistic), n_hosp (Poisson/NegBin)
# ============================================================
model_logit <- glm(event ~ age + sex + treatment + comorbidity_score,
                   data = cohort, family = binomial(link = "logit"))

out_logit_summary <- summary(model_logit)
out_logit_or      <- exp(coef(model_logit))       # odds ratios
out_logit_or_ci   <- exp(confint(model_logit))

out_logit_summary
out_logit_or
out_logit_or_ci

model_pois <- glm(n_hosp ~ age + treatment + offset(log(person_time)),
                  data = cohort, family = poisson(link = "log"))
out_pois_summary <- summary(model_pois)
out_pois_summary

library(MASS)
model_nb <- glm.nb(n_hosp ~ age + treatment + offset(log(person_time)), data = cohort)
out_nb_summary <- summary(model_nb)
out_nb_summary