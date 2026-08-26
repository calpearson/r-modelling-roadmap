# 00_build_cohort

# ============================================================
# R MODELLING ROADMAP — ONE DATASET, NAMED OUTPUT OBJECTS
# Every modelling step below assigns its result to a clearly
# named object (out_*) so you can inspect it after running,
# e.g. out_lm, out_logit, out_cox, out_rf, etc.
# A summary list `results` is built at the end.
# ============================================================

# ---- Packages you'll need (install once) ----
# install.packages(c("survival", "lme4", "MatchIt", "WeightIt",
#                     "mgcv", "glmnet", "rpart", "rpart.plot",
#                     "randomForest", "xgboost", "tidymodels", "MASS"))
# keras3 is optional (section 10) — install separately if you want it.


# ============================================================
# 0. BUILD THE COHORT (run this once, first)
# ============================================================
set.seed(42)
n <- 1000

cohort <- data.frame(
  patient_id        = 1:n,
  age               = round(rnorm(n, 60, 12)),
  sex               = factor(sample(c("M", "F"), n, replace = TRUE)),
  bmi               = round(rnorm(n, 28, 5), 1),
  comorbidity_score = pmax(0, round(rnorm(n, 5, 2)))
)

ps_true <- plogis(-2 + 0.03 * cohort$age + 0.15 * cohort$comorbidity_score)
cohort$treatment <- factor(rbinom(n, 1, ps_true), labels = c("Control", "Drug"))

cohort$sbp <- round(130 + 0.3 * cohort$age + 1.2 * cohort$bmi -
                      4 * (cohort$treatment == "Drug") + rnorm(n, 0, 8), 1)

lp_event <- -4 + 0.04 * cohort$age + 0.08 * cohort$comorbidity_score -
  0.5 * (cohort$treatment == "Drug")
cohort$event <- rbinom(n, 1, plogis(lp_event))
cohort$event_f <- factor(cohort$event, labels = c("No", "Yes"))

cohort$person_time <- round(runif(n, 0.5, 5), 1)
rate <- 0.15 * exp(0.2 * (cohort$treatment == "Drug") + 0.03 * cohort$age)
cohort$time   <- pmin(round(rexp(n, rate) * 365), round(cohort$person_time * 365))
cohort$status <- as.integer(cohort$time < round(cohort$person_time * 365))

cohort$n_hosp <- rpois(n, exp(-1.5 + 0.02 * cohort$age) * cohort$person_time)

noise <- matrix(rnorm(n * 15), n, 15)
colnames(noise) <- paste0("marker", 1:15)
cohort <- cbind(cohort, noise)

# Long (repeated-measures) version of the SAME patients — for mixed models
visits <- 4
cohort_long <- cohort[rep(1:n, each = visits),
                      c("patient_id", "age", "sex", "treatment", "comorbidity_score")]
cohort_long$visit <- rep(1:visits, times = n)
patient_re <- rep(rnorm(n, 0, 6), each = visits)
cohort_long$sbp <- with(cohort_long,
                        130 + 0.3*age - 4*(treatment == "Drug") - 1.5*visit) +
  patient_re + rnorm(n * visits, 0, 5)
cohort_long$hi_bp <- rbinom(n * visits, 1,
                            plogis(-1 + 0.02*cohort_long$age +
                                     0.15*cohort_long$visit -
                                     0.3*(cohort_long$treatment == "Drug") +
                                     patient_re/10))

out_cohort_head      <- head(cohort)
out_cohort_str       <- str(cohort)
out_cohort_long_head <- head(cohort_long)