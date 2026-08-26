# 09_tidymodels

# ============================================================
# 9. TIDYMODELS — outcome: event_f (same predictors as section 8)
# ============================================================
library(tidymodels)

split <- initial_split(cohort, prop = 0.8, strata = event_f)
train <- training(split)
test  <- testing(split)

rf_spec <- rand_forest(trees = 500) %>%
  set_engine("ranger") %>%
  set_mode("classification")

rf_wf <- workflow() %>%
  add_formula(event_f ~ age + bmi + comorbidity_score + sex + treatment) %>%
  add_model(rf_spec)

model_tidy_rf <- rf_wf %>% fit(data = train)

out_tidy_predictions <- predict(model_tidy_rf, test) %>% bind_cols(test)
out_tidy_metrics     <- out_tidy_predictions %>% metrics(truth = event_f, estimate = .pred_class)
out_tidy_predictions
out_tidy_metrics

folds <- vfold_cv(train, v = 5)
out_tidy_cv_metrics <- fit_resamples(rf_wf, folds) %>% collect_metrics()
out_tidy_cv_metrics