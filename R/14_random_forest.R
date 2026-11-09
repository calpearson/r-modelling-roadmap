# Research Question
# 
# Can we predict which patients will experience the event?

# Outcome == Event

library(randomForest)

out_rf <- randomForest(
  event_f ~ age +
    sex +
    bmi +
    comorbidity_score +
    sbp,
  data = cohort,
  importance = TRUE
)

# The questions are:
#   
#   How good is the model?
#   Which variables are most important?
#   How well does it classify patients?
#   How do predicted risks differ across patients?


# Base plots
randomForest::varImpPlot(out_rf)

# Variable importance test
importance_tbl <- importance(out_rf) |>
  as.data.frame() |>
  rownames_to_column("variable") |>
  arrange(desc(MeanDecreaseGini));importance_tbl

# Confusion matrix
pred_rf <- predict(out_rf)

table(
  Actual = cohort$event_f,
  Predicted = pred_rf
)

# ROC curve
library(pROC)

pred_prob <- predict(
  out_rf,
  type = "prob"
)[,"Yes"]

roc_obj <- roc(
  cohort$event_f,
  pred_prob
)

plot(roc_obj)
auc(roc_obj)

# Predicted Risk Distribution
cohort$predicted_risk <- pred_prob

ggplot(
  cohort,
  aes(predicted_risk)
) +
  geom_histogram(
    bins = 30,
    fill = "steelblue"
  )


# Risk Deciles Table
cohort |>
  mutate(
    risk_decile = ntile(predicted_risk, 10)
  ) |>
  group_by(risk_decile) |>
  summarise(
    observed_event_rate = mean(event)
  )

# Calibration plot
cohort |>
  mutate(
    pred_bin = ntile(predicted_risk, 10)
  ) |>
  group_by(pred_bin) |>
  summarise(
    observed = mean(event),
    predicted = mean(predicted_risk)
  )

# Partial Dependence Plot
library(pdp)

partial(
  out_rf,
  pred.var = "age",
  train = cohort
) |>
  autoplot()


