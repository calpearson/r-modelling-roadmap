# 08_tree_based_models

# ============================================================
# 8. TREE-BASED METHODS — outcome: event_f
# ============================================================
library(rpart)
library(rpart.plot)

model_tree <- rpart(event_f ~ age + bmi + comorbidity_score + sex + treatment,
                    data = cohort, method = "class")
out_tree_printcp <- printcp(model_tree)
out_tree_printcp
# rpart.plot(model_tree)  # visual, not an object

library(randomForest)
model_rf <- randomForest(event_f ~ age + bmi + comorbidity_score + sex + treatment,
                         data = cohort, ntree = 500, importance = TRUE)
out_rf_print      <- model_rf
out_rf_importance <- importance(model_rf)
out_rf_print
out_rf_importance
# varImpPlot(model_rf)  # visual, not an object

library(xgboost)
x_xgb <- model.matrix(event_f ~ age + bmi + comorbidity_score + sex + treatment,
                      data = cohort)[, -1]
y_xgb <- cohort$event
dtrain <- xgb.DMatrix(data = x_xgb, label = y_xgb)
model_xgb <- xgboost(data = dtrain, nrounds = 100, objective = "binary:logistic",
                     max_depth = 4, eta = 0.1, verbose = 0)
out_xgb_importance <- xgb.importance(model = model_xgb)
out_xgb_importance