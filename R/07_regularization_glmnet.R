# 07_regularization_glmnet

# ============================================================
# 7. REGULARIZATION (glmnet) — outcome: event, predictors incl. noise markers
# ============================================================
library(glmnet)

x_glmnet <- model.matrix(event ~ age + sex + bmi + comorbidity_score + treatment +
                           marker1 + marker2 + marker3 + marker4 + marker5 +
                           marker6 + marker7 + marker8 + marker9 + marker10 +
                           marker11 + marker12 + marker13 + marker14 + marker15,
                         data = cohort)[, -1]
y_glmnet <- cohort$event

model_lasso_cv  <- cv.glmnet(x_glmnet, y_glmnet, alpha = 1, family = "binomial")
out_lasso_coefs <- coef(model_lasso_cv, s = "lambda.min")
out_lasso_coefs   # true predictors should survive; markers should shrink to ~0

out_lasso_predict <- predict(model_lasso_cv, newx = x_glmnet[1:5, ],
                             s = "lambda.min", type = "response")
out_lasso_predict
