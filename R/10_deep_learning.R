# 10_deep_learning

# ============================================================
# 10. DEEP LEARNING (optional) — reuse x_xgb / y_xgb from section 8
# ============================================================
Requires: install.packages("keras3"); keras3::install_keras()

library(keras3)

n_train <- floor(0.8 * nrow(x_xgb))
x_train <- x_xgb[1:n_train, ]; y_train <- y_xgb[1:n_train]
x_test  <- x_xgb[(n_train+1):nrow(x_xgb), ]; y_test <- y_xgb[(n_train+1):nrow(x_xgb)]

model_keras <- keras_model_sequential() %>%
  layer_dense(units = 32, activation = "relu", input_shape = ncol(x_train)) %>%
  layer_dense(units = 16, activation = "relu") %>%
  layer_dense(units = 1, activation = "sigmoid")

model_keras %>% compile(loss = "binary_crossentropy", optimizer = "adam", metrics = "accuracy")
out_keras_fit <- model_keras %>% fit(x_train, y_train, epochs = 20, validation_split = 0.2)
out_keras_fit