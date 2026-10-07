explain_study_design <- function(
    design = c("mean_difference",
               "proportion_difference",
               "cohort_rr",
               "case_control_or"),
    alpha = 0.05,
    power = 0.80,
    effect_size = NULL,
    p1 = NULL,
    p2 = NULL,
    baseline_risk = NULL,
    rr = NULL,
    or = NULL,
    attrition = 0,
    show_table = TRUE,
    show_plot = TRUE,
    show_sap_text = TRUE
){
  
  suppressPackageStartupMessages({
    library(tidyverse)
    library(pwr)
  })
  
  design <- match.arg(design)
  
  assumptions <- character()
  
  sensitivity <- tibble()
  
  references <- c(
    "Cohen J. Statistical Power Analysis for the Behavioral Sciences. 2nd ed. 1988.",
    "Chow SC, Shao J, Wang H. Sample Size Calculations in Clinical Research. 3rd ed. 2017.",
    "Rothman KJ, Greenland S, Lash TL. Modern Epidemiology. 4th ed."
  )
  
  # ------------------------------------------------------------------------
  # Mean Difference
  # ------------------------------------------------------------------------
  
  if(design == "mean_difference"){
    
    if(is.null(effect_size)){
      stop("effect_size required.")
    }
    
    ss <- pwr.t.test(
      d = effect_size,
      power = power,
      sig.level = alpha,
      type = "two.sample"
    )
    
    n_per_group <- ceiling(ss$n)
    total_n <- n_per_group * 2
    
    assumptions <- c(
      paste0("Effect size (Cohen's d) = ", effect_size),
      "Outcome approximately normally distributed",
      "Independent groups",
      "Equal variance between groups",
      paste0("Power = ", power * 100, "%"),
      paste0("Alpha = ", alpha)
    )
    
    sensitivity <- tibble(
      effect_size = seq(
        max(0.1, effect_size * 0.5),
        effect_size * 1.5,
        length.out = 8
      )
    ) %>%
      mutate(
        n_per_group = map_dbl(
          effect_size,
          ~ ceiling(
            pwr.t.test(
              d = .x,
              power = power,
              sig.level = alpha,
              type = "two.sample"
            )$n
          )
        )
      )
    
  } else if(design == "proportion_difference"){
    
    # ----------------------------------------------------------------------
    # Difference in Proportions
    # ----------------------------------------------------------------------
    
    if(is.null(p1) | is.null(p2)){
      stop("p1 and p2 required.")
    }
    
    h <- ES.h(p1, p2)
    
    ss <- pwr.2p.test(
      h = h,
      power = power,
      sig.level = alpha
    )
    
    n_per_group <- ceiling(ss$n)
    total_n <- n_per_group * 2
    
    assumptions <- c(
      paste0("Group 1 risk = ", scales::percent(p1)),
      paste0("Group 2 risk = ", scales::percent(p2)),
      "Independent observations",
      paste0("Power = ", power * 100, "%"),
      paste0("Alpha = ", alpha)
    )
    
    sensitivity <- tibble(
      p2 = seq(
        max(0.001, p2 * 0.5),
        min(0.999, p2 * 1.5),
        length.out = 8
      )
    ) %>%
      mutate(
        h = map_dbl(p2, ~ ES.h(p1, .x)),
        n_per_group = map_dbl(
          h,
          ~ tryCatch(
            ceiling(
              pwr.2p.test(
                h = .x,
                power = power,
                sig.level = alpha
              )$n
            ),
            error = function(e) NA_real_
          )
        )
      )
    
  } else if(design == "cohort_rr"){
    
    # ----------------------------------------------------------------------
    # Cohort Study
    # ----------------------------------------------------------------------
    
    if(is.null(baseline_risk) | is.null(rr)){
      stop("baseline_risk and rr required.")
    }
    
    p1 <- baseline_risk
    p2 <- baseline_risk * rr
    
    if(p2 > 1){
      stop("RR produces risk > 1")
    }
    
    h <- ES.h(p1, p2)
    
    ss <- pwr.2p.test(
      h = h,
      power = power,
      sig.level = alpha
    )
    
    n_per_group <- ceiling(ss$n)
    total_n <- n_per_group * 2
    
    assumptions <- c(
      paste0("Baseline risk = ", scales::percent(p1)),
      paste0("Expected RR = ", rr),
      "Minimal exposure misclassification",
      "Complete follow-up",
      paste0("Power = ", power * 100, "%"),
      paste0("Alpha = ", alpha)
    )
    
    sensitivity <- tibble(
      RR = seq(
        max(1.05, rr * 0.5),
        rr * 1.5,
        length.out = 8
      )
    ) %>%
      mutate(
        risk2 = pmin(0.999, p1 * RR),
        h = map_dbl(risk2, ~ ES.h(p1, .x)),
        n_per_group = map_dbl(
          h,
          ~ {
            if(abs(.x) < 1e-8){
              NA_real_
            } else {
              tryCatch(
                ceiling(
                  pwr.2p.test(
                    h = .x,
                    power = power,
                    sig.level = alpha
                  )$n
                ),
                error = function(e) NA_real_
              )
            }
          }
        )
      )
    
  } else if(design == "case_control_or"){
    
    # ----------------------------------------------------------------------
    # Case Control
    # ----------------------------------------------------------------------
    
    if(is.null(baseline_risk) | is.null(or)){
      stop("baseline_risk and or required.")
    }
    
    p1 <- baseline_risk
    
    p2 <- (or * p1) / (1 - p1 + or * p1)
    
    h <- ES.h(p1, p2)
    
    ss <- pwr.2p.test(
      h = h,
      power = power,
      sig.level = alpha
    )
    
    n_per_group <- ceiling(ss$n)
    total_n <- n_per_group * 2
    
    assumptions <- c(
      paste0("Baseline exposure prevalence = ", scales::percent(p1)),
      paste0("Expected OR = ", or),
      "Independent sampling",
      "Minimal selection bias",
      paste0("Power = ", power * 100, "%"),
      paste0("Alpha = ", alpha)
    )
    
    sensitivity <- tibble(
      OR = seq(
        max(1.05, or * 0.5),
        or * 1.5,
        length.out = 8
      )
    ) %>%
      mutate(
        exposure2 = (OR * p1) / (1 - p1 + OR * p1),
        h = map_dbl(exposure2, ~ ES.h(p1, .x)),
        n_per_group = map_dbl(
          h,
          ~ {
            if(abs(.x) < 1e-8){
              NA_real_
            } else {
              tryCatch(
                ceiling(
                  pwr.2p.test(
                    h = .x,
                    power = power,
                    sig.level = alpha
                  )$n
                ),
                error = function(e) NA_real_
              )
            }
          }
        )
      )
    
  }
  
  # ------------------------------------------------------------------------
  # Attrition Adjustment
  # ------------------------------------------------------------------------
  
  adjusted_total_n <- ceiling(
    total_n / (1 - attrition)
  )
  
  # ------------------------------------------------------------------------
  # Interpretation
  # ------------------------------------------------------------------------
  
  interpretation <- paste0(
    "The study is designed using a two-sided significance level of ",
    alpha,
    " and statistical power of ",
    power * 100,
    "%.\n\n",
    
    "Based on the assumptions provided, approximately ",
    format(total_n, big.mark = ","),
    " participants are required.",
    
    if(attrition > 0){
      paste0(
        " Assuming an attrition rate of ",
        scales::percent(attrition),
        ", recruitment should be increased to approximately ",
        format(adjusted_total_n, big.mark = ","),
        " participants."
      )
    } else {
      ""
    },
    
    "\n\nThis estimate assumes all model assumptions are correct, the specified effect size is realistic, data quality is high, and missing data are minimal.",
    
    "\n\nIf the true effect is smaller than anticipated, a substantially larger sample may be required. Investigators should examine prior studies and conduct sensitivity analyses to assess the robustness of assumptions.",
    
    "\n\nReferences: ",
    paste(references, collapse = " | ")
  )
  
  sap_text <- NULL
  
  if(show_sap_text){
    
    sap_text <- paste0(
      "Sample size calculations were performed assuming a two-sided alpha of ",
      alpha,
      " and ",
      power * 100,
      "% statistical power. A total of ",
      format(total_n, big.mark = ","),
      " participants were estimated to be required."
    )
    
  }
  
  summary_table <- tibble(
    Design = design,
    Alpha = alpha,
    Power = power,
    Calculated_N = total_n,
    Attrition_Adjusted_N = adjusted_total_n
  )
  
  plot_obj <- NULL
  
  if(show_plot && nrow(sensitivity) > 0){
    
    x_var <- names(sensitivity)[1]
    
    plot_obj <- ggplot(
      sensitivity,
      aes(
        x = .data[[x_var]],
        y = n_per_group
      )
    ) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 3) +
      theme_bw() +
      labs(
        title = paste("Sensitivity Analysis:", design),
        x = x_var,
        y = "Required Participants Per Group"
      )
    
  }
  
  output <- list(
    interpretation = interpretation,
    assumptions = assumptions,
    results = summary_table,
    sensitivity = sensitivity,
    references = references
  )
  
  if(show_table){
    output$table <- summary_table
  }
  
  if(show_plot){
    output$plot <- plot_obj
  }
  
  if(show_sap_text){
    output$sap_text <- sap_text
  }
  
  return(output)
  
}





# Cohort Example
result <- explain_study_design(
  design = "cohort_rr",
  baseline_risk = 0.05,
  rr = 2,
  power = 0.8,
  alpha = 0.05,
  attrition = 0.10,
  show_plot = TRUE,
  show_table = TRUE
)

cat(result$interpretation)
result$assumptions
result$table
result$plot
cat(result$sap_text)
result$references



result <- explain_study_design(
  design = "mean_difference",
  effect_size = 0.5,
  power = 0.8,
  alpha = 0.05
)
cat(result$interpretation)
result$table
result$plot


result <- explain_study_design(
  design = "proportion_difference",
  p1 = 0.05,
  p2 = 0.10,
  power = 0.8,
  alpha = 0.05
)
result$table
result$plot


result <- explain_study_design(
  design = "cohort_rr",
  baseline_risk = 0.05,
  rr = 2,
  power = 0.8,
  alpha = 0.05,
  attrition = 0.10
)
result$table
result$plot
cat(result$interpretation)



result <- explain_study_design(
  design = "case_control_or",
  baseline_risk = 0.20,
  or = 2,
  power = 0.8,
  alpha = 0.05
)
result$table
result$plot

