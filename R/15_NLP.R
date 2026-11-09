# Research Question
# 
# Can we extract clinical information from free text?
  
library(tidytext)

out_terms <- cohort %>%
  unnest_tokens(word, clinical_note) %>%
  count(word, sort = TRUE)
