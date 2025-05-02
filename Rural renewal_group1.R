##################
# START OS1 (Grouping)
# This qualitative task is actually not for a script.
# Please refer to the main text for details.
#
#
# Divide and Explain: Novel Metrics and Procedures for Archetype
# Analysis in Case-Based Sustainability Research
# 
# Exploring different metrics to study the relation between
# diagnostic (design) attributes and outcome attributes
#
# Cleaned script for re-identifying archetypes from:
# Wang et al. (2019) Sustainable Rural Renewal in China: 
# Archetypical Patterns. Ecology and Society 24 (3).
# https://doi.org/10.5751/ES-11069-240332
#
# RY, May 1st

working_directory <- 'INSERT'
setwd(working_directory)
input_file <- 'Rural renewal data_group1.csv'

library(tidyverse)
library(fcaR)
source('atsel_v07_04_ke_2025-04-30.R')

#############
# read and prepare data

df <- read.csv(input_file, 
               header = TRUE, row.names = 1,
               stringsAsFactors = FALSE,
               na.strings=c('', 'NA'))
n_cases <- nrow(df)

# The outcome variable of the group 1 is O11
out_a <- 'O11' # <======= plug in the name of the outcome attribute here

out_o <- as.numeric(df[,out_a]) %>% setNames(rownames(df)) %>% as_Set()
n_out<- out_o$cardinal()
df_d <- select(df, !c(starts_with('O'))) 
n_configs_d <- n_cases - sum(duplicated(df_d))
# END data preparation

#############
# START CS1 (Formal Concept Analysis)

  # run the FCA
  context_d <- FormalContext$new(df_d)
  context_d$find_concepts()
  n_concepts <- context_d$concepts$size()
  conceptlist <- context_d$concepts
  # compute simple metrics
  size <- Matrix::colSums(conceptlist$extents()) 
  coverage <- size / n_out
  richness <- Matrix::colSums(conceptlist$intents()) 
  
  # compute consistency
  # (for Wang et al. (2019), always consistency=100% due to grouping)
  joint <- sapply(seq(1, n_concepts),
                  FUN = function(i)
                  {(out_o %&% conceptlist[i]$sub(1)$get_extent())$cardinal()})
  consistency <- joint / size
  
  # Join all metrics into one data frame and save it
  
  # to add info on intents
  get_att_string <- function(x, max_chars = 80) {
    out <- capture.output(x$sub(1)$get_intent()$print())
    oneline <- paste(out, collapse = " ")
    oneline <- gsub("\\s+", " ", oneline)
    substr(oneline, 1, max_chars)
  }
  intent_strings <- sapply(seq(1, n_concepts),
                           FUN = function(i) {get_att_string(conceptlist[i])}) %>%
    unlist()

  metrics <- data.frame(intent=intent_strings, 
                        richness=richness, size=size,
                        consistency=consistency,
                        coverage=coverage)
  rn <- rownames(metrics)
  metrics$concept=rn # join the index of concept in conceptlist to metrics
  
  write.csv2(metrics, file='Wang-2019_metrics_group1.csv')
  
  print('Summary of concepts:')
  summary(metrics)
# END CS 1

#############
# START CS2 (Concept Filter)

  # Thresholds
  cons_min <- 1.0       # <<==== automatically met due to grouping
  cov_min <- 2/n_out    # <<==== orig. study threshold: in at least 2 cases
  rho_max <- 21          # <<==== Concepts with 21 attributes are already quite complicated and specific, being nearly half of the total number of attributes (m = 47)
  rho_min <- 2          # <<==== minimal richness ensures a configuration of attributes rather than an individual attribute
  
  co_idx <- which(metrics$consistency >= cons_min
                  & metrics$coverage >= cov_min
                  & metrics$richness >= rho_min
                  & metrics$richness <= rho_max)
  co_concepts <- conceptlist[co_idx] 
  
  print('Summary of concepts passing the concept filter:')
  co_concepts$size()
  summary(metrics[co_idx,])
  
  # be careful: the concepts in co_concepts have other numbers than those in conceptlist
  # but co_idx can be used to decode: co_idx[i] contains index of concept i in conceptlist:
  # co_concepts[i]  == conceptlist[co_idx[i]]
# END CS2  
 

############
# START CS3 (Theoretical analysis)
#
# This qualitative task is actually not for a script.
# Please refer to the main text for details.