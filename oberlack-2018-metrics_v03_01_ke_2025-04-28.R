##################
#
# Divide and Explain: Novel Metrics and Procedures for Archetype
# Analysis in Case-Based Sustainability Research
# 
# Exploring different metrics to study the relation between
# diagnostic (design) attributes and outcome attributes
#
# Cleaned script for re-identifying archetypes from:
# Oberlack & Eisenack (2018) Archetypical Barriers to Adapting Water
# Governance in River Basins to Climate Change, Journal of
# Institutional Economics 14 (3): 527–55.
# DOI:10.1017/S1744137417000509
#
# KE, Apr 25

working_directory <- 'C:/Data/HU-Box/2_Publications/044_AtSelPaper/Submission-1/DataPack'
setwd(working_directory)
input_file <- 'oberlack-2018-usedcodes.csv'

library(tidyverse)
library(fcaR) # package needs to be installed
source('atsel_v07_04_ke_2025-04-30.R')

#############
# read and prepare data

{
df <- read.csv2(input_file, 
                   header = TRUE, row.names = 1,
                   stringsAsFactors = FALSE,
                   na.strings=c('', 'NA'))

n_cases <- nrow(df)

df_d <- df          # the copy df_d will be used for the FCA
rownames(df_d) <- seq(1:n_cases) # simplify row names to numbers

# set of all models, helpful for technical reasons:
out_o <- rep(1, n_cases) %>% setNames(rownames(df_d)) %>% as_Set()

# straightforward parameter:
n_out <- out_o$cardinal() 
# (for Oberlack & Eisenack (2018), n_out = n_cases by design)

# prepare filtering of concepts which are present in at least
# two different papers:
df <- df %>% rownames_to_column('paper')
df$paper <- lapply(df$paper,
                     FUN = function(x) {unlist(strsplit(x, '\\.'))[1]}) %>%
              unlist()
df$paper <- unlist(df$paper)
} # END data preparation

#############
# START CS1 (Formal Concept Analysis)

{
# run the FCA
context_d <- FormalContext$new(df_d)
context_d$find_concepts()
n_concepts <- context_d$concepts$size()
conceptlist <- context_d$concepts

# compute simple metrics
size <- Matrix::colSums(conceptlist$extents()) 
coverage <- size / n_out
richness <- Matrix::colSums(conceptlist$intents()) 

# compute lift 
pc <- size / n_cases # pc[i] rel. frequency of concept i
pa <- colSums(df_d) / n_cases # rel frequency of single attributes
M <- conceptlist$intents() == 1
lift <- sapply(seq(1, n_concepts),
                 FUN = function(i) {pc[i]/prod(pa[M[,i]])})

# compute consistency
# (for Oberlack & Eisenack (2018), always consistency=100% by design)
joint <- sapply(seq(1, n_concepts),
                FUN = function(i)
                  {(out_o %&% conceptlist[i]$sub(1)$get_extent())$cardinal()})
consistency <- joint / size

# Special for Oberlack & Eisenack (2018):
# Compute number of distinct papers being covered
cov_paper <- apply(conceptlist$extents(), MARGIN=2,
       FUN = function(x) {length(unique(df[which(x==1),'paper']))})

# Join all metrics into one data frame and save it

# to add info on intents
get_att_string <- function(x) capture.output(x$sub(1)$get_intent()$print())
intent_strings <- sapply(seq(1, n_concepts),
                FUN = function(i) {get_att_string(conceptlist[i])}) %>%
  unlist()
# manually remove last entries, uncritical dirt from get_att_string():
intent_strings <- intent_strings[1:n_concepts]

metrics <- data.frame(intent=intent_strings, 
                      richness=richness, size=size,
                      consistency=consistency,
                      coverage=coverage,
                      lift=lift,
                      cov_paper=cov_paper)
rn <- rownames(metrics)
metrics$concept=rn # join the index of concept in conceptlist to metrics

write.csv2(metrics, file='oberlack-2018_metrics.csv')

print('Summary of concepts:')
summary(metrics)
} # END CS 1

#############
# START CS2 (Concept Filter)

{
# Thresholds
cons_min <- 1.0       # <<==== automatically met be study design
cov_min <- 2/n_out    # <<==== orig. study threshold: in at least 2 models
rho_max <- 3          # <<==== orig. study has only ATs with 2 and 3 attributes
rho_min <- 2          # <<==== 
lift_min <- 1.2       # <<==== 
cov_paper_min <- 2    # Specific for Oberlack & Eisenack 2018: ATs should be from at least 2 different papers

co_idx <- which(metrics$consistency >= cons_min
                & metrics$coverage >= cov_min
                & metrics$richness >= rho_min
                & metrics$richness <= rho_max
                & metrics$cov_paper >= cov_paper_min
                & metrics$lift >= lift_min)
co_concepts <- conceptlist[co_idx] 

print('Summary of concepts passing the concept filter:')
co_concepts$size()
summary(metrics[co_idx,])

# be careful: the concepts in co_concepts have other numbers than those in conceptlist
# but co_idx can be used to decode: co_idx[i] contains index of concept i in conceptlist:
# co_concepts[i]  == conceptlist[co_idx[i]]
} # END CS2

############
# Start OS2 (Optimal Concept Selection)

{ 
# Selection criteria
N_min <- 2      # <<==== minimum number of archetypes for scanning
N_max <- 35     # <<==== maximum number of archetypes for scanning
jointcoverage_min <- .75 
D_N_thresh <- 2 # <<==== 

# run atsel() function in 'selection mode' to compare among different (N, rho)
res <- atsel(co_concepts,
             Nrange=N_min:N_max, rhorange=rho_min:rho_max, 
             strict=TRUE, verbose=FALSE)
compare <- res$comparison

compare <- compare %>% mutate(joint_coverage = nu / n_out) %>%
  mutate(D_N=NA)

# compute D_N:
for (i in 1:nrow(compare)) {
  if (compare[i,'N'] == N_min) next
  compare[i, 'D_N'] <- as.integer(round(
    filter(compare, N==compare[i,'N'] & rho==compare[i,'rho'])$nu -
    filter(compare, N==compare[i,'N']-1 & rho==compare[i,'rho'])$nu
  ))
}

write.csv2(compare, file = paste0('oberlack_selectionmode_v02.csv'))

# with rho=2, N=22, jointcoverage_min and D_N_thresh are reached
# so, run atsel() in 'inspection mode' for just these parameters

res_22_2 <- atsel(co_concepts, Nrange=22:22, rhorange=2:2, 
                  strict=TRUE, verbose=FALSE)

nrow(res_22_2$optima) # number of optima
length(unique(c(res_22_2$optima))) # number of concepts across all optima
}

############
# START CS3 (Theoretical analysis)
#
# This qualitative task is actually not for a script. Here, we compare the results
# to those from the original Oberlack & Eisenack (2018) paper.

{
# indices(in conceptlist) of the concepts in the original paper (with rho = 2)
# (collected by hand)
orig_2 <- c(89, 81, 82, 175, 95, 66, 83, 54, 28, 52, 186, 105, 98, 183, 123, 50,
            156, 109, 59, 154, 191, 33, 176, 196, 194, 187, 178, 147, 32, 143, 142)
orig_2 <- sort(orig_2)
sum(apply(conceptlist[orig_2]$extents(), 1, max)) # check: 93 models need to be covered
conceptlist[orig_2] # print those concepts for further check with original paper

# for each optimal selection: how many of the 22 concepts are in the original study?
# (since the indices of concepts in res_22_2$optima comes from co_concepts,
# the indices need to convert to those from conceptlist)

intersect_orig <- apply(res_22_2$optima, MARGIN = 1,
                        FUN=function(row) {length(intersect(orig_2, co_idx[row]))})
intersect_orig
max(intersect_orig)
best_overlap <- which(intersect_orig == max(intersect_orig))
intersect(orig_2, co_idx[res_22_2$optima[best_overlap,]])

# for each optimal selection: how many of the 22 concepts are not in the original study?
lost_concepts <- apply(res_22_2$optima, MARGIN = 1,
                        FUN=function(row) {length(setdiff(orig_2, co_idx[row]))})
lost_concepts

# for each optimal selection: how many of the 31 original archetypes are not in the 22 concepts?
new_concepts <- apply(res_22_2$optima, MARGIN = 1,
                       FUN=function(row) {length(setdiff(co_idx[row], orig_2))})
new_concepts
}
