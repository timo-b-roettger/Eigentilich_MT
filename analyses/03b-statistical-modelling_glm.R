# Load packages
library(brms)
library(tidyverse)
library(tidybayes)
library(rstudioapi)

# set the current working directory to the one where this file is
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

# Import relevant data
results_cleaned <- read_csv("../data/results-clean-06-01-23.csv")
clusters_extrapolated <- read_csv("../data/clusters-06-01-23.csv")
xDiffExplore <- read_csv("../data/x-difference-03-12-24.csv")

###########################################################
## Logistic model (bias in vertical mouse position)      ##
###########################################################

# Prepare data set for model fitting
xDiffExplore$submission_id <- as.factor(x = xDiffExplore$submission_id)
xDiffExplore$DP <- as.factor(x = xDiffExplore$DP)
xDiffExplore$group <- as.factor(x = xDiffExplore$group)
xDiffExplore$item_id <- as.factor(x = xDiffExplore$item_id)
xDiffExplore$xpos_side_num <- as.factor(x = xDiffExplore$xpos_side_num)

xtabs(~ xpos_side + DP + group, xDiffExplore)



# weakly informative priors
priors <- c(prior(cauchy(0, 0.5), class = Intercept),
            prior(cauchy(0, 0.5), class = b),
            prior(cauchy(0, 0.5), class = sd),
            prior(lkj(2), class = cor))

# Fit models (one for each relevant cluster)
MT_xdiff_glm <- brm(#xpos_side_num ~ DP * group +
                    xpos_side_num ~ DP * group * block +
                    #(1 + DP | submission_id) +
                    (1 + DP * block | submission_id) +
                    #(1 + DP * group | item_id),
                    (1 + DP * group * block | item_id),
                  data = xDiffExplore,
                  family = bernoulli, 
                  seed = 1234,
                  iter = 2000,
                  chains = 4,
                  cores = 4,
                  backend = "cmdstanr",
                  file  = "MT_xdiff_glm.RDS")

## extract posteriors
posteriors <- 
  MT_xdiff_glm |> 
  spread_draws(b_Intercept,
               b_DPindeed,
               b_groupunreliable,
               b_block,
               `b_DPindeed:groupunreliable`,
               `b_DPindeed:block`,
               `b_groupunreliable:block`,
               `b_DPindeed:groupunreliable:block`
               ) |>  
  # map onto each other
  mutate(# at 0 blocks
         actually_reliable = plogis(b_Intercept),
         indeed_reliable = plogis(b_Intercept + b_DPindeed),
         actually_unreliable = plogis(b_Intercept + b_groupunreliable),
         indeed_unreliable = plogis(b_Intercept + b_DPindeed + 
                                    b_groupunreliable + `b_DPindeed:groupunreliable`),
         # block 1
         actually_reliable_1 = plogis(b_Intercept +
                                      b_block * 1),
         indeed_reliable_1 = plogis(b_Intercept + b_DPindeed  +
                                    b_block * 1 + `b_DPindeed:block` * 1),
         actually_unreliable_1 = plogis(b_Intercept + b_groupunreliable +
                                        b_block * 1 + `b_groupunreliable:block` * 1),
         indeed_unreliable_1 = plogis(b_Intercept + b_DPindeed + 
                                      b_groupunreliable + `b_DPindeed:groupunreliable` +
                                      b_block * 1 + `b_DPindeed:block` * 1  +
                                      `b_groupunreliable:block` * 1 + `b_DPindeed:groupunreliable:block` * 1),
         # block 10
         actually_reliable_10 = plogis(b_Intercept +
                                        b_block * 10),
         indeed_reliable_10 = plogis(b_Intercept + b_DPindeed  +
                                      b_block * 10 + `b_DPindeed:block` * 10),
         actually_unreliable_10 = plogis(b_Intercept + b_groupunreliable +
                                          b_block * 10 + `b_groupunreliable:block` * 10),
         indeed_unreliable_10 = plogis(b_Intercept + b_DPindeed + 
                                        b_groupunreliable + `b_DPindeed:groupunreliable` +
                                        b_block * 10 + `b_DPindeed:block` * 10  +
                                        `b_groupunreliable:block` * 10 + `b_DPindeed:groupunreliable:block` * 10)) |>
  # make into usable data table
  select(actually_reliable_1, indeed_reliable_1, actually_unreliable_1, indeed_unreliable_1,
         actually_reliable_10, indeed_reliable_10, actually_unreliable_10, indeed_unreliable_10) |> 
  pivot_longer(cols = c(actually_reliable_1, indeed_reliable_1, actually_unreliable_1, indeed_unreliable_1,
                        actually_reliable_10, indeed_reliable_10, actually_unreliable_10, indeed_unreliable_10),
               names_to = "condition") |> 
  separate_wider_delim(condition, "_", names = c("DP", "group", "block")) |> 
  group_by(DP, group, block)  |> 
  summarise(pred_m = mean(value, na.rm = TRUE),
            pred_0025 = quantile(value, prob = 0.025),
            pred_0975 = quantile(value, prob = 0.975))
  
## Plot results
  ggplot(posteriors,
         aes(x = block, 
             y = pred_m,
             fill = DP)) +
    geom_errorbar(aes(ymin = pred_0025,
                      ymax = pred_0975),
                  color = "darkgrey",
                  width = .1) + 
    geom_line(aes(group = 1)) +
    geom_point(pch = 21,
               size = 3,
               color = "black") +
    facet_grid(DP~group) + 
    labs(y = "posterior proportion\n",
         x = "\nBlock") +
    theme_minimal() +
    theme(legend.position = "none")
 



