# Load packages
library(brms)
library(tidyverse)
library(tidybayes)
options(mc.cores = parallel::detectCores())

# Import relevant data
results_cleaned <- read_csv("../data/results-clean-06-01-23.csv")
clusters_extrapolated <- read_csv("../data/clusters-06-01-23.csv")
xDiffExplore <- read_csv("../data/x-difference-06-01-23.csv")

################################################################################
## Categorical model (proportion of data points per mouse trajectory cluster) ##
##                    see 02b-plotting-cluster-figure                         ##
################################################################################

# Function to extrapolate mouse trajectories (see 01-data-processing-and-wrangling)
extrapolate_trajectories <- function(mt_time, mt_x_norm, mt_y_norm) {
  
  t_star <- seq(0, 4000, by = 10)
  
  n_recording <- length(mt_time)
  
  mt_x_star <- rep(NA, length(t_star))
  mt_y_star <- rep(NA, length(t_star))
  
  for (i in 1:length(t_star)) {
    if (t_star[i] <= min(mt_time)) {
      # for points prior to recording MT signal
      mt_x_star[i] <- mt_x_norm[1]
      mt_y_star[i] <- mt_y_norm[1]
    } else if (t_star[i] >= max(mt_time)) {
      # for points after recording stopped MT signal
      mt_x_star[i] <- mt_x_norm[n_recording]
      mt_y_star[i] <- mt_y_norm[n_recording]
    } else if (t_star[i] %in% mt_time) {
      mt_x_star[i] <- mt_x_norm[which(t_star[i]==mt_time)]
      mt_y_star[i] <- mt_y_norm[which(t_star[i]==mt_time)]
    } else {
      # for points where we have data for MT signal
      # get the closet two times for which we have recording
      upper_index <- which(mt_time>t_star[i])[1]
      mt_time[upper_index]
      lower_index <- max(which(mt_time<t_star[i]))
      mt_time[lower_index]
      mt_x_star[i] <- mt_x_norm[lower_index] +
        (mt_x_norm[upper_index] - mt_x_norm[lower_index]) *
        (t_star[i] - mt_time[lower_index]) /
        (mt_time[upper_index] - mt_time[lower_index])
      mt_y_star[i] <- mt_y_norm[lower_index] +
        (mt_y_norm[upper_index] - mt_y_norm[lower_index]) *
        (t_star[i] - mt_time[lower_index]) /
        (mt_time[upper_index] - mt_time[lower_index])
      # x^* = (x_2-x_1) * (t^* - t_1)/(t_2 - t_1)
    }
  }
  out = tibble(
    t_star,
    mt_x_star,
    mt_y_star
  )
  return( out )
}

# Create data frame 'results_stats' 
results_stats <- results_cleaned %>%
  full_join(clusters_extrapolated, by = "mt_id") %>%
  group_by(mt_id, submission_id, DP, group, item_id, cluster, condition) %>%
  group_modify(function(d,...) {
    extrapolate_trajectories(d$mt_time, d$mt_x_norm, d$mt_y_norm)}) %>%
  ungroup() %>% 
  select(submission_id, DP, group, condition, item_id, cluster) %>% 
  distinct()

# Prepare data set for model fitting
# results_stats$submission_id <- as.factor(x = results_stats$submission_id)
# results_stats$DP <- as.factor(x = results_stats$DP)
# results_stats$group <- as.factor(x = results_stats$group)
# results_stats$item_id <- as.factor(x = results_stats$item_id)

# Fit model
# MT_cluster_cat <- brm(cluster ~ DP * group +
#                   (1 + DP | submission_id) +
#                   (1 + DP * group | item_id),
#                   data = results_stats,
#                   family = categorical(),
#                   file = "MT_cluster_cat.RDS")

# Load model object (fitted model)
MT_cluster_cat <- readRDS("MT_cluster_cat.RDS")

# Extract estimated proportions from fitted model
new_data <- results_stats
new_data$cluster_size <- 1

prob_categorical <- fitted(MT_cluster_cat, newdata = new_data) %>%
  as_tibble() %>% 
  select(contains(c("Estimate", "Q2.5", "Q97.5"))) %>%
  set_names(str_c("cluster_", 1:3),
            str_c("lower_", 1:3),
            str_c("upper_", 1:3)) %>%
  bind_cols(new_data) %>%
  group_by(DP, group) %>% 
  summarize(prob_1 = round(mean(cluster_1), digits = 2),
            prob_2 = round(mean(cluster_2), digits = 2),
            prob_3 = round(mean(cluster_3), digits = 2),
            CrI_upper_1 = round(mean(upper_1), digits = 2),
            CrI_upper_2 = round(mean(upper_2), digits = 2),
            CrI_upper_3 = round(mean(upper_3), digits = 2),
            CrI_lower_1 = round(mean(lower_1), digits = 2),
            CrI_lower_2 = round(mean(lower_2), digits = 2),
            CrI_lower_3 = round(mean(lower_3), digits = 2))

###########################################################
## Linear model (differences in vertical mouse position) ##
##        see 02c-plotting-x-difference-figure           ##
###########################################################

# Prepare data set for model fitting
# xDiffExplore$submission_id <- as.factor(x = xDiffExplore$submission_id)
# xDiffExplore$DP <- as.factor(x = xDiffExplore$DP)
# xDiffExplore$group <- as.factor(x = xDiffExplore$group)
# xDiffExplore$item_id <- as.factor(x = xDiffExplore$item_id)

# Fit models (one for each relevant cluster)
# MT_xdiff_2 <- brm(x_DPDiff ~ DP * group +
#                     (1 + DP | submission_id) +
#                     (1 + DP * group | item_id),
#                   data = filter(xDiffExplore, cluster == "2"),
#                   control = list(adapt_delta = 0.99),
#                   iter = 3000,
#                   file = "MT_xdiff_2.RDS")

# MT_xdiff_2 <- brm(x_DPDiff ~ DP * group +
#                     (1 + DP | submission_id) +
#                     (1 + DP * group | item_id),
#                   data = filter(xDiffExplore, cluster == "3"),
#                   control = list(adapt_delta = 0.99),
#                   iter = 3000,
#                   file = "MT_xdiff_3.RDS")

# Load model objects (fitted models)
MT_xdiff_2 <- readRDS("MT_xdiff_2.RDS")

MT_xdiff_3 <- readRDS("MT_xdiff_3.RDS")

# Compute empirical differences
means_xdiff <- xDiffExplore %>% 
  group_by(group, cluster, DP) %>% 
  summarize(mean = round(mean(x_DPDiff), digits = 2))

# Create vector with hypotheses to be tested
hypotheses <- 
  c("Intercept > 0",
    "Intercept + DPindeed > 0",
    "Intercept + groupunreliable > 0",
    "Intercept + DPindeed + groupunreliable + DPindeed:groupunreliable > 0")

# Cluster 2
# Extract estimated differences from fitted model
prob_linear_cluster_2 <- brms::hypothesis(MT_xdiff_2, hypotheses)
 
prob_linear_cluster_2 <- bind_rows(prob_linear_cluster_2$hypothesis) %>%
  mutate(Hypothesis = recode(Hypothesis,
                             "(Intercept) > 0" = "[reliable] Actually > 0",
                             "(Intercept+DPindeed) > 0" = "[reliable] Indeed > 0",
                             "(Intercept+groupunreliable) > 0" = "[unreliable] Actually > 0",
                             "(Intercept+DPindeed+groupunreliable+DPindeed:groupunreliable) > 0" = "[unreliable] Indeed > 0")) %>%
  reframe(across(c(2:7), ~ round(.x, digits = 2)))

# Plot posterior probability distributions
MT_xdiff_2 %>%
  spread_draws(b_Intercept,
               b_DPindeed, 
               b_groupunreliable,
               `b_DPindeed:groupunreliable`) %>%
  mutate(b_DPindeed_new = b_Intercept + b_DPindeed,
         b_groupunreliable_new = b_Intercept + b_groupunreliable,
         `b_DPindeed:groupunreliable_new` = b_Intercept +  b_DPindeed + b_groupunreliable + `b_DPindeed:groupunreliable`) %>%
  gather_draws(b_Intercept,
               b_DPindeed_new,
               b_groupunreliable_new,
               `b_DPindeed:groupunreliable_new`) %>%
  mutate(.variable = factor(x = .variable,
                            levels = c("b_DPindeed:groupunreliable_new",
                                       "b_groupunreliable_new",
                                       "b_DPindeed_new",
                                       "b_Intercept"))) %>%
  ggplot(aes(y = .variable, x = .value, fill = .variable, alpha = .variable)) +
  stat_halfeye(.width = c(.90, .95), point_interval = "median_qi") +
  scale_y_discrete(labels = c("Tatsächlich\n(unreliable)",
                              "Eigentlich\n(unreliable)",
                              "Tatsächlich\n(reliable)",
                              "Eigentlich\n(reliable)")) +
  scale_fill_manual(breaks = c("b_DPindeed:groupunreliable_new",
                               "b_groupunreliable_new",
                               "b_DPindeed_new",
                               "b_Intercept"),
                    values = c("#81B29A",
                               "#E07A5F",
                               "#81B29A",
                               "#E07A5F")) +
  scale_alpha_manual(breaks = c("b_DPindeed:groupunreliable_new",
                                "b_groupunreliable_new",
                                "b_DPindeed_new",
                                "b_Intercept"),
                     values = c("b_DPindeed:groupunreliable_new" = .3,
                                "b_groupunreliable_new" = .3,
                                "b_DPindeed_new" = 1,
                                "b_Intercept" = 1)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold",
                                  size = 20, hjust = .5,
                                  margin = margin(t = 0, r = 0, b = 15, l = 0)),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 14, color = "black",
                                   margin = margin(t = 10, r = 0, b = 15, l = 0)),
        axis.text.y = element_text(face = "bold", size = 16, color = "black", 
                                   hjust = .5,
                                   margin = margin(t = 0, r = 10, b = 0, l = 0)))

# Cluster 3
# Extract estimated differences from fitted model
prob_linear_cluster_3 <- brms::hypothesis(MT_xdiff_3, hypotheses)

prob_linear_cluster_3 <- bind_rows(prob_linear_cluster_3$hypothesis) %>%
  mutate(Hypothesis = recode(Hypothesis,
                             "(Intercept) > 0" = "[reliable] Actually > 0",
                             "(Intercept+DPindeed) > 0" = "[reliable] Indeed > 0",
                             "(Intercept+groupunreliable) > 0" = "[unreliable] Actually > 0",
                             "(Intercept+DPindeed+groupunreliable+DPindeed:groupunreliable) > 0" = "[unreliable] Indeed > 0")) %>%
  reframe(across(c(2:7), ~ round(.x, digits = 2)))

# Plot posterior probability distributions
MT_xdiff_3 %>%
  spread_draws(b_Intercept,
               b_DPindeed, 
               b_groupunreliable,
               `b_DPindeed:groupunreliable`) %>%
  mutate(b_DPindeed_new = b_Intercept + b_DPindeed,
         b_groupunreliable_new = b_Intercept + b_groupunreliable,
         `b_DPindeed:groupunreliable_new` = b_Intercept +  b_DPindeed + b_groupunreliable + `b_DPindeed:groupunreliable`) %>%
  gather_draws(b_Intercept,
               b_DPindeed_new,
               b_groupunreliable_new,
               `b_DPindeed:groupunreliable_new`) %>%
  mutate(.variable = factor(x = .variable,
                            levels = c("b_DPindeed:groupunreliable_new",
                                       "b_groupunreliable_new",
                                       "b_DPindeed_new",
                                       "b_Intercept"))) %>%
  ggplot(aes(y = .variable, x = .value, fill = .variable, alpha = .variable)) +
  stat_halfeye(.width = c(.90, .95), point_interval = "median_qi") +
  scale_y_discrete(labels = c("Tatsächlich\n(unreliable)",
                              "Eigentlich\n(unreliable)",
                              "Tatsächlich\n(reliable)",
                              "Eigentlich\n(reliable)")) +
  scale_fill_manual(breaks = c("b_DPindeed:groupunreliable_new",
                               "b_groupunreliable_new",
                               "b_DPindeed_new",
                               "b_Intercept"),
                    values = c("#81B29A",
                               "#E07A5F",
                               "#81B29A",
                               "#E07A5F")) +
  scale_alpha_manual(breaks = c("b_DPindeed:groupunreliable_new",
                                "b_groupunreliable_new",
                                "b_DPindeed_new",
                                "b_Intercept"),
                     values = c("b_DPindeed:groupunreliable_new" = .3,
                                "b_groupunreliable_new" = .3,
                                "b_DPindeed_new" = 1,
                                "b_Intercept" = 1)) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold",
                                  size = 20, hjust = .5,
                                  margin = margin(t = 0, r = 0, b = 15, l = 0)),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 14, color = "black",
                                   margin = margin(t = 10, r = 0, b = 15, l = 0)),
        axis.text.y = element_text(face = "bold", size = 16, color = "black", 
                                   hjust = .5,
                                   margin = margin(t = 0, r = 10, b = 0, l = 0)))
