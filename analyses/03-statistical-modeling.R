# Load packages
library(brms)
library(tidyverse)
library(tidybayes)

# Import relevant data
results_cleaned <- read_csv("../data/results-clean-06-01-23.csv")
xDiffExplore <- read_csv("../data/x-difference-03-12-24.csv")

######################################################
## Logistic model (bias in vertical mouse position) ##
##        see 02b-plotting-x-difference-figure      ##
######################################################

# Prepare data set for model fitting
xDiffExplore$submission_id <- as.factor(x = xDiffExplore$submission_id)
xDiffExplore$DP <- as.factor(x = xDiffExplore$DP)
xDiffExplore$group <- as.factor(x = xDiffExplore$group)
xDiffExplore$item_id <- as.factor(x = xDiffExplore$item_id)
xDiffExplore$xpos_side_num <- as.factor(x = xDiffExplore$xpos_side_num)

# Fit model
# Weakly informative priors
priors <- c(prior(cauchy(0, 0.5), class = Intercept),
            prior(cauchy(0, 0.5), class = b),
            prior(cauchy(0, 0.5), class = sd),
            prior(lkj(2), class = cor))

MT_xdiff <- brm(
  xpos_side_num ~ DP * group +
    (1 + DP | submission_id) +
    (1 + DP * group | item_id),
  data = xDiffExplore,
  family = bernoulli, 
  seed = 1234,
  iter = 2000,
  chains = 4,
  cores = 4,
  # backend = "cmdstanr",
  file  = "MT_xdiff.RDS")

# Load model object (fitted model)
# MT_xdiff <- readRDS("MT_xdiff.RDS")

# Compute empirical proportions for reference
xDiffExplore %>% 
  group_by(group, DP, xpos_side) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))

# Extract estimated probabilities from fitted model
posteriors <- MT_xdiff |> 
  spread_draws(b_Intercept,
               b_DPindeed,
               b_groupunreliable,
               `b_DPindeed:groupunreliable`) |>  
  # map onto each other
  mutate(actually_reliable = plogis(b_Intercept),
         indeed_reliable = plogis(b_Intercept + b_DPindeed),
         actually_unreliable = plogis(b_Intercept + b_groupunreliable),
         indeed_unreliable = plogis(b_Intercept + b_DPindeed + b_groupunreliable + `b_DPindeed:groupunreliable`)) |>
  # make into usable data table
  select(actually_reliable, indeed_reliable, actually_unreliable, indeed_unreliable) |>
  pivot_longer(cols = c(actually_reliable, indeed_reliable, actually_unreliable, indeed_unreliable),
               names_to = "condition") |>
  separate_wider_delim(condition, "_", names = c("DP", "group")) |>
  group_by(DP, group)  |> 
  summarise(pred_m = mean(value, na.rm = TRUE),
            pred_0025 = quantile(value, prob = 0.025),
            pred_0975 = quantile(value, prob = 0.975))

# Plot posterior probabilities
group_labels <- c("reliable" = "Reliable", "unreliable" = "Unreliable")

ggplot(posteriors,
       aes(x = DP, 
           y = pred_m,
           fill = DP)) +
  geom_errorbar(aes(ymin = pred_0025,
                    ymax = pred_0975,
                    color = DP),
                width = .1,
                alpha = .7) + 
  geom_line(aes(group = 1)) +
  geom_point(pch = 21,
             size = 3,
             color = "black", key_glyph = draw_key_blank) +
  facet_wrap( ~ group, labeller = as_labeller(group_labels)) + 
  labs(y = "Estimated proportion") +
  theme_minimal() +
  theme(legend.position = c(0.05, 1.15), legend.direction = "vertical", axis.title.x = element_blank(), 
        panel.spacing = unit(2.5, "lines"), panel.grid.major.x = element_blank(),
        legend.box.just = "right", legend.text = ggtext::element_markdown(size = 14),
        plot.margin = margin(13.2, 4.4, 4.4, 4.4), legend.spacing.y = unit(.2, "lines"),
        axis.text.y = element_text(face = "bold", size = 10,
                                   margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.text.x = element_blank(),
        strip.text.x = element_text(face = "bold", size = 20,
                                    margin = margin(t = 0, r = 0, b = 60, l = 0)),
        axis.title.y = element_text(face = "bold", size = 11,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  scale_color_manual(values = c("#E07A5F", "#81B29A")) +
  scale_fill_manual(values = c("#E07A5F", "#81B29A"),
                    labels = c("indeed" = "<b style='color:#81B29A'>Tatsächlich</b>",
                               "actually" = "<b style='color:#E07A5F'>Eigentlich</b>")) +
  guides(color = "none",
         fill = guide_legend(title = NULL, label.position = "right", byrow = TRUE))

# ggsave("estimated-prob-21-01-25.png", height = 4, width = 6, dpi = "retina", bg = "white")

# Fit model
# block effects
MT_xdiff_block <- brm(
  xpos_side_num ~ DP * group * block +
    (1 + DP * block | submission_id) +
    (1 + DP * group * block | item_id),
  data = xDiffExplore,
  family = bernoulli, 
  seed = 1234,
  iter = 2000,
  chains = 4,
  cores = 4,
  # backend = "cmdstanr",
  file  = "MT_xdiff_block.RDS")

# Load model object (fitted model)
# MT_xdiff_block <- readRDS("MT_xdiff_block.RDS")

# Compute empirical proportions for reference
xDiffExplore %>% 
  group_by(group, DP, block, xpos_side) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n)) %>% 
  filter(block %in% c(1, 10))

# Extract estimated probabilities from fitted model
posteriors_block <- MT_xdiff_block |> 
  spread_draws(b_Intercept,
               b_DPindeed,
               b_groupunreliable,
               b_block,
               `b_DPindeed:groupunreliable`,
               `b_DPindeed:block`,
               `b_groupunreliable:block`,
               `b_DPindeed:groupunreliable:block`) |>  
  # map onto each other
  mutate(# at 0 blocks
    actually_reliable = plogis(b_Intercept),
    indeed_reliable = plogis(b_Intercept + b_DPindeed),
    actually_unreliable = plogis(b_Intercept + b_groupunreliable),
    indeed_unreliable = plogis(b_Intercept + b_DPindeed + b_groupunreliable + `b_DPindeed:groupunreliable`),
    # at block 1
    actually_reliable_1 = plogis(b_Intercept + b_block * 1),
    indeed_reliable_1 = plogis(b_Intercept + b_DPindeed + b_block * 1 + `b_DPindeed:block` * 1),
    actually_unreliable_1 = plogis(b_Intercept + b_groupunreliable + b_block * 1 + `b_groupunreliable:block` * 1),
    indeed_unreliable_1 = plogis(b_Intercept + b_DPindeed + b_groupunreliable + `b_DPindeed:groupunreliable` + b_block * 1 + `b_DPindeed:block` * 1  +
                                        `b_groupunreliable:block` * 1 + `b_DPindeed:groupunreliable:block` * 1),
    # at block 10
    actually_reliable_10 = plogis(b_Intercept + b_block * 10),
    indeed_reliable_10 = plogis(b_Intercept + b_DPindeed  + b_block * 10 + `b_DPindeed:block` * 10),
    actually_unreliable_10 = plogis(b_Intercept + b_groupunreliable + b_block * 10 + `b_groupunreliable:block` * 10),
    indeed_unreliable_10 = plogis(b_Intercept + b_DPindeed + b_groupunreliable + `b_DPindeed:groupunreliable` +
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

# Plot posterior probabilities
group_labels <- c("reliable" = "Reliable", "unreliable" = "Unreliable")

x_axis_labels <- c("1.actually" = "1", "10.actually" = "10",
                   "1.indeed" = "1", "10.indeed" = "10")

ggplot(posteriors_block,
       aes(x = interaction(block, DP), 
           y = pred_m,
           fill = DP)) +
  geom_errorbar(aes(ymin = pred_0025,
                    ymax = pred_0975,
                    color = DP),
                width = .1,
                alpha = .7) + 
  geom_line(aes(group = DP), color = "black") +
  geom_point(pch = 21,
             size = 3,
             color = "black", key_glyph = draw_key_blank) +
  facet_wrap( ~ group, labeller = as_labeller(group_labels)) + 
  labs(y = "Estimated proportion") +
  theme_minimal() +
  theme(legend.position = c(0.05, 1.15), legend.direction = "vertical", axis.title.x = element_blank(), 
        panel.spacing = unit(2.5, "lines"), panel.grid.major.x = element_blank(),
        legend.box.just = "right", legend.text = ggtext::element_markdown(size = 14),
        plot.margin = margin(13.2, 4.4, 4.4, 4.4), legend.spacing.y = unit(.2, "lines"),
        axis.text.y = element_text(face = "bold", size = 10,
                                   margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.text.x = element_text(face = "bold", size = 10,
                                   margin = margin(t = 5, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(face = "bold", size = 20,
                                    margin = margin(t = 0, r = 0, b = 60, l = 0)),
        axis.title.y = element_text(face = "bold", size = 11,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  scale_color_manual(values = c("#E07A5F", "#81B29A")) +
  scale_fill_manual(values = c("#E07A5F", "#81B29A"),
                    labels = c("indeed" = "<b style='color:#81B29A'>Tatsächlich</b>",
                               "actually" = "<b style='color:#E07A5F'>Eigentlich</b>")) +
  guides(color = "none",
         fill = guide_legend(title = NULL, label.position = "right", byrow = TRUE)) +
  scale_x_discrete(labels = x_axis_labels)

# ggsave("estimated-prob-block-21-01-25.png", height = 4.5, width = 7, dpi = "retina", bg = "white")
