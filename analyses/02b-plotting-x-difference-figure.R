# Load packages
library(tidyverse)

# Import relevant data
results_cleaned <- read_csv("../data/results-clean-06-01-23.csv")
xDiffExplore <- read_csv("../data/x-difference-03-12-24.csv")

# Extrapolate trajectories from empirical trajectories
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


# Compute mean trajectories aggregated over items and participants
# in the reliable group
reliable_tra <- results_cleaned %>%
  filter(group == "reliable") %>%
  group_by(mt_id, DP) %>%
  group_modify(function(d,...) {
    extrapolate_trajectories(d$mt_time, d$mt_x_norm, d$mt_y_norm)}) %>%
  ungroup() %>%
  group_by(DP, t_star) %>%
  summarise(mean_mt_x_star = mean(mt_x_star))

# Compute mean trajectories aggregated over items and participants
# in the unreliable group
unreliable_tra <- results_cleaned %>%
  filter(group == "unreliable") %>%
  group_by(mt_id, DP) %>%
  group_modify(function(d,...) {
    extrapolate_trajectories(d$mt_time, d$mt_x_norm, d$mt_y_norm)}) %>%
  ungroup() %>%
  group_by(DP, t_star) %>%
  summarise(mean_mt_x_star = mean(mt_x_star))

# Plot differences in horizontal mouse position
# for each participant-item pair in the reliable group
xdiff_rel <- xDiffExplore %>%
  filter(group == "reliable") %>%
  group_by(DP) %>%
  mutate(mean = mean(xpos_mean)) %>%
  ungroup() %>% 
  ggplot() +
  geom_hline(aes(yintercept = mean, color = DP),
             linetype = "dashed", linewidth = .8) +
  geom_histogram(aes(y = xpos_mean, fill = DP), color = "black", position = "dodge") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA),
        axis.title.y = element_blank(), strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(), panel.spacing = unit(.8, "lines"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12, 
                                   margin = margin(t = 0, r = 5, b = 0, l = 0)),
        plot.margin = margin(8.8, 4.4, 4.4, 8.8)) +
  scale_fill_manual(values = c("#E07A5F", "#81B29A")) +
  scale_color_manual(values = c("#E07A5F", "#81B29A"))

# Plot differences in horizontal mouse position
# for each participant-item pair in the unreliable group
xdiff_unrel <- xDiffExplore %>%
  filter(group == "unreliable") %>%
  group_by(DP) %>%
  mutate(mean = mean(xpos_mean)) %>%
  ungroup() %>%
  ggplot() +
  geom_hline(aes(yintercept = mean, color = DP),
             linetype = "dashed", linewidth = .8) +
  geom_histogram(aes(y = xpos_mean, fill = DP), color = "black", position = "dodge") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_rect(fill = NA),
        axis.title.y = element_blank(), strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(), panel.spacing = unit(.8, "lines"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12, 
                                   margin = margin(t = 0, r = 5, b = 0, l = 0)),
        plot.margin = margin(8.8, 4.4, 4.4, 8.8)) +
  scale_fill_manual(values = c("#E07A5F", "#81B29A")) +
  scale_color_manual(values = c("#E07A5F", "#81B29A"))

# Compute mean onsets of relevant linguistic material in the stimulus
# to be used as marks in the plot
plot_mean_onsets <- results_cleaned %>%
  summarize(mean_DP_onset = mean(DP_onset),
            mean_disambiguation_onset = mean(disambiguation_onset),
            mean_answer_duration = mean(answer_duration))

# Plot highlighted mean trajectories aggregated over items and participants
# in the reliable group
mt_rel_high <- reliable_tra %>%
  mutate(highlight = case_when((plot_mean_onsets$mean_DP_onset[1] + 150 < t_star & plot_mean_onsets$mean_disambiguation_onset[1] + 150 > t_star | t_star == plot_mean_onsets$mean_disambiguation_onset[1] | t_star == plot_mean_onsets$mean_DP_onset[1]) & DP == "indeed" ~ "on_indeed",
                               (plot_mean_onsets$mean_DP_onset[1] + 150 < t_star & plot_mean_onsets$mean_disambiguation_onset[1] + 150 > t_star) & DP == "actually" ~ "on_actually",
                               (plot_mean_onsets$mean_DP_onset[1] + 150 > t_star | plot_mean_onsets$mean_disambiguation_onset[1] + 150 < t_star) & DP == "indeed" ~ "off_indeed",
                               (plot_mean_onsets$mean_DP_onset[1] + 150 > t_star | plot_mean_onsets$mean_disambiguation_onset[1] + 150 < t_star) & DP == "actually" ~ "off_actually")) %>% 
  select(t_star, mean_mt_x_star, highlight) %>%
  ggplot() +
  geom_rect(
    data = tibble(
      `Time window` = c("Particle onset", "Disambiguation") %>%
        factor(ordered = T, levels = c("Particle onset",
                                       "Disambiguation")),
      start = c(plot_mean_onsets$mean_DP_onset[1] + 150,
                plot_mean_onsets$mean_disambiguation_onset[1]),
      end = c(plot_mean_onsets$mean_disambiguation_onset[1],
              plot_mean_onsets$mean_disambiguation_onset[1] + 150)),
    aes(NULL, NULL, xmin = start, xmax = end, fill = `Time window`),
    ymin = -.25,
    ymax = 1,
    colour = "black",
    linewidth = .3,
    alpha = .7
  ) +
  scale_fill_manual(
    values = c("#bdbdbd", "#636363")
  ) +
  geom_path(aes(x = t_star, y = mean_mt_x_star, color = highlight, group = DP),
            linewidth = 1.2, lineend = "round", linejoin = "bevel", key_glyph = draw_key_blank) +
  theme_minimal() +
  ggtitle("Reliable") +
  theme(legend.position = c(.8, 1.4), legend.direction = "vertical",
        panel.border = element_rect(fill = NA), legend.text = ggtext::element_markdown(size = 16),
        axis.title.x = element_blank(), strip.text.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
        axis.title.y = element_blank(), panel.spacing = unit(.8, "lines"),
        axis.text.y = element_text(face = "bold", size = 10, 
                                   margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(face = "bold", size = 8, 
                                   margin = margin(t = 5, r = 0, b = 5, l = 0)),
        plot.title = element_text(face = "bold", size = 24, hjust = .5,
                                   margin = margin(t = 0, r = 0, b = 70, l = 0)),
        plot.margin = margin(4.4, 8.8, 4.4, 4.4), legend.spacing.y = unit(.2, "lines")) +
  scale_y_continuous(breaks = seq(-1, 1, .25)) +
  scale_x_continuous(breaks = c(plot_mean_onsets$mean_DP_onset[1] + 150, plot_mean_onsets$mean_disambiguation_onset[1] + 150),
                     labels = c(signif(plot_mean_onsets$mean_DP_onset[1] + 150, 2), signif(plot_mean_onsets$mean_disambiguation_onset[1] + 150, 2))) +
  scale_color_manual(values = c("#f5d7cf", "#d9e7e0", "#E07A5F", "#81B29A"),
                     labels = c("on_indeed" = "<b style='color:#81B29A'>Tatsächlich</b>",
                                "off_indeed" = "",
                                "on_actually" = "<b style='color:#E07A5F'>Eigentlich</b>",
                                "off_actually" = "")) +
  guides(fill = "none",
         color = guide_legend(title = NULL, label.position = "right", byrow = TRUE)) +
  ggpubr::geom_bracket(xmin = plot_mean_onsets$mean_DP_onset[1],
                       xmax = plot_mean_onsets$mean_disambiguation_onset[1] + 300,
                       y.position = 1.15, size = 1.8, color = "grey",
                       tip.length = c(.1, .1), label = "") +
  coord_cartesian(ylim = c(-.25, 1), clip = "off") +
  geom_label(aes(x = (plot_mean_onsets$mean_disambiguation_onset[1]) - (plot_mean_onsets$mean_DP_onset[1] + 150), y = 1.4),
           label.padding = unit(.4, "lines"), size = 14/ .pt, color = "white", fill = "black", label = "Critical window")
  

# Plot highlighted mean trajectories aggregated over items and participants
# in the unreliable group
mt_unrel_high <- unreliable_tra %>%
  mutate(highlight = case_when((plot_mean_onsets$mean_DP_onset[1] + 150 < t_star & plot_mean_onsets$mean_disambiguation_onset[1] + 150 > t_star | t_star == plot_mean_onsets$mean_disambiguation_onset[1] | t_star == plot_mean_onsets$mean_DP_onset[1]) & DP == "indeed" ~ "on_indeed",
                               (plot_mean_onsets$mean_DP_onset[1] + 150 < t_star & plot_mean_onsets$mean_disambiguation_onset[1] + 150 > t_star) & DP == "actually" ~ "on_actually",
                               (plot_mean_onsets$mean_DP_onset[1] + 150 > t_star | plot_mean_onsets$mean_disambiguation_onset[1] + 150 < t_star) & DP == "indeed" ~ "off_indeed",
                               (plot_mean_onsets$mean_DP_onset[1] + 150 > t_star | plot_mean_onsets$mean_disambiguation_onset[1] + 150 < t_star) & DP == "actually" ~ "off_actually")) %>% 
  select(t_star, mean_mt_x_star, highlight) %>%
  ggplot() +
  geom_rect(
    data = tibble(
      `Time window` = c("Particle onset", "Disambiguation") %>%
        factor(ordered = T, levels = c("Particle onset",
                                       "Disambiguation")),
      start = c(plot_mean_onsets$mean_DP_onset[1] + 150,
                plot_mean_onsets$mean_disambiguation_onset[1]),
      end = c(plot_mean_onsets$mean_disambiguation_onset[1],
              plot_mean_onsets$mean_disambiguation_onset[1] + 150)),
    aes(NULL, NULL, xmin = start, xmax = end, fill = `Time window`),
    ymin = -.25,
    ymax = 1,
    colour = "black",
    linewidth = .3,
    alpha = .7
  ) +
  scale_fill_manual(
    values = c("#bdbdbd", "#636363")
  ) +
  geom_path(aes(x = t_star, y = mean_mt_x_star, color = highlight, group = DP),
          linewidth = 1.2, lineend = "round", linejoin = "bevel", key_glyph = draw_key_blank) +
  theme_minimal() +
  ggtitle("Unreliable") +
  theme(legend.position = "none", legend.direction = "vertical",
        panel.border = element_rect(fill = NA), legend.text = ggtext::element_markdown(size = 18),
        axis.title.x = element_blank(), strip.text.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
        axis.title.y = element_blank(), panel.spacing = unit(.8, "lines"),
        axis.text.y = element_text(face = "bold", size = 10, 
                                   margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(face = "bold", size = 8, 
                                   margin = margin(t = 5, r = 0, b = 5, l = 0)),
        plot.title = element_text(face = "bold", size = 24, hjust = .5,
                                  margin = margin(t = 0, r = 0, b = 70, l = 0)),
        plot.margin = margin(4.4, 8.8, 4.4, 4.4), legend.spacing.y = unit(.5, "lines")) +
  scale_y_continuous(breaks = seq(-1, 1, .25)) +
  scale_x_continuous(breaks = c(plot_mean_onsets$mean_DP_onset[1] + 150, plot_mean_onsets$mean_disambiguation_onset[1] + 150),
                     labels = c(signif(plot_mean_onsets$mean_DP_onset[1] + 150, 2), signif(plot_mean_onsets$mean_disambiguation_onset[1] + 150, 2))) +
  scale_color_manual(values = c("#f5d7cf", "#d9e7e0", "#E07A5F", "#81B29A"),
                     labels = c("on_indeed" = "<b style='color:#81B29A'>Tatsächlich</b>",
                                "off_indeed" = "",
                                "on_actually" = "<b style='color:#E07A5F'>Eigentlich</b>",
                                "off_actually" = "")) +
  guides(fill = "none",
         color = guide_legend(title = NULL, label.position = "right", byrow = TRUE)) +
  coord_cartesian(ylim = c(-.25, 1), clip = "off")
  
# Combine highlighted mean trajectories plots
mt_high <- ggpubr::ggarrange(mt_rel_high, mt_unrel_high)

# Combine distributions plots
xdiff <- ggpubr::ggarrange(xdiff_rel, xdiff_unrel)

# Combine highlighted mean trajectories plot (top) with mouse position distributions plot (bottom)
ggpubr::ggarrange(mt_high, xdiff, nrow = 2, heights = c(.8, .6))

# ggsave("xdiff-full-21-01-25.png", height = 6, width = 7, dpi = "retina", bg = "white")

# Plot differences in horizontal mouse position
# for each participant-item pair in the reliable group
# across blocks
xdiff_rel_block <- xDiffExplore %>%
  filter(group == "reliable") %>%
  group_by(DP, block) %>%
  mutate(mean = mean(xpos_mean)) %>%
  ungroup() %>% 
  ggplot() +
  ggtitle("Reliable") +
  geom_hline(aes(yintercept = mean, color = DP),
             linetype = "dashed", linewidth = .8) +
  geom_histogram(aes(y = xpos_mean, fill = DP), color = "black", position = "dodge") +
  facet_grid( ~ block) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 24, hjust = .5,
                                  margin = margin(t = 0, r = 0, b = 20, l = 0)),
        panel.border = element_rect(fill = NA),
        axis.title.y = element_blank(), strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(), panel.spacing = unit(.8, "lines"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12, 
                                   margin = margin(t = 0, r = 5, b = 0, l = 0)),
        plot.margin = margin(8.8, 4.4, 4.4, 8.8)) +
  scale_fill_manual(values = c("#E07A5F", "#81B29A")) +
  scale_color_manual(values = c("#E07A5F", "#81B29A"))

# ggsave("xdiff-rel-block-13-12-24.png", height = 6, width = 11, dpi = "retina", bg = "white")

# Plot differences in horizontal mouse position
# for each participant-item pair in the unreliable group
# across blocks
xdiff_unrel_block <- xDiffExplore %>%
  filter(group == "unreliable") %>%
  group_by(DP, block) %>%
  mutate(mean = mean(xpos_mean)) %>%
  ungroup() %>%
  ggplot() +
  ggtitle("Unreliable") +
  geom_hline(aes(yintercept = mean, color = DP),
             linetype = "dashed", linewidth = .8) +
  geom_histogram(aes(y = xpos_mean, fill = DP), color = "black", position = "dodge") +
  facet_grid( ~ block) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 24, hjust = .5,
                                  margin = margin(t = 0, r = 0, b = 20, l = 0)),
        panel.border = element_rect(fill = NA),
        axis.title.y = element_blank(), strip.text.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(), panel.spacing = unit(.8, "lines"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 12, 
                                   margin = margin(t = 0, r = 5, b = 0, l = 0)),
        plot.margin = margin(8.8, 4.4, 4.4, 8.8)) +
  scale_fill_manual(values = c("#E07A5F", "#81B29A")) +
  scale_color_manual(values = c("#E07A5F", "#81B29A"))

# ggsave("xdiff-unrel-block-13-12-24.png", height = 6, width = 11, dpi = "retina", bg = "white")
