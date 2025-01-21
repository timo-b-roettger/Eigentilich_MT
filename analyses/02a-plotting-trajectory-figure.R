# Load packages
library(tidyverse)
library(ggtext)
library(ggborderline)

# Import relevant data
results_cleaned <- read_csv("../data/results-clean-06-01-23.csv")

# Extrapolate trajectories from empirical trajectories
time_window_size = 150

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

# Compute mean item trajectories aggregated over participants
plot_tra <- results_cleaned %>%
  group_by(mt_id, DP, group, item_id, condition) %>%
  group_modify(function(d,...) {
    extrapolate_trajectories(d$mt_time, d$mt_x_norm, d$mt_y_norm)}) %>%
  ungroup() %>%
  group_by(DP, group, condition, item_id, t_star) %>%
  summarise(mean_mt_x_star = mean(mt_x_star))
  
# Compute mean onsets of relevant linguistic material in the stimulus
# to be used as marks in the plot
plot_mean_onsets <- results_cleaned %>%
  summarize(mean_DP_onset = mean(DP_onset),
            mean_disambiguation_onset = mean(disambiguation_onset),
            mean_answer_duration = mean(answer_duration))

# Compute group means
mean_tra <- results_cleaned %>%
  group_by(group, mt_id, DP) %>%
  group_modify(function(d,...) {
    extrapolate_trajectories(d$mt_time, d$mt_x_norm, d$mt_y_norm)}) %>%
  ungroup() %>%
  group_by(group, DP, t_star) %>%
  summarise(mean_mt_x_star = mean(mt_x_star))

group_labels <- c("reliable" = "Reliable", "unreliable" = "Unreliable")

# Plot mean item trajectories aggregated over participants
# with group meas overlaid
plot_tra %>%
  ggplot() +
  # rectangle shading time zone
  geom_rect(
    data = tibble(
      `Time window` = c("Das ist", "PARTICLE", "NOUN") %>%
        factor(ordered = T, levels = c("Das ist", "PARTICLE",
                                       "NOUN")),
      start = c(0,
                plot_mean_onsets$mean_DP_onset[1],
                plot_mean_onsets$mean_disambiguation_onset[1]),
      end = c(plot_mean_onsets$mean_DP_onset[1],
              plot_mean_onsets$mean_disambiguation_onset[1],
              plot_mean_onsets$mean_answer_duration[1])
    ),
    aes(NULL, NULL, xmin = start, xmax = end, fill = `Time window`),
    ymin = -.5,
    ymax = 1,
    colour = "black",
    # size = .7,
    alpha = .7
  ) +
  scale_fill_brewer(
    palette = "Greys"
  ) +
  # horizontal lines indicating time zones for means
  geom_vline(
    data = tibble(
      xintercept = with(plot_mean_onsets, c(
        mean(mean_DP_onset) + time_window_size,
        mean(mean_DP_onset) - time_window_size,
        mean(mean_disambiguation_onset) + time_window_size,
        mean(mean_disambiguation_onset) - time_window_size
      ))
    ),
    aes(xintercept = xintercept), alpha = 0.7, linetype = "longdash"
  ) +
  # actual mouse path
  geom_path(aes(x = t_star, y = mean_mt_x_star, color = DP, group = interaction(item_id, DP)),
            linewidth = .5, lineend = "round", linejoin = "bevel", alpha = .3, key_glyph = draw_key_blank) +
  ggborderline::geom_borderline(data = mean_tra,
                                aes(x = t_star, y = mean_mt_x_star, color = DP),
                                linewidth = 1, borderwidth = .3, lineend = "round", linejoin = "bevel", key_glyph = draw_key_blank) +
  theme_minimal() +
  facet_wrap( ~ group, labeller = as_labeller(group_labels)) +
  theme(legend.position = c(.38, 1.15), legend.direction = "vertical",
        axis.title.x = element_blank(), legend.text = ggtext::element_markdown(size = 14),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
        axis.title.y = element_blank(), legend.box.just = "right", panel.spacing = unit(2.5, "lines"),
        axis.text.y = element_text(face = "bold", size = 8, 
                                   margin = margin(t = 0, r = 5, b = 0, l = 0)), 
        axis.text.x = element_text(face = "bold", size = 8, 
                                   margin = margin(t = 5, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(face = "bold", size = 20,
                                    margin = margin(t = 0, r = 0, b = 60, l = 0)),
        plot.margin = margin(13.2, 4.4, 8.8, 4.4), legend.spacing.y = unit(.2, "lines")) +
  scale_y_continuous(breaks = seq(-1, 1, .25)) +
  scale_x_continuous(breaks = seq(0, 4000, 500)) +
  scale_color_manual(values = c("#E07A5F", "#81B29A"),
                     labels = c("indeed" = "<b style='color:#81B29A'>Tatsächlich</b>",
                                "actually" = "<b style='color:#E07A5F'>Eigentlich</b>")) +
  guides(fill = "none",
         color = guide_legend(title = NULL, label.position = "right", byrow = TRUE)) +
  coord_cartesian(ylim = c(-.5, 1), clip = "off") +
  geom_segment(data = data.frame(group = "reliable"),
           aes(x = plot_mean_onsets$mean_DP_onset[1],
           xend = plot_mean_onsets$mean_DP_onset[1],
           y = 1.35,
           yend = 1.05), color = '#bdbdbd', lwd = 2,
           arrow = arrow(length = unit(0.2, "inches"), ends = "last")) +
  geom_label(data = data.frame(group = "reliable"), size = 10/ .pt,
           aes(x = plot_mean_onsets$mean_DP_onset[1],
           y = 1.45), color = "white", fill = "#bdbdbd", label = "PARTICLE") +
  geom_segment(data = data.frame(group = "reliable"),
           aes(x = plot_mean_onsets$mean_disambiguation_onset[1],
           xend = plot_mean_onsets$mean_disambiguation_onset[1],
           y = 1.35,
           yend = 1.05), color = '#636363', lwd = 2.3,
           arrow = arrow(length = unit(0.2, "inches"), ends = "last")) +
  geom_label(data = data.frame(group = "reliable"), size = 10/ .pt,
           aes(x = plot_mean_onsets$mean_disambiguation_onset[1],
           y = 1.45), color = "white", fill = "#636363", label = "NOUN")

# ggsave("MT-main-21-01-25.png", height = 4, width = 8, dpi = "retina", bg = "white")
