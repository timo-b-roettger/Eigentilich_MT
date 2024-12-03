# Load packages
library(tidyverse)

# Import relevant data
results_cleaned <- read_csv("../data/results-clean-06-01-23.csv")
clusters_extrapolated <- read_csv("../data/clusters-06-01-23.csv")

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

# Compute number of trajectories per cluster
counts <- results_cleaned %>%
  full_join(clusters_extrapolated, by = "mt_id") %>%
  group_by(mt_id, submission_id, DP, group, item_id, cluster, condition) %>%
  group_modify(function(d,...) {
    extrapolate_trajectories(d$mt_time, d$mt_x_norm, d$mt_y_norm)}) %>%
  ungroup() %>% 
  group_by(DP, group, cluster) %>%
  summarize(n_cluster = n()) %>%
  ungroup() %>%
  group_by(DP, group) %>%
  mutate(perc_clust = round(n_cluster/ sum(n_cluster), digits = 2))

group_labels <- c("reliable" = "Reliable", "unreliable" = "Unreliable")

# Plot proportion of trajectories per cluster
results_cleaned %>%
  full_join(clusters_extrapolated, by = "mt_id") %>%
  group_by(mt_id, DP, group, cluster, condition) %>%
  group_modify(function(d,...) {
    extrapolate_trajectories(d$mt_time, d$mt_x_norm, d$mt_y_norm)}) %>%
  ungroup() %>%
  group_by(cluster, DP, group) %>% 
  summarize(n_cluster = n()) %>%
  ungroup() %>% 
  group_by(DP, group) %>%
  mutate(perc_clust = n_cluster/ sum(n_cluster)) %>%
  ungroup() %>%
  ggplot() +
  geom_col(aes(x = cluster, y = perc_clust, fill = DP), color = "black", key_glyph = draw_key_blank) +
  facet_grid(DP ~ group, labeller = as_labeller(group_labels)) +
  theme_minimal() +
  theme(legend.position = c(1.15, .5), legend.text = ggtext::element_markdown(size = 12),
        axis.title.y = element_blank(), strip.text.y = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(), panel.spacing.y = unit(1.2, "lines"),
        axis.title.x = element_blank(), panel.spacing.x = unit(2.5, "lines"),
        axis.text.x = element_text(face = "bold", size = 10,
                                   margin = margin(t = 5, r = 0, b = 0, l = 0)),
        strip.text.x = element_text(face = "bold", size = 14,
                                    margin = margin(t = 0, r = 0, b = 20, l = 0)),
        axis.text.y = element_text(face = "bold", size = 8, 
                                   margin = margin(t = 0, r = 5, b = 0, l = 0)),
        plot.margin = margin(4.4, 88, 4.4, 4.4), legend.spacing.y = unit(.01, "lines")) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_fill_manual(values = c("#E07A5F", "#81B29A"),
                     labels = c("indeed" = "<b style='color:#81B29A'>Tatsächlich</b>",
                                "actually" = "<b style='color:#E07A5F'>Eigentlich</b>")) +
  scale_x_discrete(breaks = seq(1, 4, 1),
                   labels = c("1" = "Left",
                              "2" = "Middle",
                              "3" = "Right",
                              "4" = " ")) +
  guides(color = "none",
         fill = guide_legend(title = NULL, label.position = "right", byrow = TRUE)) +
  coord_cartesian(xlim = c(1, 3), clip = "off")

# ggsave("clusters-28-04-22.png", height = 4, width = 5, dpi = "retina", bg = "white")