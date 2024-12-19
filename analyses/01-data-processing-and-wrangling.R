# Load packages
library(mousetrap)
library(tidyverse)
library(rstudioapi)

# Import raw data outputted from experiment ------------------------------------

# set the current working directory to the one where this file is
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

raw_data <- read_csv("../data/results_75_disc-part-proc_VMS.csv") %>% 
  filter(trialType == 'main', condition != "filler") %>% 
  #create manual id for each trial 
  mutate(
    mt_id = str_c(submission_id, "_", trialNumber)
  )

# Add acoustic landmarks -------------------------------------------------------

raw_data <- raw_data %>% 
  mutate(
    audio_file_name_minimal = str_sub(answer_file, start = 3, end = -5)
  )

# Load information about temporal landmarks in all answer file recordings

audio_landmarks <- read_csv("../data/acousticLandmarks.csv") %>% 
  # only critical trials matter 
  filter(! str_detect(stimulus, "f_")) %>% 
  # mutate/rename the relevant information
  #    cast time info to miliseconds
  mutate(
    audio_file_name_minimal = stimulus,
    DP_onset = adverb_onset_rel * 1000,
    disambiguation_onset = referent_onset_rel * 1000,
    answer_duration = duration * 1000
  ) %>% 
  # select only the relevant (onset relative) information
  select(
    audio_file_name_minimal,
    DP_onset,
    disambiguation_onset,
    answer_duration
  ) %>% 
  filter(audio_file_name_minimal %in% raw_data$audio_file_name_minimal)

raw_data <- full_join(raw_data, audio_landmarks, by = "audio_file_name_minimal")

# Compute non-normalized trajectories and absolute timestamps from scratch -----

extract_mt_triples = function(mt_time,mt_x,mt_y) {
  unroll <- function(x) {as.numeric(unlist(stringr::str_split(x, "\\|")))}
  tibble(
    mt_time = unroll(mt_time),
    mt_x = unroll(mt_x),
    mt_y = unroll(mt_y),
    c = 1:length(mt_time)
  ) %>% 
    arrange(-c) %>% 
    distinct(mt_time, .keep_all = TRUE) %>% 
    arrange(c) %>% 
    select(-c) %>% 
    mutate(
      mt_nobs = length(mt_x),
      mt_step = 1:mt_nobs[1],
      mt_time_start = as.numeric(mt_time[1]),
      mt_time_end  = as.numeric(mt_time[mt_nobs]),
      mt_time_duration = mt_time_end - mt_time_start,
      mt_x_end    = as.numeric(mt_x[mt_nobs]),
      mt_x_start  = as.numeric(mt_x[1]),
      mt_y_end    = as.numeric(mt_y[mt_nobs]),
      mt_y_start  = as.numeric(mt_y[1]),
      # symmetrize to the left
      mt_x = ifelse(mt_x_end > mt_x_start, mt_x, -mt_x),
      mt_x_end = as.numeric(mt_x[mt_nobs]),
      mt_x_start = as.numeric(mt_x[1]),
      # # all trajectories go up
      mt_y = ifelse(mt_y_end > mt_y_start, mt_y, -mt_y),
      mt_y_end    = as.numeric(mt_y[mt_nobs]),
      mt_y_start  = as.numeric(mt_y[1]),
      # 'space-normalization' -> all traject. start at (0,0) & end at (1,1)
      mt_x_norm = (mt_x - mt_x_start)/abs(mt_x_end - mt_x_start),
      mt_y_norm = (mt_y - mt_y_start)/abs(mt_y_end - mt_y_start)
    )
}

raw_mt_data <- raw_data %>% 
  filter(! is.na(mt_time)) %>% 
  mutate(
    mt_nobs = stringr::str_count(mt_x, "\\|")
    # mt_time_start = 
  ) %>% 
  select(
    mt_id,
    response_time,
    mt_nobs,
    mt_time,
    mt_x,
    mt_y,
    mt_start_time, 
    experiment_start_time, 
    experiment_end_time, 
    experiment_duration, 
    mt_time
  ) %>% 
  unique()

raw_mt_data_long <- raw_mt_data %>% 
  group_by(mt_id) %>% 
  nest(mt_data_wide = c(2:10)) %>%
  summarise(
    mt_data_long <- map(mt_data_wide, function(d) {
      extract_mt_triples(d$mt_time, d$mt_x, d$mt_y)
    })
  ) %>%
  unnest(c(2)) %>% 
  ungroup()

# Importing & processing with mousetrap ----------------------------------------

# Create a mousetrap object
mtdata <- mt_import_mousetrap(
  raw_data,
  xpos_label = "mt_x", # see x,y,time column names in the _magpie output
  ypos_label = "mt_y",
  timestamps_label = "mt_time",
  mt_id_label = "mt_id",
  split = "\\|", # see sample-separating character used in _magpie output
  verbose = TRUE
)

# Create "nobs" (per-trial sample count ("Number of OBServations"))
mtdata <- mt_count(mtdata, save_as = "data")

# Flip all trajectories to the left & up
mtdata <- mt_remap_symmetric(mtdata, remap_xpos = "left", remap_ypos = "up")
# Transfer all trajectories to a unit space & align their start & end points
mtdata <- mt_align_start_end(mtdata, start = c(0,0), end = c(1,1))

# Calculate distance, velocity, and acceleration
mtdata <- mt_derivatives(mtdata)
# Calculate mouse-tracking measures (MAD, AUC etc.)
mtdata <- mt_measures(mtdata)

# Add time-normalized version of the trajectories to the array
# These are most often used for:
# a) plotting averaged trajectories
# b) the comparison of x,y positions at different relative time points
mtdata <- mt_time_normalize(mtdata, save_as = "tn_trajectories")

# Joining/exporting everything to a data frame

mtdata_measures <- mtdata$data %>%
  full_join(mtdata$measures, by = "mt_id")

mtdata_raw_trajectories <- as.data.frame.table(mtdata$trajectories) %>%
  spread(key = Var3, value = Freq) %>%
  select(-Var2) %>%
  rename("mt_id" = Var1)

results_mt <- 
  full_join(mtdata_measures, mtdata_raw_trajectories, by = "mt_id") %>% 
  filter(!is.na(xpos))

mtdata_tn_trajectories <- as.data.frame.table(mtdata$tn_trajectories) %>%
  spread(key = Var3, value = Freq) %>%
  select(-Var2) %>%
  rename("mt_id" = Var1)

results_tn <-
  full_join(mtdata_measures, mtdata_tn_trajectories, by = "mt_id")

# Cross-reference data frames --------------------------------------------------

results <- cbind(
  results_mt %>% arrange(mt_id, timestamps),
  raw_mt_data_long %>% arrange(mt_id, mt_time) %>% select(-mt_id)
)

results <- results %>%
  mutate(response_tc = case_when(
    condition == "filler" ~ ifelse(MT_target == response, "target", "competitor"),
    condition == "reliable" & DP == "actually" ~ ifelse(MT_target != response, "target", "competitor"),
    condition == "reliable" & DP == "indeed" ~ ifelse(MT_target == response, "target", "competitor"),
    condition == "unreliable" & DP == "actually" ~ ifelse(MT_target == response, "target", "competitor"),
    condition == "unreliable" & DP == "indeed" ~ ifelse(MT_target != response, "target", "competitor"))
  )


# Extrapolate trajectories from empirical trajectories -------------------------

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

resultsExtrapolatedTime <- results %>% group_by(mt_id) %>% 
  group_modify(function(d,...) {extrapolate_trajectories(d$mt_time, d$mt_x_norm, d$mt_y_norm)}) %>%  
  ungroup()

# Cleaning the data -----------------------------------------------------------
message('exclusions')

message(
  "Number of trials with total RT longer than 1sec + answer duration: ", 
  results %>% filter(RT > answer_duration + 1000) %>% select(mt_id) %>% unique() %>% nrow())

round((results %>% filter(RT > answer_duration + 1000) %>% select(mt_id) %>% unique() %>% nrow())/ (results %>% select(mt_id) %>% unique() %>% nrow()), digits = 2)

message(
  'Number of trials with at most 10 mouse-coordinate observations: ',
  results %>% filter(mt_nobs <= 10) %>% select(mt_id) %>% unique() %>% nrow()
)

round((results %>% filter(mt_nobs <= 10) %>% select(mt_id) %>% unique() %>% nrow())/ (results %>% select(mt_id) %>% unique() %>% nrow()), digits = 2)

message(
  'Number of trials with exceptionally large x-pos excess: ',
  results %>% group_by(mt_id) %>% 
    mutate(outOfXRange = !(max(mt_x_norm) < 1.2 & min(mt_x_norm) > -1.2)) %>% 
    ungroup() %>% 
    filter(outOfXRange == TRUE) %>% select(mt_id) %>% unique() %>% nrow()
)

message(
  'Number of trials that did not end selecting the correct picture: ',
  results %>% 
    filter(response_tc != "target") %>% select(mt_id) %>% unique() %>% nrow()
)

round((results %>% filter(response_tc != "target") %>% select(mt_id) %>% unique() %>% nrow())/ (results %>% select(mt_id) %>% unique() %>% nrow()), digits = 2)

message(
  'Number of trials in the unreliable condition: ',
  results %>% 
    filter(condition == "unreliable") %>% select(mt_id) %>% unique() %>% nrow()
)

message(
  'Number of early onset trials: ',
  results %>%
    group_by(mt_id) %>%
    mutate(early_onset = mt_time_end < DP_onset) %>%
    ungroup() %>% 
    filter(early_onset == TRUE) %>% select(mt_id) %>% unique() %>% nrow()
)

message(
  'Number of early onsetters: ',
  results %>%
    group_by(mt_id) %>%
    mutate(early_onset = mt_time_end < DP_onset) %>% 
    group_by(submission_id) %>% 
    mutate(early_onsetter = mean(early_onset, na.rm = T) > 0.3) %>%
    filter(early_onsetter == TRUE) %>% select(submission_id) %>% unique() %>% nrow()
)


results_cleaned <- results %>%
  group_by(mt_id) %>%
  mutate(early_onset = mt_time_end < DP_onset) %>%
  group_by(submission_id) %>% 
  mutate(early_onsetter = mean(early_onset, na.rm = T) > 0.3) %>%
  filter(early_onsetter == FALSE) %>%
  # filter(early_onset == FALSE) %>% 
  filter(RT < answer_duration + 1000) %>% 
  filter(mt_nobs > 10) %>% 
  filter(response_tc == "target") %>% 
  filter(condition == "reliable") %>% 
  group_by(mt_id) %>% 
  mutate(outOfXRange = !(max(mt_x_norm) < 1.2 & min(mt_x_norm) > -1.2)) %>% 
  ungroup() %>% 
  filter(outOfXRange == FALSE)

message(
  'Total number of trials excluded with these criteria: ',
  setdiff(results$mt_id, results_cleaned$mt_id) %>% unique %>% length()
)

round((setdiff(results$mt_id, results_cleaned$mt_id) %>% unique %>% length())/ (results$mt_id %>% unique %>% length()), digits = 2)

# applying the same exclusion criteria to the extrapolated data
resultsExtrapolatedTime_cleaned <- resultsExtrapolatedTime %>% 
  filter(mt_id %in% (results_cleaned$mt_id %>% unique()))

# produce data set with cleaned data for plotting
# write.csv(results_cleaned, "../data/results-clean-06-01-23.csv", row.names = FALSE, fileEncoding = 'UTF-8')


# Compute bias in horizontal mouse position within the particle window  ------

extrapolate <- function(x,t,t_star){
  upper_index <- which(t>t_star)[1]
  lower_index <- max(which(t<t_star))
  
  x_star <- x[lower_index] +
    (x[upper_index] - x[lower_index]) *
    (t_star - t[lower_index]) /
    (t[upper_index] - t[lower_index])
  return(x_star)
}

# Plot distribution of measured horizontal mouse positions per trial

group_labels <- c("reliable" = "Reliable", "unreliable" = "Unreliable")

resultsExtrapolatedTime_cleaned |> 
  full_join(
    results_cleaned |> 
      select(mt_id, submission_id, group, DP, item_id, DP_onset, disambiguation_onset, block, trialNumber) %>% unique(),
    by = "mt_id"
  ) |> 
  mutate(xpos_critical = ifelse(t_star > (DP_onset + 150) & (disambiguation_onset + 150) > t_star,
                                1, 0)) |> 
  filter(xpos_critical == 1) |> 
  group_by(DP, group, submission_id, item_id, block, trialNumber) |>
  summarize(xpos_points = n()) |>
  ggplot() +
  facet_grid(DP ~ group, labeller = as_labeller(group_labels)) +
  geom_histogram(aes(xpos_points, fill = DP), color = "black", key_glyph = draw_key_blank) +
  theme_minimal() +
  theme(legend.position = c(1.15, .5), legend.text = ggtext::element_markdown(size = 12),
        axis.title.y = element_blank(), strip.text.y = element_blank(), legend.title = element_blank(),
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
    scale_fill_manual(values = c("#E07A5F", "#81B29A"),
                      labels = c("indeed" = "<b style='color:#81B29A'>Tatsächlich</b>",
                                 "actually" = "<b style='color:#E07A5F'>Eigentlich</b>"))
  
# ggsave("n-points-xpos-16-12-24.png", height = 4, width = 5, dpi = "retina", bg = "white")
 
# Compute mean horizontal mouse position within critical window per trial  
xDiffExplore <-
resultsExtrapolatedTime_cleaned |>
  full_join(
    results_cleaned |>
      select(mt_id, submission_id, group, DP, item_id, DP_onset, disambiguation_onset, block, trialNumber) %>% unique(),
    by = "mt_id"
  ) |>
  mutate(xpos_critical = ifelse(t_star > (DP_onset + 150) & (disambiguation_onset + 150) > t_star,
                                1, 0)) |> 
  filter(xpos_critical == 1) |> 
  group_by(DP, group, submission_id, item_id, block, trialNumber) |> 
  summarize(xpos_mean = mean(mt_x_star, na.rm = TRUE)) |>
  mutate(xpos_side = ifelse(xpos_mean > 0, "target", "competitor"),
         xpos_side_num = ifelse(xpos_mean > 0, 1, 0))


# Produce data set with mean horizontal positions for plotting
# write.csv(xDiffExplore, "../data/x-difference-03-12-24.csv", row.names = FALSE, fileEncoding = 'UTF-8')
