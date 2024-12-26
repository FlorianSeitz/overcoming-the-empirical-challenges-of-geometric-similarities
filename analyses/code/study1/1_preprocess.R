# ==========================================================================
# Analyses Study 1: Preprocesses data
# ==========================================================================

rm(list = ls(all = TRUE))

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # changes working directory to folder where this script is

# ==========================================================================
# Reads data and performs preprocessing
# ==========================================================================
files <- list.files(path = "../../../data/raw/main_study1_experiment/", pattern = "", full.names = TRUE)
dt <- rbindlist(lapply(files, fread, fill = TRUE, colClasses = list("character" = c("stim", "color_presentation_order", "feature_presentation_order"))))
demographics <- fread("../../../data/raw/main_study1_demographics/main_study1_demographics.csv")

# Makes stimuli start with 0
dt[, stimulus := str_sub(paste0("00", as.numeric(stim) - 111), -3, -1)]
dt[, paste0("f", 1:3) := tstrsplit(stimulus, "")]

# Checks if time limit is exceeded in > 50% of test trials for a given test stimulus
valid_subj_ids <- dt[time_pressure_cond == TRUE & block == "test", .(perc_too_slow = mean(too_slow)), by = list(subj_id, stim)][, all(perc_too_slow <= .50), by = subj_id]$subj_id

# Checks if log reaction time is below M-3SD of log reaction times of learning phase in > 50% of test trials for a given test stimulus
log_times <- dt[time_pressure_cond == FALSE & block == "training" & time > 0, .(mean_log_learn = mean(log(time)),
                                                                                sd_log_learn = sd(log(time))), by = list(subj_id)]
log_times <- merge(dt[time_pressure_cond == FALSE & block == "test" & time > 0, .(log_test = log(time)), by = list(subj_id, stim)], log_times)
include <- c(valid_subj_ids, log_times[, .(perc_too_slow = mean(log_test < mean_log_learn - 3*sd_log_learn)), by = list(subj_id, stim)][, all(perc_too_slow <= .50), by = subj_id]$subj_id)

# Drops invalid subj_ids
dt <- dt[subj_id %in% include, ]

# Participants who stated task was rather unclear or absolutely unclear
exclude <- c(4, 9, 16, 32, 34, 40, 46, 54) # had to be done by hand as subject number was not automatically saved in questionnaire
dt <- dt[!subj_id %in% exclude]

# Prepares demographics data
demographics <- demographics[!1:2, !3:17]
demographics <- demographics[!grep("un", task_clear)]

# Changes subject ids to go from 1 to x
ids <- dt[, .(subj_id = sort(unique(subj_id)))][, .(subj_id, subj = order(subj_id))][, .(subj_id, subj = paste0("s", str_sub(paste0("00", subj), -3, -1)))]
dt <- merge(dt, ids)

# Old or novel stimuli
dt[, phase := ifelse(block == "training", "learning", "transfer")]
dt[, type := ifelse(stim %in% unique(stim[block == "training"]), "old", "new")]
dt <- merge(dt, unique(dt[!(phase == "transfer" & type == "old"), .(stim, c = true_cat)]), by = "stim", sort = FALSE)
dt[, stimulus := paste0(stimulus, "-", c)]

# Experimental condition
dt[, cond := ifelse(time_pressure_cond == TRUE, "tp", "no_tp")]

# Category-label association
dt[, c_key := ifelse(cat0_is_right_key == 1, "0-right", "0-left")]
dt <- dt[, .(subj, cond, phase, trial, type, stimulus, f1, f2, f3, c, resp = response, time, tp = time_pressure_in_ms, col_order = color_presentation_order, f_order = feature_presentation_order, c_key)]

# fwrite(dt, "../../../data/processed/categorization.csv")
