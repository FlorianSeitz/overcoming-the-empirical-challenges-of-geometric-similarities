# ==========================================================================
# Analyses Study 2: Preprocesses data
# ==========================================================================

rm(list = ls(all = TRUE))

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # changes working directory to folder where this script is
demo_15_30 <- fread("../../../data/raw/main_study3_demographics/main_study3_demographics.csv")
demo_50 <- fread("../../../data/raw/main_study2_demographics/main_study2_demographics.csv")

# ==========================================================================
# Defines to-be-excluded participants who said task was unclear
# ==========================================================================
exclude <- c(demo_15_30[grepl("un", task_clear), as.numeric(sub(".*=", "", Q_URL))],
             demo_50[grepl("un", task_clear), as.numeric(sub(".*=", "", Q_URL))] + 200)

# ==========================================================================
# Reads data from 50% condition
# ==========================================================================
files <- list.files(path = "../../../data/raw/main_study2_experiment/", pattern = "", full.names = TRUE)
dt50 <- rbindlist(lapply(files, fread, fill = TRUE, colClasses = list("character" = c("type", "color_presentation_order", "shape_presentation_order"))))
dt50[, cond := "50"]
dt50[, subj_id := subj_id + 200]

# ==========================================================================
# Reads data from 30% condition
# ==========================================================================
files <- list.files(path = "../../../data/raw/main_study3b_experiment/", pattern = "", full.names = TRUE)
dt30 <- rbindlist(lapply(files, fread, fill = TRUE, colClasses = list("character" = c("type", "color_presentation_order", "shape_presentation_order"))))
dt30[, cond := "30"]

# ==========================================================================
# Reads data from 15% condition
# ==========================================================================
files <- list.files(path = "../../../data/raw/main_study3a_experiment/", pattern = "", full.names = TRUE)
dt15 <- rbindlist(lapply(files, fread, fill = TRUE, colClasses = list("character" = c("type", "color_presentation_order", "shape_presentation_order"))))
dt15[, cond := "15"]

# ==========================================================================
# Combines conditions and performs preprocessing
# ==========================================================================
dt <- rbind(dt50, dt30, dt15)
dt <- dt[!subj_id %in% exclude]

# Extracts feature values
dt[, paste0("f", 1:3, "l") := tstrsplit(stim_left, "")]
dt[, paste0("f", 1:3, "r") := tstrsplit(stim_right, "")]

# Rescales responses to lie between 0 and 1
dt[!is.na(response), resp := (response - min(response))/max(response - min(response)), by = subj_id]

dt[, trial := 1:.N, by = subj_id]
dt[, phase := ifelse(block == "familiarization", "no_tp",
                     ifelse(block == "familiarization_time_pressure", "tp",
                            ifelse((time_pressure_first == TRUE & block == "test_1") | (time_pressure_first == FALSE & block == "test_2"), "tp", "no_tp")))]
dt[, type := ifelse(type == "A", "one-large", ifelse(type == "B", "all-small", ifelse(type == "C", "two-large", ifelse(type == "D", "two-small", ifelse(type == "I", "all-equal", ifelse(type == "V", "all-different", "familiarization"))))))]

ids <- dt[, .(subj_id = sort(unique(subj_id))), by = cond][, .(subj_id, subj = order(subj_id)), by = cond][, .(cond, subj_id, subj = paste0("s", ifelse(cond == 50, subj + 200, ifelse(cond == 30, subj + 100, str_sub(paste0("00", subj), -3, -1)))))]
dt <- merge(dt, ids)

dt <- dt[, .(subj, cond = paste0("tp_", cond), phase, trial, type, stimulus = paste0(stim_left, "-", stim_right), f1l, f2l, f3l, f1r, f2r, f3r, raw = response, resp, contemplation_time = contemplation_time, contemplation_tp = contemplation_time_pressure, resp_time = reaction_time, resp_tp = reaction_time_pressure, col_order = color_presentation_order, shape_order = shape_presentation_order, tp_first = time_pressure_first)]
# fwrite(dt, "../../../data/processed/similarity.csv")

demo_15_30[, subj_id := as.numeric(sub(".*=", "", Q_URL))]
demo_50[, subj_id := as.numeric(sub(".*=", "", Q_URL)) + 200]
demographics <- rbind(demo_15_30, demo_50)
demographics <- merge(demographics, ids)
demographics <- demographics[, c(33, 19:26)]
# fwrite(demographics, "../../../data/processed/similarity_demographics.csv")
