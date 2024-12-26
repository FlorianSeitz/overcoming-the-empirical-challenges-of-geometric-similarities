# ==========================================================================
# Analyses Study 3: Runs a general linear model with logit link
# ==========================================================================

# ==========================================================================
# Prepares packages and data
# ==========================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, lme4, lmerTest, emmeans, pbkrtest, multcomp, stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # change this to the folder where this script is lying; works only in RStudio

dt <- fread("../../../data/processed/similarity.csv")[type %in% c("one-large", "all-small", "two-large", "two-small")]

dt[, subj := as.factor(subj)]
dt[, time_pressure := ifelse(phase == "no_tp", "none", str_sub(cond, -2, -1))]
dt[, time_pressure := as.factor(time_pressure)]
dt[, type := factor(type, levels = c("two-large", "one-large", "all-small", "two-small"), labels = c("two-large", "one-large", "all-small", "two-small"))]

dcast(dt[, mean(resp, na.rm = T), by = .(time_pressure, type, subj)][, round(mean(V1, na.rm = T), 2), by = .(time_pressure, type)], type ~ time_pressure)
dcast(dt[, sd(resp, na.rm = T), by = .(time_pressure, type, subj)][, round(mean(V1, na.rm = T), 2), by = .(time_pressure, type)], type ~ time_pressure)

# ==========================================================================
# Performs generalized linear model on raw responses
# ==========================================================================
contrasts(dt$time_pressure) <- contr.treatment(4, base = 1)
contrasts(dt$type) <- contr.treatment(4, base = 1)

full_model <- lmer(resp ~ time_pressure * type + (1|subj), data = dt) 
restricted_model <- lmer(resp ~ time_pressure + type + (1|subj), data = dt) 
anova(restricted_model, full_model, refit = FALSE)

# main effects
emm1 <- emmeans(full_model, ~ type | time_pressure, pbkrtest.limit = 8206)
pairs(emm1, adjust = "holm") # all pairwise comparisons 

emm2 <- emmeans(full_model, ~ time_pressure | type, pbkrtest.limit = 8206)
pairs(emm2, adjust = "holm") # all pairwise comparisons 

# interactions
contrasts <- contrast(emm1, interaction = c(type = "pairwise", type = "pairwise"), by = NULL, adjust = "holm")
rel_contrasts <- summary(contrasts[c(6, 12, 18, 29, 35, 41, 51, 57, 63, 72, 78, 84, 92, 98, 104, 111, 117, 123, 129, 135, 146, 152, 162, 168, 177, 183, 191, 197, 204, 210, 216, 227, 237, 246, 254, 261), ])

# ==========================================================================
# Performs generalized linear model on response variability (see Table 9)
# ==========================================================================
dt_sub <- dt[, .(sd = sd(resp, na.rm = T)), by = .(subj, time_pressure, type)]
contrasts(dt_sub$time_pressure) <- contr.treatment(4, base = 4)
contrasts(dt_sub$type) <- contr.treatment(4, base = 4)

m1 <- dt_sub[, lmer(sd ~ type * time_pressure + (1|subj))]
m2 <- dt_sub[, lmer(sd ~ type + time_pressure + (1|subj))]
anova(m1, m2, refit = FALSE)

round(summary(m2)$coefficients, 2)

emm <- emmeans(m2, ~ time_pressure | type)
pairs(emm, adjust = "holm") # all pairwise comparisons --> checks for choice inconsistency

# ==========================================================================
# Performs generalized linear model on closeness to scale center
# ==========================================================================
dt <- fread("../../../data/processed/similarity.csv")
dt[type == "all-different", type := paste0(type, "-", abs(as.numeric(str_sub(stimulus, 1, 3)) - as.numeric(str_sub(stimulus, -3, -1)))/111)]
dt[, time_pressure := ifelse(phase == "no_tp", "none", str_sub(cond, -2, -1))]
dt[, subj := as.factor(subj)]
dt[, time_pressure := as.factor(time_pressure)]
dt[, condition_is_15 := factor(time_pressure == "15", labels = c("TRUE", "FALSE"), levels = c("TRUE", "FALSE"))]
dt[, type := factor(type, levels = c("two-large", "one-large", "all-small", "two-small", "all-equal", "all-different-1", "all-different-3", "all-different-4"), labels = c("two-large", "one-large", "all-small", "two-small", "all-equal", "all-different-1", "all-different-3", "all-different-4"))]
dt[, mean := 1 - abs(resp - .5)]
contrasts(dt$time_pressure) <- contr.treatment(4, base = 4)
contrasts(dt$type) <- contr.treatment(8, base = 4)

full_model <- lmer(mean ~ time_pressure * type + (1|subj), data = dt) 
restricted_model <- lmer(mean ~ time_pressure + type + (1|subj), data = dt) 
anova(restricted_model, full_model, refit = FALSE)

emm2 <- emmeans(full_model, ~ time_pressure | type, pbkrtest.limit = 8206)
pairs(emm2, adjust = "holm") # all pairwise comparisons 