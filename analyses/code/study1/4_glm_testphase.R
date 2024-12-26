# ==========================================================================
# Analyses Study 1: Runs a general linear model with logit link
# ==========================================================================

# ==========================================================================
# Prepares packages and data
# ==========================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, lme4, emmeans, stringr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # change this to the folder where this script is lying; works only in RStudio

dt <- fread("../../../data/processed/categorization.csv")[type == "new" & !is.na(resp)]
dt[, subj := as.factor(subj)]
dt[, cond := factor(cond, levels = c("tp", "no_tp"), labels = c("tp", "no_tp"))]
dt[, test_stim_type := as.factor(ifelse(stimulus %in% c("003-NA", "100-NA"), str_sub(stimulus, 1, 3), "221-331"))]

contrasts(dt$test_stim_type) <- contr.sum
contrasts(dt$cond) <- contr.sum

# ==========================================================================
# Performs generalized linear modelling (see Table 3)
# ==========================================================================
full_model <- glmer(resp ~ cond * test_stim_type + (1|subj), data = dt, family = binomial)
round(summary(full_model)$coefficients, 2) # Table 3

restricted_model <- glmer(resp ~ cond + test_stim_type + (1|subj), data = dt, family = binomial)

anova(restricted_model, full_model) # tests for necessity of interaction term

# ==========================================================================
# Makes pairwise comparisons
# ==========================================================================
emm1 <- emmeans(full_model, ~ test_stim_type | cond)
pairs(emm1, adjust = "holm")

emm2 <- emmeans(full_model, ~ cond |  test_stim_type)
pairs(emm2, adjust = "holm", type = "response") 