# ==========================================================================
# Analyses Study 2: Runs a general linear model with identity link
# ==========================================================================

# ==========================================================================
# Prepares packages and data
# ==========================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, lme4, lmerTest, emmeans, stringr, pbkrtest)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # change this to the folder where this script is lying; works only in RStudio

dt <- fread("../../../data/processed/similarity.csv")[cond == "tp_50"][grep("large|small", type)]
dt[, subj := as.factor(subj)]
dt[, phase := factor(phase, levels = c("tp", "no_tp"), labels = c("tp", "no_tp"))]
dt[, type := factor(type, levels = c("two-large", "one-large", "all-small", "two-small"), labels = c("two-large", "one-large", "all-small", "two-small"))]

contrasts(dt$type) <- contr.sum
contrasts(dt$phase) <- contr.sum

# ==========================================================================
# Performs generalized linear modelling (see Table 6)
# ==========================================================================
full_model <- lmer(resp ~ phase * type + (1|subj), data = dt) 
restricted_model <- lmer(resp ~ phase + type + (1|subj), data = dt) 

anova(restricted_model, full_model, refit = FALSE)

round(summary(full_model)$coefficients, 2) # equals estimates in Table 6

# ==========================================================================
# Makes pairwise comparisons
# ==========================================================================
emm1 <- emmeans(full_model, ~ type |  phase, pbkrtest.limit = 8206)
pairs(emm1, adjust = "holm") # all pairwise comparisons 

# ==========================================================================
# Makes exploratory standard deviation analyses
# ==========================================================================
dt_sd <- dt[, .(sd = sd(resp, na.rm = T)), by = .(subj, phase, type)]
dt_sd <- dcast(dt_sd, formula = subj + type ~ phase, value.var = "sd")
dt_sub[, t.test(x = tp, y = no_tp, paired = TRUE)]
