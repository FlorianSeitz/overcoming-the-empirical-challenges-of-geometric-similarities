# ==========================================================================
# Analyses Study 3: Analyses the results of the cognitive model fiting
# ==========================================================================

# ==========================================================================
# Parameters that need to be specified
models <- c("gcm", "gcm_disc", "gcm_unidim", "gcm_disc_unidim") # specify the models to analyze
crit <- .90 # pairwise model comparison criterion
assign <- .70 # model assignment threshold
# ==========================================================================

# ==========================================================================
# Prepares packages and data
# ==========================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, emmeans, cognitiveutils, RVAideMemoire, pwr, esc, splitstackshape, gtools, lme4, ez)
devtools::load_all("~/R/cognitivemodels/")
# pacman::p_load_gh("JanaJarecki/cognitivemodels")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # change wd to where this script lies (RStudio needed!)
source("setup_models.R")

d <- fread("../../../data/processed/similarity.csv")[type != "famialiarization"]
res15 <- fread("../../other/similarity_pars_likelihoods_cond15.csv")
res30 <- fread("../../other/similarity_pars_likelihoods_cond30.csv")
res50 <- fread("../../other/similarity_pars_likelihoods_cond50.csv")
res <- rbind(res15, res30, res50)

# ==========================================================================
# Extracts the models' goodness-of-fit (gof, i.e., log likelihood)
# ==========================================================================
res <- res[, .(sigma = median(sigma), ll = median(ll_test)), by = .(subj, cond = time_pressure, model)]
gofs <- dcast(res, subj + cond ~ model, value.var = "ll")
gofs[, random := 0] # because density for a uniform model is 1

models <- c(models[models != "random"], "random") # includes random model in analyses

# ==========================================================================
# Individual analyses: calculates weights, makes pairwise model comparisons
# ==========================================================================
gofs[, best_model := names(which.max(.SD)), by = .(subj, cond)]
weights <- gofs[, exp(.SD - max(.SD)), .SDcols = models, by = .(subj, cond, best_model)]
weights[, (models) := .SD/ rowSums(.SD), .SDcols = models] # equal to AIC weights

# computes number of participants best described by each model (see Appendix Table C1)
model_distr <- weights[rowSums((weights > assign)[, models]) > 0, .N, by = .(best_model, cond)]

# defines which participants are described by which model in which time pressure condition
described <- weights[rowSums((weights > assign)[, models]) > 0, .(subj, cond, model = best_model)]
described_pars <- res[described[model != "random"], , on = c("subj", "cond", "model")]
described_pars[, cond := factor(cond, levels = c("none", "50", "30", "15"))]

# ==========================================================================
# Runs linear model on parameter estimates for sigma
# ==========================================================================
model <- lmer(sigma ~ cond + (1|subj), data = described_pars) 
emm2 <- emmeans(model, ~ cond, pbkrtest.limit = 8206)
pairs(emm2, adjust = "holm") # all pairwise comparisons 
