# ==========================================================================
# Analyses Study 2: Analyses the results of the cognitive model fiting
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
pacman::p_load(data.table, cognitiveutils, RVAideMemoire, pwr, esc, splitstackshape, gtools, lme4, ez)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # change wd to where this script lies (RStudio needed!)
source("utils_ranking.R")

res <- fread("../../other/similarity_pars_likelihoods_cond50.csv")
pars <- melt(res, id.vars = c("subj", "time_pressure", "cv", "model", "ll_test"), variable.name = "par")
pars[!grepl("unidim", model), median(value), by = .(subj, time_pressure, par)][, .(mean(V1), median(V1), sd(V1)), by = .(time_pressure, par)][, round(.SD, 2), by = .(time_pressure, par)]

# ==========================================================================
# Extracts the models' goodness-of-fit (gof, i.e., log likelihood)
# ==========================================================================
gofs <- res[, .(ll = median(ll_test)), by = .(subj, time_pressure, model)]
gofs[, .(m = mean(ll), md = median(ll), sd = sd(ll)), by = .(time_pressure, model)][, round(.SD, 2), by = .(time_pressure, model)]
gofs <- dcast(gofs, subj + time_pressure ~ model, value.var = "ll")
gofs[, random := 0] # because density for a uniform model is 1

models <- c(models[models != "random"], "random") # includes random model in analyses

# ==========================================================================
# Calculates weights, makes pairwise model comparisons (see Table 7)
# ==========================================================================
gofs[, best_model := names(which.max(.SD)), by = .(subj, time_pressure)]
weights <- gofs[, exp(.SD - max(.SD)), .SDcols = models, by = .(subj, time_pressure, best_model)]
weights[, (models) := .SD/ rowSums(.SD), .SDcols = models] # equal to AIC weights

# computes model ranks and mean model evidence strengths (see Table 7)
ranking(weights, models, crit, exclude = T) # ranks on individual level
ranking(weights[, lapply(.SD, mean), .SDcols = models, by = .(time_pressure)], models, crit) # ranks on aggregate level
weights[, lapply(.SD, function(x) round(mean(x), 2)), .SDcols = models, by = .(time_pressure)] # mean evidence strength

# computes number of participants best described by each model (see Appendix Table C1)
model_distr <- weights[rowSums((weights > assign)[, models]) > 0, .N, by = .(best_model, time_pressure)]
model_distr[, chisq.test(N), by = time_pressure] # chi-squared test
model_distr[time_pressure == "50", .(N = sum(N)), by = best_model == "gcm"][, multinomial.multcomp(N, p.method = "holm")]
model_distr[time_pressure == "none", .(N = sum(N)), by = best_model == "gcm"][, multinomial.multcomp(N, p.method = "holm")]
model_distr[time_pressure == "50", esc_bin_prop(N[best_model == "gcm"]/64, 64, sum(N[best_model != "gcm"])/64, 64, es.type = "or")]
model_distr[time_pressure == "none", esc_bin_prop(N[best_model == "gcm"]/64, 64, sum(N[best_model != "gcm"])/64, 64, es.type = "or")]

model_distr_wide <- dcast(model_distr, time_pressure ~ best_model)
fisher.test(model_distr_wide[, -1])

# performs exploratory choice sensitivity analyses
described <- weights[rowSums((weights > assign)[, models]) > 0, .(subj, time_pressure, model = best_model)]
described_pars <- res[described[model != "random"], , on = c("subj", "time_pressure", "model")]
described_pars <- dcast(described_pars[, median(sigma), by = .(subj, time_pressure)], subj ~ time_pressure, value.var = "V1")
colnames(described_pars)[2] <- "time_pressure"
described_pars[, t.test(time_pressure, none, paired = F)]

described_pars <- res[described[subj %in% described[model != "random", .N, by = subj][N == 2, subj]], , on = c("subj", "time_pressure", "model")]
described_pars <- dcast(described_pars[, median(sigma), by = .(subj, time_pressure)], subj ~ time_pressure, value.var = "V1")
colnames(described_pars)[2] <- "time_pressure"
described_pars[, t.test(time_pressure, none, paired = T)]

# ==========================================================================
# Additional goodness of fit indices (see Table 7)
# ==========================================================================
add_gofs <- fread("../../other/similarity_fit_measures_cond50.csv")
add_gofs <- melt(add_gofs, measure.vars = c("mape", "mae", "mse", "ll"), variable.name = "gof")
add_gofs[gof == "mape", round(mean(value), 2) / 10^13, by = .(model, gof, time_pressure)]
add_gofs[gof != "mape", round(mean(value), 2), by = .(model, gof, time_pressure)]
