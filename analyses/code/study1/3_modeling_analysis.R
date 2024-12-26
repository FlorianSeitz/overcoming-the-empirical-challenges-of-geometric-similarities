# ==========================================================================
# Analyses Study 1: Analyses the results of the cognitive model fitting
# ==========================================================================

# ==========================================================================
# Parameters that need to be specified
refit_models <- FALSE # specify if you want to refit the models
models <- c("gcm", "gcm_disc", "gcm_unidim", "gcm_disc_unidim", "rule_seitz2021") # specify the models to analyze
crit <- .90 # pairwise model comparison criterium
assign <- .70 # model assignment threshold
exploratory <- FALSE # if false makes main analyses, if true makes exploratory analyses (with rule_seitz2021)
# ==========================================================================

# ==========================================================================
# Prepares packages and data
# ==========================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, cognitiveutils, RVAideMemoire, pwr, esc, splitstackshape, gtools, lme4, ez)
devtools::load_all("~/R/cognitivemodels/")
# pacman::p_load_gh("FlorianSeitz/cognitivemodels")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # change wd to where this script lies (RStudio needed!)
source("setup_models.R")
source("utils_ranking.R")
d <- fread("../../../data/processed/categorization.csv", key = "subj")[type == "new"]
ifelse(refit_models, 
       source("2_cognitive_modeling.R"), 
       fits <- readRDS("../../other/categorization_cognitive_models.rds"))

# ==========================================================================
# Gives a summary of parameter estimates
# ==========================================================================
pars <- fits[model %in% models, V2[[1]]$parm, by = list(subj, model)]
pars_long <- melt(pars, id.vars = c("subj", "model"), variable.name = "par")
pars_long[!grepl("unidim", model), .(m = mean(value), md = median(value), sd = sd(value)), by = par][, round(.SD, 2), by = par]
pars_long[grepl("unidim", model) & grepl("w", par), .(sum = sum(value)), by = .(model, par)]

# ==========================================================================
# Makes model predictions
# ==========================================================================
preds <- fits[model %in% models, V2[[1]]$predict(newdata = d[subj]), by = .(subj, model)]
preds[, trial := d[subj, "trial"], by = .(subj, model)]
preds <- as.data.table(tidyr::pivot_wider(preds, names_from = "model", values_from = "V1"))
preds[, random := 0.5] # prediction of random model
d <- merge(d, preds, by = c("subj", "trial"))

# ==========================================================================
# Calculates the models' goodness-of-fit (gof, i.e., log likelihood)
# ==========================================================================
if(exploratory) { # includes random model in analyses
  models <- c(models[models != "random"], "random") 
} else { # exludes rule_seitz2021 model, but includes random model in analyses
  models <- c(models[!models %in% c("rule_seitz2021", "random")], "random") 
}

gofs <- d[, lapply(.SD, function(x) {
  gof(obs = resp, pred = x, type = "loglik", response = "disc", na.rm = TRUE)
}), .SDcols = models, by = .(subj, cond)]

# ==========================================================================
# Performs analyses: calculates weights, makes pairwise model comparisons
# ==========================================================================
gofs[, best_model := names(which.max(.SD)), by = .(subj, cond)]
weights <- gofs[, exp(.SD - max(.SD)), .SDcols = models, by = .(subj, cond, best_model)]
weights[, (models) := .SD/ rowSums(.SD), .SDcols = models] # equal to AIC weights

# computes model ranks (see Table 4)
ranking(weights, models, crit, exclude = TRUE)
ranking(weights[, lapply(.SD, mean), .SDcols = models, by = .(cond)], models, crit)

# computes number of participants best described by each model (see Appendix Table C1)
model_distr <- weights[rowSums((weights > assign)[, models]) > 0, .N, by = .(best_model, cond)]
model_distr[, chisq.test(N), by = cond] # chi-squared test
model_distr[cond == "tp", multinomial.multcomp(N, p.method = "holm")]
model_distr[cond == "no_tp", multinomial.multcomp(N, p.method = "holm")]

model_distr_wide <- dcast(model_distr, cond ~ best_model)
fisher.test(model_distr_wide[, -1])
disc_mink <- model_distr[best_model != "random", .(N = sum(N)), by = .(cond, disc = grepl("disc", best_model))]
fisher.test(dcast(disc_mink, cond ~ disc)[, -1])

# ==========================================================================
# Additional fit indices (see Table 7)
# ==========================================================================
d_long <- melt(d, measure.vars = list(models), variable.name = "model", value.name = "pred")
add_gofs <- d_long[!is.na(resp), .(mape = MAPE(obs = resp, pred = pred, response = 'disc'),
                                   arg = mean(choicerule(x = pred, type = "argmax")[, 1] == resp),
                                   mse = mean((resp - pred)^2, na.rm = T),
                                   ll = gof(obs = resp, pred = pred, type = "loglik", response = 'disc', na.rm = TRUE)), by = .(model, subj, cond)]
add_gofs <- melt(add_gofs, measure.vars = c("mape", "arg", "mse", "ll"), variable.name = "gof")
add_gofs[, .(m = mean(value), 
             md = median(value), 
             sd = sd(value)), by = .(model, gof, cond)][, round(.SD, 2), by = .(model, gof, cond)]

# ==========================================================================
# Choice inconsistency analyses (exploratory)
# ==========================================================================
d2 <- fread("../../../data/processed/categorization.csv", key = "subj")[phase != "learning"]
models <- models[grepl("gcm", models)]
fits[model %in% models, V2[[1]]$choicerule <- NULL, by = .(subj, model)] # removes choice rule
preds_no_cr <- fits[model %in% models, V2[[1]]$predict(newdata = d2[subj]), by = .(subj, model)] # makes predictions (V1)
preds_no_cr[, c("phase", "type", "resp") := d2[subj, .(phase, type, resp)], by = .(subj, model)]

pacman::p_load(doFuture)
registerDoFuture()
plan(multisession, workers = 4L)  ## on MS Windows
setkey(preds_no_cr, "subj")  
tau <- foreach(x = unique(preds_no_cr$subj),
               .combine = "rbind",
               .inorder = FALSE, 
               .packages = c("data.table", "cognitivemodels", "devtools")) %dopar% {
                 devtools::load_all("~/R/cognitivemodels/")
                 preds_no_cr[.(x), .(
                   test = softmax(resp ~ V1, data = .SD, options = list(solver = "solnp", lb = c(tau = .1)))$par), by = .(subj, model)]
               }   
tau[, test := as.double(test)]
tau <- merge(weights[weights[, rowSums(.SD > assign) > 0, .SDcols = models], .(subj, cond, model = best_model)], tau)
tau <- merge(tau, pars[, .(subj, model, training = tau)])
tau[, t.test(log(test)[cond == "tp"], log(test)[cond == "no_tp"])]
tau <- as.data.table(pivot_longer(tau, cols = c("training", "test"), names_to = "phase"))
tau[, .(m = mean(value), md = median(value), sd = sd(value)), by = .(cond, phase)][, round(.SD, 2), by = .(cond, phase)]
