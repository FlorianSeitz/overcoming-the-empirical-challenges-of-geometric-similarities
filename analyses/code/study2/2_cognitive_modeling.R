# ==========================================================================
# Analyses Study 2: Fits cognitive models
# ==========================================================================

# ==========================================================================
tp_cond <- "50" # specify time pressure condition (only 50 in Study 2)
# ==========================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, modelr, purrr, devtools, cognitiveutils, doRNG, Rsolnp, stringr)

parallel <- TRUE # fit on a parallel machine (Unix) or single core
if (parallel == TRUE) pacman::p_load(doFuture)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
d <- fread("../../../data/processed/similarity.csv")[!type == "familiarization"]
d[, type_trial := 1:.N, by = .(paste0(type, trial), subj, phase)]
d[, time_pressure := ifelse(phase == "no_tp", "none", str_sub(cond, -2, -1))]
d <- d[cond == paste0("tp_", tp_cond)]

# ==========================================================================
# Makes objective function and equality constraints for attention weights
# ==========================================================================
fit_fun <- function(pars, d, ll_bool = NULL, predict = FALSE, disc = FALSE, unidim = NULL) {
  stim_l <- d[, paste0("f", 1:3, "l")]
  stim_r <- d[, paste0("f", 1:3, "r")]
  
  diff <- abs(stim_l - stim_r)
  if(disc) diff <- stim_l != stim_r
  
  w <- pars[grepl("^w", names(pars))]
  if(!is.null(unidim)) w <- diag(3)[unidim, ]
  
  pred <- exp(-pars["c"] * c(w %*% t(diff)))
  if(predict) return(pred)
  obs <- d$resp
  ll <- dnorm(x = obs, mean = pred, sd = pars["sigma"], log = TRUE)
  ll[is.na(ll)] <- 0
  return(-sum(ll_bool * ll))
}

eq_fun <- function(pars, d, ll_bool, predict, disc, unidim) {
  return(sum(pars[grepl("^w", names(pars))]))
}

is_fit <- function(x, comb) {
  return(as.numeric(x[, (type == "one-large" & type_trial %in% comb[, 1]) | 
                        (type == "all-small" & type_trial %in% comb[, 2]) |
                        (type == "two-large" & type_trial %in% comb[, 3]) |
                        (type == "two-small" & type_trial %in% comb[, 4])]))
}

# ==========================================================================
# Specifies to-be-fitted models
# ==========================================================================
source("setup_models.R")
model_list <- list(
  gcm = GCM,                          # gcm (Minkowski)
  gcm_disc = GCM_disc,                # gcm (Discrete)
  gcm_unidim = GCM_unidim,            # gcm (Minkowski, unidimensional)
  gcm_disc_unidim = GCM_disc_unidim   # gcm (Discrete, unidimensional)
)

# ==========================================================================
# Makes cross-validation data sets cv_data for each participant
# ==========================================================================
combs <- combn(1:4, 2) # which stimuli of a given type to pick for fitting
i_combs <- as.matrix(expand.grid(1:ncol(combs), 1:ncol(combs), 1:ncol(combs), 1:ncol(combs))) # indices for combs for each type A-D
its <- expand.grid(id = d[, unique(time_pressure)], cv = 1:nrow(i_combs))
colnames(its)[1] <- c("time_pressure")

# ==========================================================================
# Fits models
# ==========================================================================
start_time <- proc.time()
if (parallel == TRUE) {
  registerDoFuture()
  plan(multisession, workers = 4L)  ## on MS Windows
  fits <- foreach(x = unique(d$subj),
                  .combine = "rbind",
                  .inorder = FALSE, 
                  .packages = c("data.table", "modelr", "devtools", "Rsolnp"),
                  .export = c("model_list", "its", "i_combs", "combs", 
                              "GCM", "GCM_disc", "GCM_unidim", "GCM_disc_unidim",
                              "fit_fun", "eq_fun")) %dorng% {
                    source("setup_models.R", local = TRUE)
                    curr_its <- as.data.table(its)

                    do.call("rbind", lapply(1:nrow(curr_its), function (i) {
                      it <- curr_its[i, ]
                      i_comb <- i_combs[it$cv, ]
                      comb <- combs[, i_comb]
                      curr_d <- d[(subj == x) & (time_pressure == it$time_pressure), ]
                      ll_bool <- is_fit(x = curr_d, comb = comb)
                      curr_d[, .(
                        cv = paste0(i_comb, collapse = ""),
                        model = names(model_list),
                        res = map(model_list, exec, dt = .SD, comb = comb, ll_bool = ll_bool)), by = .(subj, time_pressure)]
                    }))
                  }
} 
stop_time <- proc.time()
timetaken(start_time)

saveRDS(fits, gsub("X", tp_cond, "../../other/similarity_cognitive_models_X.rds"))
fwrite(fits[, res[[1]][[1]], by = .(subj, time_pressure, cv, model)], gsub("X", tp_cond, "../../other/similarity_pars_likelihoods_condX.csv"))

# ==========================================================================
# Access model fitting results the following way:
# ==========================================================================
# fits[, res[[1]][[1]], by = .(subj, time_pressure, cv, model)] # parameters + gof
# fits[, res[[1]][[2]], by = .(subj, time_pressure, cv, model)] # predictions
