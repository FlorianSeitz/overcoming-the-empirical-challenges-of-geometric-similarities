# ==========================================================================
# Analyses Experiment: 2. Fits cognitive models
# ==========================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, modelr, purrr, devtools, cognitiveutils)
devtools::load_all("~/R/cognitivemodels/")
parallel <- TRUE # fit on a parallel machine (Unix) or single core
if (parallel == TRUE) pacman::p_load(doFuture)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
d <- fread("../../../data/processed/categorization.csv")[phase == "learning"]

# ==========================================================================
# Specifies to-be-fitted multidimensional models
# ==========================================================================
source("setup_models.R")
model_list <- list(
  gcm = GCM,           # multidimensional gcm (Minkowski similarity)
  gcm_disc = GCM_disc  # multidimensional gcm (discrete similarity)
)

# Opt-in: appends exploratory decision tree (first feature 2, then feature 1 or 3)
# model_list <- c(model_list, rule_seitz2021 = RULE_seitz2021)

# ==========================================================================
# Fits models
# ==========================================================================
if (parallel == TRUE) {
  # registerDoMC(cores = detectCores())
  registerDoFuture()
  plan(multisession, workers = 4L)  # on MS Windows
  setkey(d, "subj")  
  fits <- foreach(x = unique(d$subj),
                  .combine = "rbind",
                  .inorder = FALSE, 
                  .packages = c("data.table", "modelr", "devtools"),
                  .export = c("model_list", "GCM", "GCM_disc", "RULE_seitz2021")) %dopar% {
                    devtools::load_all("~/R/cognitivemodels/")
                    source("setup_models.R", local = TRUE)
                    d[.(x), .(
                      model = names(model_list),
                      map(model_list, exec, dt = .SD)), by = subj]
                  }   
} else {
  fits <- d[, .(
    model = names(model_list),
    fit = map(model_list, exec, dt = .SD)),
    by = subj]
}

# saveRDS(fits, "../../other/categorization_cognitive_models.rds")

# ==========================================================================
# Specifies to-be-fitted unidimensional models and to-be-fixed parameters
# ==========================================================================
d <- fread("../../../data/processed/categorization.csv")[!is.na(resp) & type != "new"]
# fits <- readRDS("../../other/categorization_cognitive_models.rds")
pars <- fits[model == "gcm", V2[[1]]$parm, by = subj]

model_list <- list(
  gcm_unidim = GCM_unidim,           # unidimensional gcm (Minkowski similarity)
  gcm_disc_unidim = GCM_disc_unidim  # unidimensional gcm (discrete similarity)
)

# ==========================================================================
# Fits unidimensional models
# ==========================================================================
skip <- 69; end <- 92 # used for sourcing below (keep and ignore)
if (parallel == TRUE) {
  # registerDoMC(cores = detectCores())
  registerDoFuture()
  plan(multisession, workers = 4L)  ## on MS Windows
  setkey(d, "subj")  
  fits_unidim <- foreach(x = unique(d$subj),
                  .combine = "rbind",
                  .inorder = FALSE, 
                  .packages = c("data.table", "cognitivemodels", "cognitiveutils", "modelr", "devtools"),
                  .export = c("model_list", "pars", "GCM_unidim", "GCM_disc_unidim")) %dopar% {
                    devtools::load_all("~/R/cognitivemodels/")
                    source("setup_models.R", local = TRUE)
                    d[.(x), .(
                      model = names(model_list),
                      map(model_list, exec, dt = .SD, fix = as.list(pars[subj == .(x), 5:8]))), by = subj]
                  }   
  
} else {
  fits_unidim <- d[, .(
    model = names(model_list),
    fit = map(model_list, exec, dt = .SD)),
    by = subj]
}

# ==========================================================================
# Repeats fitting for invalid subjects (w2 = 1) using all trials
# ==========================================================================
invalid_subjs <- fits_unidim[, is.na(V2[[1]]$gofvalue), by = subj][V1 == TRUE, subj]
fits_unidim_valid <- fits_unidim[!subj %in% invalid_subjs]

d <- fread("../../../data/processed/categorization.csv")[!is.na(resp) & subj %in% invalid_subjs]
lines <- paste(scan("2-cognitive-modeling.R", what=character(), skip = skip, nlines = end-skip, sep = '\n'), collapse = "\n")
source(textConnection(lines)) # Or execute lines manually again, if it does not work.
fits_unidim <- rbind(fits_unidim_valid, fits_unidim)

# ==========================================================================
# Makes unidimensional models have only learning stimuli as exemplars
# ==========================================================================
d <- fread("../../../data/processed/categorization.csv")[phase == "learning"]
preds <- fits_unidim[, V2[[1]]$parm, by = .(subj, model)]
fits_unidim[, V2 := sapply(1:nrow(preds), function(i) {
  gcm_unidim(formula = resp ~ f1 + f2 + f3,
             cat = ~ c,
             data = d[subj == preds[i, subj]],
             metric = ifelse(preds[i, model] == "gcm_unidim", "minkowski", "discrete"),
             fixed = as.list(preds[i, 3:9]),
             choicerule = "softmax")})]


fits <- rbind(fits, fits_unidim)[order(subj)]

saveRDS(fits, "../../other/categorization_cognitive_models.rds")
