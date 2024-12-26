# ==========================================================================
# Setup: Cognitive models
# Author: Florian I. Seitz
# ==========================================================================
# library(cognitivemodels) # from janajarecki/github


## Set up the models -------------------------------------------------------
# GCM: standard generalized context model (with Minkowski similarity)
GCM <- function(dt, comb, ll_bool, disc = FALSE) {
  m <- solnp(pars = c(w1 = 1/3, w2 = 1/3, w3 = 1/3, c = 2.5, sigma = .5), 
             fun = fit_fun, 
             LB = c(w1 = 0, w2 = 0, w3 = 0, c = .0001, sigma = 0), 
             UB = c(w1 = 1, w2 = 1, w3 = 1, c = 5, sigma = 1), 
             eqfun = eq_fun, 
             eqB = 1,
             d = dt,
             disc = disc,
             ll_bool = ll_bool)
  return(list(list(w1 = m$pars[["w1"]],
                   w2 = m$pars[["w2"]],
                   w3 = m$pars[["w3"]],
                   c = m$pars[["c"]],
                   sigma = m$pars[["sigma"]],
                   ll_train = -1 * tail(m$values, 1),
                   ll_test = -1 * fit_fun(m$pars, dt, (1-ll_bool) * as.numeric(!dt$type %in% c("all-equal", "all-different")), disc = disc),
                   ll_test_filler = -1 * fit_fun(m$pars, dt, 1-ll_bool, disc = disc)),
              list(type = dt$type,
                   stimulus = dt$stimulus,
                   preds = fit_fun(m$pars, dt, predict = TRUE, disc = disc),
                   obs = dt$resp,
                   bool_train = ll_bool)))
}

# GCM with discrete metric
GCM_disc <- function(dt, comb, ll_bool) {
  GCM(dt = dt, comb = comb, ll_bool = ll_bool, disc = TRUE)
}

# unidimensional GCM
GCM_unidim <- function(dt, comb, ll_bool, disc = FALSE) {
  gof <- Inf
  for(i in 1:3) {
    m <- solnp(pars = c(c = 2.5, sigma = .5), 
               fun = fit_fun, 
               LB = c(c = .0001, sigma = 0), 
               UB = c(c = 5, sigma = 1), 
               d = dt,
               disc = disc,
               unidim = i,
               ll_bool = ll_bool)
    if(tail(m$values, 1) < gof) {
      par <- c(diag(3)[i, ], m$pars); names(par)[1:3] <- paste0("w", 1:3)
      gof <- tail(m$values, 1)
    }
  }
  return(list(list(w1 = par[["w1"]],
                   w2 = par[["w2"]],
                   w3 = par[["w3"]],
                   c = par[["c"]],
                   sigma = par[["sigma"]],
                   ll_train = -1 * tail(m$values, 1),
                   ll_test = -1 * fit_fun(par, dt, (1-ll_bool) * as.numeric(!dt$type %in% c("all-equal", "all-different")), disc = disc),
                   ll_test_filler = -1 * fit_fun(par, dt, 1-ll_bool, disc = disc)),
              list(type = dt$type,
                   stimulus = dt$stimulus,
                   preds = fit_fun(par, dt, predict = TRUE, disc = disc),
                   obs = dt$resp,
                   bool_train = ll_bool)))
}

# unidimensional GCM with discrete metric
GCM_disc_unidim <- function(dt, comb, ll_bool) {
  GCM_unidim(dt = dt, comb = comb, ll_bool = ll_bool, disc = TRUE)
}
