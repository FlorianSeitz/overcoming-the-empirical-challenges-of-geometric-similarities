# ==========================================================================
# Setup: Cognitive models
# Author: Florian I. Seitz
# ==========================================================================

## Set up the models -------------------------------------------------------
# GCM: standard generalized context model (with Minkowski similarity)
# GCM_disc: generalized context model with the discrete similarity
# GCM with fixed bias (b0 = b1), metric exponent (r), and similarity exponent (q); first block discounted (dc)

GCM <- function(dt, metr = "minkowski", fix = list(r = 1, p = "r"), discount = nrow(dt) - 100) {
  gcm(formula = resp ~ f1 + f2 + f3, 
      cat = ~ c, 
      data = dt, 
      discount = discount, 
      metric = metr, 
      fix = fix, 
      choicerule = "softmax"
  )
}

# GCM with discrete metric
GCM_disc <- function(dt) {
  GCM(dt = dt, metr = "discrete")
}

# unidimensional GCM
GCM_unidim <- function(dt, metr = "minkowski", fix) {
  gcm_unidim(formula = resp ~ f1 + f2 + f3, 
             cat = ~ c, 
             data = dt, 
             discount = nrow(dt[phase == "learning"]), 
             metric = metr, 
             fixed = fix, 
             choicerule = "softmax"
  )
}

# unidimensional GCM with discrete metric
GCM_disc_unidim <- function(dt, fix) {
  GCM_unidim(dt = dt, metr = "discrete", fix = fix)
}

# Decision tree (attend to second feature first, then to one of the remaining features)
RULE_seitz2021 <- function(dt, discount = nrow(dt) - 100) {
  rule_seitz2021(resp ~ f1 + f2 + f3, 
                 data = dt,
                 choicerule = "epsilon",
                 discount = discount
  )
}