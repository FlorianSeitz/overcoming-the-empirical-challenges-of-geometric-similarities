# ==========================================================================
# Simulation Category Variability - Bias and Category-Specific Sensitivity
# ==========================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, purrr, Rsolnp, doRNG)
parallel <- TRUE # fit on a parallel machine (Unix) or single core
if (parallel == TRUE) pacman::p_load(doFuture)

# ==========================================================================
set.seed(123)
n <- 100
l_max <- 5; g_max <- 5

m_l <- 0
sd_ls <- c(.01, .02, .05, .1)
m_h <- 1
sd_h <- .2
trans <- data.table(f1 = rep(seq(0, 1, .1), each = 20), c = NA)
# ==========================================================================

# ==========================================================================
# Specifies goodness-of-fit function
# ==========================================================================
fits <- function(pars, train, trans, gof = c("ll", "acc", "err", "both", "none")) {
  l1 <- pars["l1"]
  l0 <- pars["l0"]
  b1 <- pars["b1"]
  b0 <- pars["b0"]
  g1 <- pars["g1"]
  g0 <- pars["g0"]
  preds <- apply(trans, 1, function(i) {
    train[, .(s = sum(exp(-ifelse(c == 1, l1, l0) * abs(f1 - i[["f1"]])))), by = c][, 1 / (1 + (b0/b1) * ((s[c == 0])^g0/(s[c == 1])^g1))]
    # train[, .(s = sum(exp(-ifelse(c == 1, l1, l0) * abs(f1 - i[["f1"]])))), by = c][, b1 * s[c == 1]^g1 / (b1 * s[c == 1]^g1 + b0 * s[c == 0]^g0)]
  })
  
  if(gof == "none") return(preds)
  if(gof == "ll") return(-sum(dbinom(x = trans[, c], prob = preds, size = 1, log = TRUE)))
  if(gof == "acc") return(100 * mean(round(preds) == trans[, c]))
  if(gof == "err") return(100 * (1 - mean(round(preds) == trans[, c])))
  if(gof == "both") return(100 * (1 - mean(round(preds) == trans[, c])) - sum(dbinom(x = trans[, c], prob = preds, size = 1, log = TRUE)))
}

# ==========================================================================
# Specifies sampling function
# ==========================================================================
samples <- function(n, m_l, sd_l, m_h, sd_h) {
  c0 <- seq(m_l - sd_l, m_l, length.out = n)
  c1 <- seq(m_h, m_h + sd_h, length.out = n)
  
  return( data.table(f1 = c(c0, c1), c = rep(c(0, 1), each = n)) )
}

# ==========================================================================
# Specifies models (lambda and bias that can be free or constrained)
# ==========================================================================
LC_BC_GC <- function(train, trans, gof, eq = function(pars, train, trans, gof) {c(sum(pars[3:4]), diff(pars[1:2]), diff(pars[3:4]), diff(pars[5:6]))}, eb = c(1, 0, 0, 0)) {
  m <- solnp(pars = c(l1 = 1, l0 = 1, b1 = .5, b0 = .5, g1 = 1, g0 = 1), fun = fits, 
             LB = c(.001, .001, 0, 0, 0, 0), UB = c(5, 5, 1, 1, 5, 5),
             eqfun = eq, eqB = eb, train = train, trans = trans[abs(f1 - .5) > .01], gof = gof) 

  list(
    c(# m$pars,
      gof = tail(m$values, 1),
      oos = fits(pars = m$pars, train = train, trans = trans[abs(f1 - .5) < .01], gof = "ll"),
      pred = unique(fits(pars = m$pars, train = train, trans = trans[abs(f1 - .5) < .01], gof = "none")),
      resp = trans[abs(f1 - .5) < .01, mean(c)])
  )
}

LF_BC_GC <- function(train, trans, gof) {
  eq <- function(pars, train, trans, gof) {c(sum(pars[3:4]), diff(pars[3:4]), diff(pars[5:6]))}
  eb <- c(1, 0, 0)
  LC_BC_GC(train = train, trans = trans, gof = gof, eq = eq, eb = eb)
}

LC_BF_GC <- function(train, trans, gof) {
  eq <- function(pars, train, trans, gof) {c(sum(pars[3:4]), diff(pars[1:2]), diff(pars[5:6]))}
  eb <- c(1, 0, 0)
  LC_BC_GC(train = train, trans = trans, gof = gof, eq = eq, eb = eb)
}

LC_BC_GF <- function(train, trans, gof) {
  eq <- function(pars, train, trans, gof) {c(sum(pars[3:4]), diff(pars[1:2]), diff(pars[3:4]))}
  eb <- c(1, 0, 0)
  LC_BC_GC(train = train, trans = trans, gof = gof, eq = eq, eb = eb)
}

model_list <- list(
  # lc_bc_gc = LC_BC_GC,
  lf_bc_gc = LF_BC_GC,
  lc_bf_gc = LC_BF_GC,
  lc_bc_gf = LC_BC_GF
)

# ==========================================================================
# Fits models
# ==========================================================================
# res <- as.data.table(expand.grid(
#   sd_l = c(.1), 
#   true = c("lf_bc_gc", "lc_bf_gc", "lc_bc_gf"),
#   rep = 101:105))
# res[, id := as.character(1:.N)]

start_time <- proc.time()
res <- as.data.table(do.call("rbind", lapply(751:nrow(res), function(i) {
  print(timetaken(start_time))
  print(paste0("-------------- Iter: ", i, "/", nrow(res), " --------------"))
  res_i <- res[i, ]
  pars <- unlist(res_i[, .(l1, l0, b1, b0, g1, g0)])
  # pars <- c(l0 = runif(1, max = l_max), g0 = runif(1, max = g_max))
  # if (res_i$true == "lc_bc_gc") pars <- c(pars, b0 = .5, l1 = pars[["l0"]], g1 = pars[["g0"]])
  # if (res_i$true == "lf_bc_gc") pars <- c(pars, b0 = .5, l1 = runif(1, max = pars["l0"]), g1 = pars[["g0"]])
  # if (res_i$true == "lc_bf_gc") pars <- c(pars, b0 = runif(1, max = .5), l1 = pars[["l0"]], g1 = pars[["g0"]])
  # if (res_i$true == "lc_bc_gf") pars <- c(pars, b0 = .5, l1 = pars[["l0"]], g1 = runif(1, min = pars[["g0"]], max = g_max))
  # pars <- round(c(pars, b1 = 1 - pars[["b0"]]), 2)
  print(pars)
  
  train <- samples(n, m_l, res_i$sd_l, m_h, sd_h)
  preds <- fits(pars = pars, train = train, trans = trans, gof = "none")
  trans$c <- rbinom(length(preds), 1, preds)
    
  res[id == i, .(
    sd_l, true,
    fit = names(model_list),
    res = map(model_list, exec, train = train, trans = trans, gof = "ll")), by = id]
})))
timetaken(start_time)


# saveRDS(res, "~/Projects/review-geometric/category-variability-recovery.rds")
res <- res[, as.list(res[[1]][[1]]), by = .(id, sd_l, true, fit)]
res[, true := factor(true, levels = c("lc_bf_gc", "lc_bc_gf", "lf_bc_gc"))]
res[, fit := factor(fit, levels = c("lc_bf_gc", "lc_bc_gf", "lf_bc_gc"))]
res <- res[sd_l != .04]
res <- res[order(true, fit)]
res[, round(summary(resp), 2)]
res[, round(mean(resp), 2), by = true]
res[, round(mean(pred), 2), by = .(fit)]
res[, round(mean(pred), 2), by = .(true, fit)]

res[, round(mean(abs(pred - resp)), 2)]
res[, round(mean(abs(pred - resp)), 2), by = true == fit]
res[, round(mean(abs(pred - resp)), 2), by = sd_l]

res[, round(mean(oos), 2), by = .(fit)]
res[, round(mean(oos), 2), by = .(true, fit)]

res[, .SD[which.min(oos)], by = id][, mean(true == fit), by = .(true, sd_l)]
res[, .SD[diff(sort(oos)[1:2]) > .1], by = id][, .SD[which.min(oos)], by = id][, mean(true == fit), by = .(true, sd_l)]
res[, .SD[which.min(oos)], by = id][, .N, by = .(sd_l, true, fit)][, .(fit = fit, N / sum(N)), by = .(sd_l, true)][order(sd_l, true, fit)]
res[, waic := exp(-oos) / sum(exp(-oos)), by = .(id, true)]

res[, round(mean(waic), 2), by = .(fit)]
res[, round(mean(waic), 2), by = .(true, fit)][order(true, fit)]

# pars <- res[, as.list(fit[[1]][[1]]), by = .(id, sd_l, distr, rep, model)]
# pars <- melt(pars, measure.vars = c("l1", "l0", "b1", "b0", "g1", "g0"), variable.name = "par")
# pars <- dcast(pars[, mean(value), by = .(sd_l, distr, model, par)], sd_l + model ~ par + distr)
# 
# gofs <- res[, .(gof = fit[[1]][[2]], acc = fit[[1]][[3]]), by = .(id, sd_l, distr, rep, model)]
# gofs <- melt(gofs, measure.vars = c("gof", "acc"), variable.name = "var")
# gofs <- dcast(gofs[, .(mean(value)), by = .(sd_l, distr, model, var)], sd_l + model ~ var + distr)
# gofs <- merge(pars, gofs)
# 
# gofs <- gofs[, round(.SD, 2), by = .(sd_l, model)]
# gofs <- gofs[order(model, sd_l)]
# 
# gofs[, !grepl("unif", colnames(gofs)), with = FALSE][, c(1, 5:8, 3, 4, 9, 10)]
# xtable::xtable(gofs[, !grepl("unif", colnames(gofs)), with = FALSE][, c(1, 6, 5, 8, 7, 4, 3, 9, 10)], type = "latex")
# xtable::xtable(gofs[, !grepl("unif|model", colnames(gofs)), with = FALSE], type = "latex")


par_grid <- as.data.table(expand.grid(l1 = .5:4.5, l0 = .5:4.5, b1 = c(.5, .7, .9), g1 = .5:4.5, g0 = .5:4.5))
par_grid[, b0 := 1-b1]
par_grid <- par_grid[g1 >= g0 & l1 <= l0]
par_grid <- par_grid[(l1 != l0) & (b1 == b0) & (g1 == g0) | # category-specific sensitivity
                       (l1 == l0) & (b1 != b0) & (g1 == g0) | # category-specific biases
                       (l1 == l0) & (b1 == b0) & (g1 != g0) # category-specific scaling
                     ]
par_grid[, true := ifelse(l1 != l0, "lf_bc_gc", ifelse(b1 != b0, "lc_bf_gc", "lc_bc_gf"))]
par_grid[, id := 1:.N]

res <- as.data.table(expand.grid(
  id = 1:nrow(par_grid),
  sd_l = c(.002, .01, .02, .04, .1, .2)))
res <- merge(res, par_grid, sort = F)
res[, id := as.character(1:.N)]
