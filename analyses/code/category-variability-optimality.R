# ==========================================================================
# Simulation Category Variability - Bias and Category-Specific Sensitivity
# ==========================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, purrr, Rsolnp, doRNG)
parallel <- TRUE # fit on a parallel machine (Unix) or single core
if (parallel == TRUE) pacman::p_load(doFuture)

# ==========================================================================
n <- 200

m_l <- 0
m_h <- 1
sd_h <- 1
# ==========================================================================

# ==========================================================================
# Specifies goodness-of-fit function
# ==========================================================================
fits <- function(pars, d, gof = c("ll", "acc", "err", "both", "none")) {
  l1 <- pars["l1"]
  l0 <- pars["l0"]
  b1 <- pars["b1"]
  b0 <- pars["b0"]
  g1 <- pars["g1"]
  g0 <- pars["g0"]
  train <- d[["train"]]
  trans <- d[["trans"]]
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
samples <- function(n, m_l, sd_l, m_h, sd_h, distr, seed) {
  set.seed(seed)
  if (distr == "norm") {
    c0 <- rnorm(n, mean = m_l, sd = sd_l)
    c1 <- rnorm(n, mean = m_h, sd = sd_h)
  } else {
    c0 <- m_l + runif(n, -sd_l * sqrt(3), sd_l * sqrt(3))
    c1 <- m_h + runif(n, -sd_h * sqrt(3), sd_h * sqrt(3))
  }
  # return(data.table(f1 = c(c0, c1), c = rep(c(0, 1), each = n)))
  train <- data.table(f1 = c(c0[1:(n/2)], c1[1:(n/2)]), c = rep(c(0, 1), each = n/2))
  trans <- data.table(f1 = c(c0[-c(1:(n/2))], c1[-c(1:(n/2))]), c = rep(c(0, 1), each = n/2))
  return(list(train = train, trans = trans))
}

# ==========================================================================
# Specifies models (lambda and bias that can be free or constrained)
# ==========================================================================
LC_BC_GC <- function(d, gof, eq = function(pars, d, gof) {c(sum(pars[3:4]), diff(pars[1:2]), diff(pars[3:4]), diff(pars[5:6]))}, eb = c(1, 0, 0, 0)) {
  m <- solnp(pars = c(l1 = 1, l0 = 1, b1 = .5, b0 = .5, g1 = 1, g0 = 1), fun = fits, 
             LB = c(.001, .001, 0, 0, 0, 0), UB = c(5, 5, 1, 1, 5, 5),
             eqfun = eq, eqB = eb, d = d, gof = gof) 
  list(
    pars = m$pars,
    gof = tail(m$values, 1),
    acc = fits(pars = m$pars, d = d, gof = "acc")
  )
}

LF_BC_GC <- function(d, gof) {
  eq <- function(pars, d, gof) {c(sum(pars[3:4]), diff(pars[3:4]), diff(pars[5:6]))}
  eb <- c(1, 0, 0)
  LC_BC_GC(d = d, gof = gof, eq = eq, eb = eb)
}

LC_BF_GC <- function(d, gof) {
  eq <- function(pars, d, gof) {c(sum(pars[3:4]), diff(pars[1:2]), diff(pars[5:6]))}
  eb <- c(1, 0, 0)
  LC_BC_GC(d = d, gof = gof, eq = eq, eb = eb)
}

LC_BC_GF <- function(d, gof) {
  eq <- function(pars, d, gof) {c(sum(pars[3:4]), diff(pars[1:2]), diff(pars[3:4]))}
  eb <- c(1, 0, 0)
  LC_BC_GC(d = d, gof = gof, eq = eq, eb = eb)
}

model_list <- list(
  lc_bc_gc = LC_BC_GC,
  lf_bc_gc = LF_BC_GC,
  lc_bf_gc = LC_BF_GC,
  lc_bc_gf = LC_BC_GF
)

# ==========================================================================
# Fits models
# ==========================================================================
res <- as.data.table(expand.grid(sd_l = c(.01, .05, .1, .5, 1), distr = c("norm"), rep = 151:200)) # distr = c("norm", "unif")
res[, id := as.character(1:.N)]
# res <- res[rep %between% c(201, 220)]

x1 <- Sys.time()
if (parallel == TRUE) {
  # registerDoMC(cores = detectCores())
  registerDoFuture()
  plan(multisession, workers = 4L)  ## on MS Windows
  res <- foreach(x = unique(res$id),
                  .combine = "rbind",
                  .inorder = FALSE, 
                  .packages = c("data.table", "purrr", "Rsolnp"),
                  .export = c("model_list", "fits", "samples", "LC_BC_GC", "LF_BC_GC", "LC_BF_GC", "LC_BC_GF")) %dorng% {
                    res[id == x, .(
                      sd_l, distr, rep,
                      model = names(model_list),
                      fit = map(model_list, exec, 
                                d = samples(n, m_l, sd_l, m_h, sd_h, distr, seed = rep), 
                                gof = "ll")), by = id]
                  }   
} else {
  res <- res[, .(
    sd_l, distr, rep,
    model = names(model_list),
    fit = map(model_list, exec, 
              d = samples(n, m_l, sd_l, m_h, sd_h, distr, seed = rep), 
              gof = "ll")), by = id]
}
x2 <- Sys.time()
x2 - x1

# saveRDS(res, "~/Projects/review-geometric/category-variability-simulation-new-c5-2.rds")
res[, fit[[1]][[1]], by = .(id, sd_l, distr, rep, model)]

pars <- res[, as.list(fit[[1]][[1]]), by = .(id, sd_l, distr, rep, model)]
pars <- melt(pars, measure.vars = c("l1", "l0", "b1", "b0", "g1", "g0"), variable.name = "par")
pars <- dcast(pars[, mean(value), by = .(sd_l, distr, model, par)], sd_l + model ~ par + distr)

gofs <- res[, .(gof = fit[[1]][[2]], acc = fit[[1]][[3]]), by = .(id, sd_l, distr, rep, model)]
gofs <- melt(gofs, measure.vars = c("gof", "acc"), variable.name = "var")
gofs <- dcast(gofs[, .(mean(value)), by = .(sd_l, distr, model, var)], sd_l + model ~ var + distr)
gofs <- merge(pars, gofs)

gofs <- gofs[, round(.SD, 2), by = .(sd_l, model)]
gofs <- gofs[order(model, sd_l)]

gofs[, !grepl("unif", colnames(gofs)), with = FALSE][, c(1, 5:8, 3, 4, 9, 10)]
xtable::xtable(gofs[, !grepl("unif", colnames(gofs)), with = FALSE][, c(1, 6, 5, 8, 7, 4, 3, 9, 10)], type = "latex")
xtable::xtable(gofs[, !grepl("unif|model", colnames(gofs)), with = FALSE], type = "latex")















microbenchmark(lf_bc = LF_BC(train = samples(n, m_l, sd_l, m_h, sd_h, distr = "norm", seed = sample(101:999, 1)), 
                             trans = samples(n, m_l, sd_l, m_h, sd_h, distr = "norm", seed = sample(101:999, 1)), 
                             gof, eq = function(pars, train, trans, gof) {c(sum(pars[3:4]), diff(pars[3:4]))}, eb = c(1, 0)),
               lf_bf = LF_BF(train = samples(n, m_l, sd_l, m_h, sd_h, distr = "norm", seed = sample(101:999, 1)), 
                             trans = samples(n, m_l, sd_l, m_h, sd_h, distr = "norm", seed = sample(101:999, 1)), 
                             gof),
               lc_bc = LC_BC(train = samples(n, m_l, sd_l, m_h, sd_h, distr = "norm", seed = sample(101:999, 1)), 
                             trans = samples(n, m_l, sd_l, m_h, sd_h, distr = "norm", seed = sample(101:999, 1)), 
                             gof),
               lc_bf = LC_BF(train = samples(n, m_l, sd_l, m_h, sd_h, distr = "norm", seed = sample(101:999, 1)), 
                             trans = samples(n, m_l, sd_l, m_h, sd_h, distr = "norm", seed = sample(101:999, 1)), 
                             gof), 
               times = 4)


























