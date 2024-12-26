# ==========================================================================
# Reanalysis Wills Experiment 1
# ==========================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, dplyr, ggplot2, tidyr, yarrr, purrr, ICSNP, Rfast)
theme_set(theme_bw())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # change this to the folder where this script is lying; works only in RStudio

dt <- fread("../../../data/processed/wills_ex1data.txt") # resp: 1 = top, 2 = left, 3 = right; triad: see exe1code.txt
stim <- fread("../../../data/processed/wills_ex1code.txt") # each triplet has six triads (arrangements)
# Strategies Wills: unidimensional responses (hull, sail), overall similarity (os), identical attribute (ident)

# ==========================================================================
# Extract feature values
# ==========================================================================
triads <- stim[, c("triad", "top1", "left2", "right3"), with = FALSE]
triads <- as.data.table(pivot_longer(triads, !triad, names_to = "stim"))
triads[, c("f1", "f2") := list(f1 = floor(value / 10), f2 = value %% 10)]

# ==========================================================================
# Extracts feature values
# ==========================================================================
diffs <- triads[, .(pair = c("top-left", "left-right", "top-right"), 
                    f1 = abs(c(diff(f1), diff(f1, lag = 2))), 
                    f2 = abs(c(diff(f2), diff(f2, lag = 2)))), by = triad]

# ==========================================================================
# Computes pairwise similarities and response probabilities
# ==========================================================================
similarity <- function(x, w, c, p, r) {
  exp(-c * (x^r %*% w)^(p/r))
}

x <- as.matrix(diffs[, .(f1, f2)])
c <- 1
diffs[, ":=" (disc = similarity(x != 0, w = c(.5, .5), c = c, p = 1, r = 1),
              cb = similarity(x, w = c(.5, .5), c = c, p = 1, r = 1),
              eucl = similarity(x, w = c(.5, .5), c = c, p = 1, r = 2),
              cb_f1 = similarity(x, w = c(1, 0), c = c, p = 1, r = 1),
              cb_f2 = similarity(x, w = c(0, 1), c = c, p = 1, r = 1),
              disc_f1 = similarity(x != 0, w = c(1, 0), c = c, p = 1, r = 1),
              disc_f2 = similarity(x != 0, w = c(0, 1), c = c, p = 1, r = 1))]
diffs[, id := ifelse(triad %in% diffs[f1 == 0, triad], cb_f1, cb_f2)]

preds <- diffs[, lapply(.SD, function(x) x / sum(x)), by = triad, .SDcols = -c("pair", "f1", "f2")]
preds[, stim := rep(c("right", "top", "left"), length.out = .N)]
preds[, resp := as.numeric(factor(stim, levels = c("top", "left", "right"), labels = 1:3))]

# ==========================================================================
# Computes goodness-of-fit measures
# ==========================================================================
dt <- merge(dt, preds, by = c("triad", "resp"), sort = FALSE)
models <- c("disc", "cb", "eucl", "cb_f1", "cb_f2", "disc_f1", "disc_f2", "id")

gofs_ind <- dt[, lapply(.SD, function(x) -sum(log(x))), by = .(subj, cond), .SDcols = models]
gofs_ind[, ":=" (cb_fx = ifelse(cb_f1 < cb_f2, cb_f1, cb_f2), # optimal unidim model per participant
                 disc_fx = ifelse(disc_f1 < disc_f2, disc_f1, disc_f2))]
gofs_ind <- merge(gofs_ind, dt[, .N, by = cond], by = "cond")
models <- c(models, "cb_fx", "disc_fx")
gofs_ind <- gofs_ind[, lapply(.SD, function(x) x / N), by = .(subj, cond), .SDcols = models]

# ==========================================================================
# Adds prediction of best fitting unidimensional model per participant
# ==========================================================================
dt <- merge(dt, gofs_ind[, .(cb_fx = ifelse(cb_f1 < cb_f2, "cb_f1", "cb_f2"),
                             disc_fx = ifelse(disc_f1 < disc_f2, "disc_f1", "disc_f2")), by = .(subj, cond)])
dt[, ":=" (cb_fx = as.numeric(apply(dt, 1, function(x) x[x["cb_fx"]])),
           disc_fx = as.numeric(apply(dt, 1, function(x) x[x["disc_fx"]])))]

models <- c("cb", "id", "cb_fx", "disc")
names(models) <- c("MULTI-MINK", "UNI-MINK (identity)", "UNI-MINK (unidimensional)", "MULTI-DISC")

# ==========================================================================
# Computes binary goodness-of-fit measures (chosen vs. not chosen)
# ==========================================================================
gofs_ind_bin <- dt[, lapply(.SD, function(x) sum(log(x))), by = .(subj, cond), .SDcols = models]
gofs_ind_bin <- merge(gofs_ind_bin, gofs_ind_bin[, .(model = names(.SD)[which.max(.SD)]), .SDcols = models, by = .(cond, subj)])
melt(gofs_ind_bin, id.vars = c("subj", "cond", "model"))[, sort(value, decreasing = T)[2] / max(value), by = .(subj, cond, model)][, mean(V1), by = .(cond, model)]
weights <- gofs_ind_bin[, exp(.SD - max(.SD)), .SDcols = models, by = .(subj, cond, model)]
weights[, (models) := .SD/ rowSums(.SD), .SDcols = models] # equal to AIC weights
melt(weights, id.vars = c("subj", "cond", "model"))[model == variable, mean(value), by = .(cond, model)]

gofs_distr <- gofs_ind_bin[, .(n = .N), by = list(cond, model)]
gofs_distr <- merge(gofs_distr, gofs_ind_bin[, .(N_cond = .N), by = list(cond)], by = "cond")
gofs_distr <- merge(gofs_distr, gofs_ind_bin[, .(N_model = .N), by = list(model)], by = "model")
gofs_distr[, ":=" (Model = factor(model, levels = models, labels = c(names(models))),
                   Condition = factor(cond, levels = rev(unique(cond)), labels = rev(unique(cond))),
                   our_model = (model == "disc"))]

plots_gofs <- function(dt, x, y, fill, color = NULL, pos = "stack", facet = NULL, limits = NULL, add = c(0, 0)) {
  p <- ggplot(data = dt, aes_string(x = x, y = y, fill = fill, color = color)) + 
    geom_bar(stat = "identity", position = pos, size = 2) +
    scale_fill_manual(values = c(grey(.9), grey(0.6), grey(0.3), grey(0))) +
    # scale_fill_manual(values = c(grey(.9), grey(0.675), grey(0.45), grey(0.225), grey(0))) +
    scale_color_manual(values = c("FALSE" = alpha("black", 0), "TRUE" = grey(0))) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(limits = limits, expand = expansion(add = add)) +
    labs(y = paste0(ifelse(grepl("/", y), "Share", "Number"), " of Described Participants"),
         x = "Time Pressure (in ms)") + 
    theme(axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = "bottom") +
    guides(color = "none")
  if(!is.null(facet)) p <- p + facet_wrap(~get(facet), nrow = 1) + theme(axis.text.x = element_blank())
  p
}

plots_gofs(dt = gofs_distr, x = "Condition", y = "n/N_cond", fill = "Model")
# ggsave("../../figures/figure2.png", width = 6.4, height = 3.2)

# ==========================================================================
# Analyzes goodness-of-fit measures
# ==========================================================================
weights <- gofs_ind_bin[, exp(.SD - max(.SD)), by = list(subj, cond, model), .SDcols = models]
weights[, (models) := round(.SD/ rowSums(.SD), 4), by = list(subj, cond, model)]

weights_above_crit <- weights[, base::max(.SD) >= .70, by = list(cond, subj, model)][, V1]
weights[weights_above_crit, .N, by = list(cond, model)][order(cond)]
