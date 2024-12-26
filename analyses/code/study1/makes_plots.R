# ==========================================================================
# Makes evidence strength (AIC weights) figure
# ==========================================================================

# ==========================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # change wd to where this script lies (RStudio needed!)
source("3-modeling-analysis.R")
models <- c("gcm", "gcm_disc", "gcm_unidim", "gcm_disc_unidim", "random")

pacman::p_load(data.table, ggplot2, ggpattern, ggnewscale)
theme_set(theme_bw())
# ==========================================================================

# ==========================================================================
# Makes evidence strength plot (Figure 6)
# ==========================================================================
weights <- gofs[, exp(.SD - max(.SD)), .SDcols = models, by = .(subj, cond)]
weights[, (models) := .SD/ rowSums(.SD), .SDcols = models] # equal to AIC weights
weights[, best_model_exists := rowSums(.SD > assign) > 0, .SDcols = models]
weights[, mds := cmdscale(dist(.SD), k = 1), .SDcols = models, by = .(cond, best_model_exists)]
setorderv(weights, cols = c("cond", "best_model_exists", "mds"), order = c(-1, -1, 1))
weights[, ord := 1:.N, by = .(cond)]
weights[best_model_exists == FALSE, ord := ord + 2 + as.numeric(cond)]
weights[cond == "no_tp" & ord %between% c(1, 30), ord := c(9:30, 1:8)]
weights[cond == "tp" & ord %between% c(15, 27), ord := 27:15]

weights <- melt(weights, id.vars = c("subj", "cond", "ord", "mds", "best_model_exists"), variable.name = "model", value.name = "weight")
weights[, model := factor(model, levels = c("gcm", "gcm_unidim", "gcm_disc", "gcm_disc_unidim", "random"), labels = c("MULTI-MINK", "UNI-MINK", "MULTI-DISC", "UNI-DISC", "RANDOM"))]
weights[, cond := factor(cond, levels = c("no_tp", "tp"), labels = c("No time pressure", "Time pressure"))]
dt_text <- weights[, .(f1 = mean(ord), f2 = 1.05), by = .(cond, best_model_exists)]
dt_text[, text := c("classified as belonging to a model", "not cl.", "classified as belonging to a model", "not cl.")]

ggplot(weights, aes(x = ord, y = weight)) +
  geom_bar_pattern(aes(pattern = model, fill = model),
                   stat = "identity", position = "fill",
                   color = "black", 
                   pattern_fill = "white",
                   pattern_angle = 45,
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_key_scale_factor = 0.6) +
  scale_pattern_manual(values = c("none", "stripe", "none", "stripe", "none"), name = "Model") +
  scale_fill_manual(values = c(grey(.4), grey(.4), grey(.7), grey(.7), grey(1)), name = "Model") +
  facet_wrap(~cond, nrow = 2) +
  xlab("Participants") +
  ylab("Evidence strength") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank()) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 34)) +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1), expand = c(0, 0), limits = c(0, 1.08)) +
  geom_text(data = dt_text, mapping = aes(x = f1, y = f2, label = text), size = 2.2, hjust = .5, lineheight = .7)

ggsave("../../figures/figure6.png", height = 3.5, width = 7)

# ==========================================================================
# Makes transfer phase plot (Figure 5)
# ==========================================================================
d[, type := ifelse(stimulus == "100-NA", "100", ifelse(stimulus == "003-NA", "003", "221, 231, 321, 331"))]
d_sub <- melt(d, id.vars = c("subj", "cond", "stimulus", "type"), measure.vars = c(models[grepl("^gcm", models)]), variable.name = "model", value.name = "preds")
d_sub <- merge(d_sub, pars[, .(subj, model, paste0(model, "_", w1, w2, w3))])
d_sub[!grepl("unidim", model), V3 := model]
d_sub <- d_sub[model != "gcm_disc_unidim"]
d_sub[, model := NULL]
colnames(d_sub)[colnames(d_sub) == "V3"] <- "model"
preds <- d_sub[, .(m = mean(preds), sd = sd(preds)), by = .(model, type, subj)]
preds <- preds[, .(m = mean(m), sd = mean(sd), se = sd(m)/sqrt(length(unique(subj)))), by = .(model, type)]
preds[, model := factor(model, levels = c("gcm", "gcm_unidim_100", "gcm_unidim_001", "gcm_disc"), labels = c("MULTI-MINK", "UNI-MINK2", "UNI-MINK", "MULTI-DISC"))]
preds[, type := factor(type, levels = c("003", "221, 231, 321, 331", "100"))]
obs <- d[, .(m = mean(resp, na.rm = T), sd = sd(resp, na.rm = T)), by = .(time_pressure = cond, type, subj)]
obs <- obs[, .(m = mean(m), sd = mean(sd), se = sd(m)/sqrt(length(unique(subj)))), by = .(time_pressure, type)]
obs[, time_pressure := factor(time_pressure, levels = c("no_tp", "tp"), labels = c("No time pressure", "Time pressure"))]
obs[, type := factor(type, levels = c("003", "221, 231, 321, 331", "100"))]

dt_plot <- rbind(preds, obs, fill = TRUE)
dt_plot[, x := model]
dt_plot[is.na(model), x := time_pressure]
dt_plot[, x := factor(x, levels = c("No time pressure", "Time pressure", "MULTI-MINK", "UNI-MINK", "UNI-MINK2", "MULTI-DISC"))]
dt_plot[, var := ifelse(is.na(model), se, 0)]

ggplot(dt_plot[!x == "UNI-MINK2"], aes(type)) +
  geom_errorbar(aes(ymin = m - var, ymax = m + var, group = x), color = "black", position = position_dodge(.8), width = .3) +
  geom_point(aes(y = m, fill = time_pressure, alpha = x, group = x), shape = 21, size = 4, position = position_dodge(.8)) +
  scale_fill_manual(values = c("black", "white"), name = "Condition", breaks = c("No time pressure", "Time pressure")) +
  scale_alpha_manual(values = c(rep(1, 2), rep(0, 3)), guide = "none") +
  new_scale("fill") +
  new_scale("alpha") +
  geom_bar(aes(y = m * as.numeric(!is.na(model)), fill = model, alpha = x, group = x), color = "black", stat = "identity", position = "dodge", width = .8) +
  scale_fill_manual(values = c(grey(.5), grey(.7), grey(.9)), name = "Model", breaks = c("MULTI-MINK", "UNI-MINK", "MULTI-DISC")) +
  scale_alpha_manual(values = c(rep(0, 2), rep(1, 3)), guide = "none") +
  labs(x = "Stimulus", y = "Classification P(A|Stimulus)", shape = "Model") +
  theme(axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(), 
        axis.line.y = element_line()) +
  scale_y_continuous(breaks = c(0, .5, 1), expand = c(0, 0), limits = c(0, 1))

ggsave("../../figures/figure5.png", width = 7, height = 3)
