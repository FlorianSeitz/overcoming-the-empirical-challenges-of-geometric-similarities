# ==========================================================================
# Makes evidence strength (AIC weights) figure
# ==========================================================================

# ==========================================================================
# Parameters that need to be specified
models <- c("gcm", "gcm_disc", "gcm_unidim", "gcm_disc_unidim", "random") # specify the models to analyze
assign <- .70 # model assignment threshold
# ==========================================================================

# ==========================================================================
# Prepares packages and data
# ==========================================================================
pacman::p_load(data.table, ggplot2, ggpattern, ggnewscale)
theme_set(theme_bw())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # change wd to where this script lies (RStudio needed!)
res <- fread("../../other/similarity_pars_likelihoods_cond50.csv")

# ==========================================================================
# Makes evidence strength plot (Figure 9)
# ==========================================================================
gofs <- res[, .(ll = median(ll_test)), by = .(subj, time_pressure, model)]
gofs <- dcast(gofs, subj + time_pressure ~ model, value.var = "ll")
gofs[, random := 0] # because density for a uniform model is 1

gofs[, best_model := names(which.max(.SD)), by = .(subj, time_pressure)]
weights <- gofs[, exp(.SD - max(.SD)), .SDcols = models, by = .(subj, time_pressure)]
weights[, (models) := .SD/ rowSums(.SD), .SDcols = models] # equal to AIC weights

weights[, best_model_exists := rowSums(.SD > assign) > 0, .SDcols = models]
weights[, mds := cmdscale(dist(.SD), k = 1), .SDcols = models, by = .(time_pressure, best_model_exists)]
setorderv(weights, cols = c("time_pressure", "best_model_exists", "mds"), order = c(-1, -1, 1))
weights[, ord := 1:.N, by = .(time_pressure)]
weights[best_model_exists == FALSE, ord := ord + 2]
weights[time_pressure == "none" & ord %between% c(42, 58), ord := 58:42]
weights[time_pressure == "50" & ord %between% c(50, 51), ord := 51:50]
weights[time_pressure == "50" & ord %between% c(51, 61), ord := 61:51]

weights <- melt(weights, id.vars = c("subj", "time_pressure", "ord", "mds", "best_model_exists"), variable.name = "model", value.name = "weight")
weights[, model := factor(model, levels = c("gcm", "gcm_unidim", "gcm_disc", "gcm_disc_unidim", "random"), labels = c("MULTI-MINK", "UNI-MINK", "MULTI-DISC", "UNI-DISC", "RANDOM"))]
weights[, time_pressure := factor(time_pressure, levels = c("none", "50"), labels = c("No time pressure", "Time pressure"))]
dt_text <- data.table(text = c("classified as belonging to a model", "not cl.", "classified as belonging to a model", "not cl."), 
                      time_pressure = c("No time pressure", "No time pressure", "Time pressure", "Time pressure"),
                      f1 = c(29.5, 63.5, 31, 65), f2 = 1.05)

ggplot(weights, aes(x = ord, y = weight)) +
  # geom_bar(aes(fill = model), position = "fill", stat = "identity", color = "black", size = 0.25) +
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
  # scale_fill_manual(values = c("#424242", "#424242", "#BDBDBD", "#BDBDBD", "white"), name = "Model") +
  facet_wrap(~time_pressure, nrow = 2) +
  xlab("Participants") +
  ylab("Evidence strength") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank()) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 67)) +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1), expand = c(0, 0), limits = c(0, 1.08)) +
  geom_text(data = dt_text, mapping = aes(x = f1, y = f2, label = text), size = 2.2, hjust = .5, lineheight = .7)

ggsave("../../figures/figure9.png", height = 3.5, width = 7)

# ==========================================================================
# Makes test phase plot (Figure 8)
# ==========================================================================
preds <- fread("../../other/similarity_predictions_cond50.csv")
preds <- preds[, .(m = mean(m), sd = mean(sd), se = sd(m)/sqrt(length(unique(subj)))), by = .(model, type)]
preds[, model := factor(model, levels = c("gcm", "gcm_unidim", "gcm_disc", "gcm_disc_unidim"), labels = c("MULTI-MINK", "UNI-MINK", "MULTI-DISC", "UNI-DISC"))]
preds[, type := factor(type, levels = c("two-large", "one-large", "all-small", "two-small"), labels = c("two-large", "one-large", "all-small", "two-small"))]

obs <- fread("../../../data/processed/similarity.csv")[cond == "tp_50"][grep("small|large", type)]
obs <- obs[, .(m = mean(resp, na.rm = T), sd = sd(resp, na.rm = T)), by = .(time_pressure = phase, type, subj)]
obs <- obs[, .(m = mean(m), sd = mean(sd), se = sd(m)/sqrt(length(unique(subj)))), by = .(time_pressure, type)]
obs[, time_pressure := factor(time_pressure, levels = c("no_tp", "tp"), labels = c("No time pressure", "Time pressure"))]
obs[, type := factor(type, levels = c("two-large", "one-large", "all-small", "two-small"), labels = c("two-large", "one-large", "all-small", "two-small"))]

dt_plot <- rbind(preds, obs, fill = TRUE)
dt_plot[, x := model]
dt_plot[is.na(model), x := time_pressure]
dt_plot[, x := factor(x, levels = c("No time pressure", "Time pressure", "MULTI-MINK", "UNI-MINK", "MULTI-DISC", "UNI-DISC"))]
dt_plot[, var := ifelse(is.na(model), sd, 0)]

ggplot(dt_plot[!x == "UNI-DISC"], aes(type)) +
  geom_errorbar(aes(ymin = m - var, ymax = m + var, group = x), color = "black", position = position_dodge(.8), width = .3) +
  geom_point(aes(y = m, fill = time_pressure, alpha = x, group = x), shape = 21, size = 4, position = position_dodge(.8)) +
  scale_fill_manual(values = c("black", "white"), name = "Condition", breaks = c("No time pressure", "Time pressure")) +
  scale_alpha_manual(values = c(rep(1, 2), rep(0, 3)), guide = "none") +
  new_scale("fill") +
  new_scale("alpha") +
  geom_bar(aes(y = m * as.numeric(!is.na(model)), fill = model, alpha = x, group = x), color = "black", stat = "identity", position = "dodge", width = .8) +
  scale_fill_manual(values = c(grey(.5), grey(.7), grey(.9)), name = "Model", breaks = c("MULTI-MINK", "UNI-MINK", "MULTI-DISC")) +
  scale_alpha_manual(values = c(rep(0, 2), rep(1, 3)), guide = "none") +
  labs(x = "Stimulus pair type", y = "Similarity", shape = "Model") +
  theme(axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        axis.line.x = element_line(), 
        axis.line.y = element_line()) +
  # scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, .5, 1), expand = c(0, 0), limits = c(0, 1))

ggsave("../../figures/figure8.png", width = 7, height = 3)

