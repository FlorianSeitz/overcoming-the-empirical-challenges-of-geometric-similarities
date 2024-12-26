library(data.table)
library(combinat)

source("Models/GCM_Prot.R", chdir = TRUE)
source("VariableDefinition.R", chdir = TRUE)

learningset <- do.call(make.learningset, list(a.design.ID = 121, n.Cat1.IDs, n.dim, n.val, n.stim, n.transfer, n.learning.stim, n.Cat1, all.stim))$learningset
learningset <- learningset[, c(4, 1, 2, 3, 6)]
learningset[, 2:4] <- learningset[, 2:4] + 1

learn_test <- learningset[learningset$ID %in% c("002", "012", "101", "111"), 1:4]

learningset$ID <- apply(learningset[, 2:4], 1, paste, collapse = "")
learningset$crit <- -1
#learningset$Cat[learningset$Cat == 1] <- "R"
#learningset$Cat[learningset$Cat == 0] <- "L"

fwrite(learningset, "C:/Users/Sylvia/Documents/Masterarbeit/Experiment_Python/Categorization Experiment/exemplars_training.csv", sep = ";", row.names = FALSE)

remaining <- as.data.table(do.call(make.learningset, list(a.design.ID = 121, n.Cat1.IDs, n.dim, n.val, n.stim, n.transfer, n.learning.stim, n.Cat1, all.stim))$testset)
test.ids <- c("003", "100", "231", "321", "221", "331")
testset <- remaining[remaining$ID %in% test.ids, ]
testset <- testset[, c(4, 1, 2, 3)]

testset[, 2:4] <- testset[, 2:4] + 1
testset <- rbind(testset, learn_test)
testset$ID <- apply(testset[, 2:4], 1, paste, collapse = "")
testset$crit <- -1

fwrite(testset, "C:/Users/Sylvia/Documents/Masterarbeit/Experiment_Python/Categorization Experiment/exemplars_test.csv", sep = ";", row.names = FALSE)

# remaining <- remaining[remaining$ID %in% test.ids == FALSE, ]
# remaining <- remaining[, c(4, 1, 2, 3)]
# remaining[, 2:4] <- remaining[, 2:4] + 1
# remaining$ID <- apply(remaining[, 2:4], 1, paste, collapse = "")
# 
# res <- NULL
# for(j in 1:3){
#   ws <- c(0, 0, 0)
#   ws[j] <- 1
#   for(i in 1:nrow(remaining)) {
#     values <- unlist(remaining[i, 2:4])
#     disc <- Predict(values, learningset, Metric = "Attr", Model = "GCM", c = 1, w = c(0.33, 0.33, 0.34), p = 1, r = 1)
#     mink <- Predict(values, learningset, Metric = "Eucl", Model = "GCM", c = 1, w = ws, p = 1, r = 1)
#     res <- rbind(res, data.table(remaining$ID[i], disc, mink, j))
#   }
# }
# 
# res[, pred.diff := abs(disc-mink)]
# res[, order(pred.diff), by = j]
# filler.ids <- c(res[j == 1 & pred.diff >= 0.16, V1], res[j == 2 & pred.diff >= 0.08, V1], res[j == 3 & pred.diff >= 0.16, V1])
# 
# # Filler items
# filler <- remaining[remaining$ID %in% filler.ids, ]
# filler$crit <- -1
# 
# fwrite(filler, "C:/Users/Sylvia/Documents/Masterarbeit/Experiment_Python/Categorization Experiment/exemplars_filler.csv", sep = ";", row.names = FALSE)

# Plots: UD(X) vs. DISC-EQ vs. DISC-UD(X) --> Learningset diskriminiert für UD(2) und DISC-EQ