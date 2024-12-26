rm(list = ls(all = T)) # Empty workspace
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, magic)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # change wd to where this script lies (RStudio needed!)

n.val = 5
stimuli <- expand.grid(Dim1 = 1:n.val, Dim2 = 1:n.val, Dim3 = 1:n.val)
stimuli$ID <- apply(stimuli, 1, paste0, collapse = "")
stimuli$crit <- -1
stimuli <- stimuli[, c(4, 1, 2, 3, 5)]

latin_square <- rlatin(1, size = 3)
id_1 <- vector()
id_2 <- vector()
type <- vector()

for(i in 1:3){
  # Type A
  type <- c(type, "A")
  stim <- c(sample(c(1, 5), 2), sample(2:4, 1))
  id_1 <- c(id_1, paste0(stim[latin_square[i, ]], collapse = ""))
  stim[1] <- 6 - stim[1]
  id_2 <- c(id_2, paste0(stim[latin_square[i, ]], collapse = ""))

  # Type B
  type <- c(type, "B")
  n_pos_changes <- sample(1:2, 1)
  dim_changes <- c(rep(1, n_pos_changes), rep(-1, 3-n_pos_changes))
  dim_changes <- dim_changes[sample(1:3, 3)]
  stim <- sample(2:4, size = 3, replace = TRUE)
  id_1 <- c(id_1, paste0(stim, collapse = ""))
  id_2 <- c(id_2, paste0(stim + dim_changes, collapse = ""))
  
  # Type C
  type <- c(type, "C")
  stim <- c(sample(c(1, 5), 2), sample(2:4, 1))
  print(stim)
  id_1 <- c(id_1, paste0(stim[latin_square[i, ]], collapse = ""))
  stim[1:2] <- 6 - stim[1:2]
  id_2 <- c(id_2, paste0(stim[latin_square[i, ]], collapse = ""))

  # Type D
  type <- c(type, "D")
  dim_changes <- c(sample(c(-1, 1), 2), 0)
  stim <- c(sample(2:4, size = 2, replace = TRUE), sample(1:5, 1))
  id_1 <- c(id_1, paste0(stim[latin_square[i, ]], collapse = ""))
  id_2 <- c(id_2, paste0((stim + dim_changes)[latin_square[i, ]], collapse = ""))
}

# Type I
for(i in 1:4) {
  identical <- paste0(sample(1:5, size = 3, replace = TRUE), collapse = "")
  id_1 <- c(id_1, identical)
  id_2 <- c(id_2, identical)
  type <- c(type, "I")
}

# Type V
id_1 <- c(id_1, rep("111", 3), rep("555", 2))
id_2 <- c(id_2, rep(c("222", "444", "555"), length.out = 5))
type <- c(type, rep("V", 5))

# Makes test stimuli
stim_1 <- NULL
stim_2 <- NULL
for(i in 1:length(id_1)) {
  stim_1 <- rbind(stim_1, stimuli[stimuli$ID == id_1[i], ])
  stim_2 <- rbind(stim_2, stimuli[stimuli$ID == id_2[i], ])
}

stim_2$type = stim_1$type = type

# Checks if all worked out
all(id_1 == stim_1$ID)
all(id_2 == stim_2$ID)

fwrite(stim_1, "stimuli.1.csv", sep = ";", row.names = FALSE, showProgress = TRUE)
fwrite(stim_2, "stimuli.2.csv", sep = ";", row.names = FALSE, showProgress = TRUE)

# Familiarization
fam.id <- stimuli$ID[stimuli$ID %in% c(stim_1$ID, stim_2$ID) == FALSE]
fam.stimuli <- expand.grid(fam.id, fam.id, stringsAsFactors = FALSE)

fam.stimuli.left <- matrix(unlist(strsplit(fam.stimuli[, 1], split = "")), ncol = 3, byrow = T)
fam.stimuli.right <- matrix(unlist(strsplit(fam.stimuli[, 2], split = "")), ncol = 3, byrow = T)
fam.pairs <- rbind(fam.stimuli[sample(which(rowSums(fam.stimuli.left == fam.stimuli.right) == 0), 10), ],
                   fam.stimuli[sample(which(rowSums(fam.stimuli.left == fam.stimuli.right) == 1), 10), ],
                   fam.stimuli[sample(which(rowSums(fam.stimuli.left == fam.stimuli.right) == 2), 10), ],
                   fam.stimuli[sample(which(rowSums(fam.stimuli.left == fam.stimuli.right) == 3), 4), ])

fam.stim.left <- as.data.frame(t(sapply(fam.pairs[, 1], function(x) stimuli[stimuli$ID == x, ])), row.names = FALSE)
fam.stim.right <- as.data.frame(t(sapply(fam.pairs[, 2], function(x) stimuli[stimuli$ID == x, ])), row.names = FALSE)

fam.stim.right$type = fam.stim.left$type = "fam"

fwrite(fam.stim.left, "familiarization.stimuli.left.csv", sep = ";", row.names = FALSE, showProgress = TRUE)
fwrite(fam.stim.right, "familiarization.stimuli.right.csv", sep = ";", row.names = FALSE, showProgress = TRUE)
