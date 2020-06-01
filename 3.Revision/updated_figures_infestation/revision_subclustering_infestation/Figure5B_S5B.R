library(reshape2)
library(scales)
library(ggplot2)
library(plyr)
library(RColorBrewer)
label <- readRDS('tmp/infested_label.Rds')
head(label)
#dir.create('stats')


### subclustering ###
label_summ <- label[, c(6, 7)]
label_summ <- droplevels(label_summ)
head(label_summ)

result <- data.frame(matrix(nrow = 2, ncol = 12))
colnames(result) <- c('timepoint', levels(label_summ$labelTransfer))
result$timepoint <- levels(label_summ$timepoint)

for (new in levels(label_summ$labelTransfer)) {
  tmp <- subset(label_summ, labelTransfer == new)
  result[, new] <- summary(tmp$timepoint)
}

for (new in levels(label_summ$labelTransfer)) {
  result[, new] <- round(result[, new]/sum(result[, new])*100, digits = 2)
} 

result_m <- melt(result)
result_m$timepoint <- factor(result_m$timepoint, levels = levels(label_summ$timepoint))
head(result_m)


p <- ggplot(result_m ,aes(x = variable, y = value, fill = timepoint)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('dodgerblue', 'red2')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
ggsave('stats/subclustering.proportions.pdf', units = 'cm', width = 10, height = 10)



### cell types ###
label_summ <- label[, c(6, 8)]
label_summ <- droplevels(label_summ)
head(label_summ)

result <- data.frame(matrix(nrow = 2, ncol = 7))
colnames(result) <- c('timepoint', levels(label_summ$labelTransfer_simple))
result$timepoint <- levels(label_summ$timepoint)

for (new in levels(label_summ$labelTransfer_simple)) {
  tmp <- subset(label_summ, labelTransfer_simple == new)
  result[, new] <- summary(tmp$timepoint)
}

for (new in levels(label_summ$labelTransfer_simple)) {
  result[, new] <- round(result[, new]/sum(result[, new])*100, digits = 2)
} 

result_m <- melt(result)
result_m$timepoint <- factor(result_m$timepoint, levels = levels(label_summ$timepoint))
head(result_m)


p <- ggplot(result_m ,aes(x = variable, y = value, fill = timepoint)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('dodgerblue', 'red2')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
ggsave('stats/celltypes.proportions.pdf', units = 'cm', width = 10, height = 10)



########################
### Normalized count ###
########################
real_normal_med <- 1328
sc_normal_total <- rowSums(result[, -c(1)])[1]
normal_scaleFactor <- sc_normal_total/real_normal_med

real_infested_med <- 1786.5
sc_infested_total <- rowSums(result[, -c(1)])[2]
infested_scaleFactor <- sc_infested_total/real_infested_med
########################
########################
########################

### subclustering ###
label_summ <- label[, c(6, 7)]
label_summ <- droplevels(label_summ)
head(label_summ)

result <- data.frame(matrix(nrow = 2, ncol = 12))
colnames(result) <- c('timepoint', levels(label_summ$labelTransfer))
result$timepoint <- levels(label_summ$timepoint)

for (new in levels(label_summ$labelTransfer)) {
  tmp <- subset(label_summ, labelTransfer == new)
  result[, new] <- c(summary(tmp$timepoint)[1]/normal_scaleFactor, summary(tmp$timepoint)[2]/infested_scaleFactor)
}
result

result_m <- melt(result)
result_m$timepoint <- factor(result_m$timepoint, levels = levels(label_summ$timepoint))
head(result_m)
p <- ggplot(result_m ,aes(x = variable, y = value, fill = timepoint)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('dodgerblue', 'red2')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
ggsave('stats/subclustering.normalizedCount.pdf', units = 'cm', width = 10, height = 10)

p <- ggplot(subset(result_m, variable == 'PH 4' | variable == 'PM 1') ,aes(x = variable, y = value, fill = timepoint)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('dodgerblue', 'red2')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
ggsave('stats/subclustering.normalizedCount_pt1.pdf', units = 'cm', width = 7, height = 10)

p <- ggplot(subset(result_m, variable != 'PH 4' & variable != 'PM 1') ,aes(x = variable, y = value, fill = timepoint)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('dodgerblue', 'red2')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
ggsave('stats/subclustering.normalizedCount_pt2.pdf', units = 'cm', width = 9, height = 10)



### cell types ###
label_summ <- label[, c(6, 8)]
label_summ <- droplevels(label_summ)
head(label_summ)

result <- data.frame(matrix(nrow = 2, ncol = 7))
colnames(result) <- c('timepoint', levels(label_summ$labelTransfer_simple))
result$timepoint <- levels(label_summ$timepoint)

for (new in levels(label_summ$labelTransfer_simple)) {
  tmp <- subset(label_summ, labelTransfer_simple == new)
  result[, new] <- c(summary(tmp$timepoint)[1]/normal_scaleFactor, summary(tmp$timepoint)[2]/infested_scaleFactor)
}
result

result_m <- melt(result)
result_m$timepoint <- factor(result_m$timepoint, levels = levels(label_summ$timepoint))
head(result_m)
p <- ggplot(result_m ,aes(x = variable, y = value, fill = timepoint)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('dodgerblue', 'red2')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
ggsave('stats/celltypes.normalizedCount.pdf', units = 'cm', width = 10, height = 10)

p <- ggplot(subset(result_m, variable == 'PH' | variable == 'PM') ,aes(x = variable, y = value, fill = timepoint)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('dodgerblue', 'red2')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
ggsave('stats/celltypes.normalizedCount_pt1.pdf', units = 'cm', width = 7, height = 10)

p <- ggplot(subset(result_m, variable != 'PH' & variable != 'PM') ,aes(x = variable, y = value, fill = timepoint)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('dodgerblue', 'red2')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, colour = 'black'),
        axis.ticks = element_line(colour = 'black')) +
  labs(title = '', x = '', y = ''); p
ggsave('stats/celltypes.normalizedCount_pt2.pdf', units = 'cm', width = 9, height = 10)








